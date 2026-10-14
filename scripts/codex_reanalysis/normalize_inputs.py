"""Normalise delivered FASTQ inputs that were not filtered like the cohort.

Two deliverables were rebuilt after corruption and lost the cohort's minimum
read-length floor. Mapping every animal on the same read population requires a
derived, length-filtered FASTQ for those samples only.

The delivery root is immutable, so derivatives are written under the reanalysis
output root and registered with the hash of the delivered source, the hash of
the derivative, and the exact number of reads and bases removed. Samples that
already meet the floor are recorded as ``not_required`` rather than omitted, so
the record states what happened to all thirteen.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance


class NormalizationError(RuntimeError):
    """Raised when an input cannot be normalised into a trustworthy derivative."""


NORMALIZED_FIELDS = (
    "deliverable_id",
    "sample_id",
    "action",
    "source_relative_path",
    "source_sha256",
    "normalized_relative_path",
    "normalized_sha256",
    "normalized_byte_size",
    "minimum_read_length",
    "reads_in",
    "reads_out",
    "reads_removed",
    "bases_in",
    "bases_out",
    "bases_removed",
)


def plan_normalization(
    qc_rows: Sequence[Mapping[str, str]], minimum_read_length: int
) -> list[dict[str, object]]:
    """Decide, for every audited sample, whether a length filter is required."""
    if minimum_read_length < 1:
        raise NormalizationError("minimum read length must be positive")
    plan: list[dict[str, object]] = []
    for row in qc_rows:
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            raise NormalizationError("audit QC row has no sample identifier")
        raw = str(row.get("min_read_length", "")).strip()
        try:
            observed = int(raw)
        except ValueError as error:
            raise NormalizationError(
                f"audit QC row for {sample_id} has an unreadable minimum read length: {raw!r}"
            ) from error
        plan.append(
            {
                "sample_id": sample_id,
                "observed_min_read_length": observed,
                "action": "length_filter" if observed < minimum_read_length else "not_required",
            }
        )
    return plan


def _records(handle) -> Iterable[tuple[str, str, str, str]]:
    record_number = 0
    while True:
        header = handle.readline()
        if header == "":
            return
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        record_number += 1
        if not sequence or not separator or not quality:
            raise NormalizationError(f"truncated FASTQ record {record_number}")
        if not header.startswith("@") or not separator.startswith("+"):
            raise NormalizationError(f"malformed FASTQ record {record_number}")
        if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
            raise NormalizationError(f"quality length mismatch in FASTQ record {record_number}")
        yield header, sequence, separator, quality


def filter_fastq_by_length(
    source: Path,
    destination: Path,
    minimum_read_length: int,
    compressor: str | None = None,
) -> dict[str, int]:
    """Write only reads at or above the floor, keeping kept records verbatim.

    The destination is published only after the whole source has been read, so a
    broken source can never leave a plausible-looking partial derivative.
    """
    if minimum_read_length < 1:
        raise NormalizationError("minimum read length must be positive")
    source_path = Path(source)
    final_path = Path(destination)
    partial = final_path.with_suffix(f"{final_path.suffix}.partial")
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.unlink(missing_ok=True)

    counts = {"reads_in": 0, "reads_out": 0, "bases_in": 0, "bases_out": 0}
    writer_process = None
    try:
        with gzip.open(source_path, "rt", encoding="ascii", errors="strict", newline="") as reader:
            if compressor:
                handle = partial.open("wb")
                writer_process = subprocess.Popen(
                    [compressor, "-c"], stdin=subprocess.PIPE, stdout=handle, text=True
                )
                handle.close()
                writer = writer_process.stdin
            else:
                writer = gzip.open(partial, "wt", encoding="ascii", newline="")
            try:
                for header, sequence, separator, quality in _records(reader):
                    length = len(sequence.rstrip("\r\n"))
                    counts["reads_in"] += 1
                    counts["bases_in"] += length
                    if length >= minimum_read_length:
                        counts["reads_out"] += 1
                        counts["bases_out"] += length
                        writer.write(header)
                        writer.write(sequence)
                        writer.write(separator)
                        writer.write(quality)
            finally:
                writer.close()
        if writer_process is not None and writer_process.wait() != 0:
            raise NormalizationError("external compressor failed")
    except (gzip.BadGzipFile, EOFError, UnicodeDecodeError) as error:
        partial.unlink(missing_ok=True)
        raise NormalizationError("malformed or truncated gzip FASTQ stream") from error
    except BaseException:
        if writer_process is not None and writer_process.poll() is None:
            writer_process.kill()
            writer_process.wait()
        partial.unlink(missing_ok=True)
        raise

    counts["reads_removed"] = counts["reads_in"] - counts["reads_out"]
    counts["bases_removed"] = counts["bases_in"] - counts["bases_out"]
    partial.replace(final_path)
    return counts


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path, label: str) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file():
        raise NormalizationError(f"{label} is absent: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise NormalizationError(f"{label} has no header")
        return [dict(row) for row in reader]


def _verify_source(path: Path, expected_size: int, expected_sha256: str) -> None:
    if not path.is_file():
        raise NormalizationError(f"delivered FASTQ is absent: {path}")
    if path.stat().st_size != expected_size:
        raise NormalizationError(f"delivered FASTQ byte size does not match the audit: {path.name}")
    if sha256_file(path).lower() != expected_sha256.lower():
        raise NormalizationError(f"delivered FASTQ SHA-256 does not match the audit: {path.name}")


def run_normalization(
    *,
    deliver_root: Path,
    audit_dir: Path,
    output_dir: Path,
    minimum_read_length: int,
    compressor: str | None,
) -> dict[str, object]:
    """Normalise every sample that needs it and register the result."""
    audit = Path(audit_dir).resolve()
    output = Path(output_dir).resolve()
    _validate_output_root(output)
    output.mkdir(parents=True, exist_ok=True)

    qc = read_tsv(audit / "fastq_qc.tsv", "audit QC table")
    manifest = {row["sample_id"]: row for row in read_tsv(audit / "input_manifest.tsv", "input manifest")}
    plan = plan_normalization(qc, minimum_read_length)

    rows: list[dict[str, object]] = []
    for entry in plan:
        sample_id = str(entry["sample_id"])
        source_row = manifest.get(sample_id)
        if source_row is None:
            raise NormalizationError(f"input manifest has no row for {sample_id}")
        relative = str(source_row["input_relative_path"])
        source_path = (Path(deliver_root).resolve() / relative).resolve()
        record: dict[str, object] = {
            "deliverable_id": source_row["deliverable_id"],
            "sample_id": sample_id,
            "action": entry["action"],
            "source_relative_path": relative,
            "source_sha256": source_row["sha256"],
            "minimum_read_length": minimum_read_length,
        }
        if entry["action"] == "not_required":
            record.update(
                {
                    "normalized_relative_path": "",
                    "normalized_sha256": "",
                    "normalized_byte_size": "",
                    "reads_in": "",
                    "reads_out": "",
                    "reads_removed": "",
                    "bases_in": "",
                    "bases_out": "",
                    "bases_removed": "",
                }
            )
            rows.append(record)
            continue

        _verify_source(source_path, int(source_row["byte_size"]), str(source_row["sha256"]))
        destination = output / f"{sample_id}.min{minimum_read_length}.fq.gz"
        counts = filter_fastq_by_length(source_path, destination, minimum_read_length, compressor)
        record.update(
            {
                # Recorded relative to this manifest so the register stays valid
                # wherever the output directory is mounted or copied.
                "normalized_relative_path": destination.name,
                "normalized_sha256": sha256_file(destination),
                "normalized_byte_size": destination.stat().st_size,
                **counts,
            }
        )
        rows.append(record)

    manifest_path = output / "normalized_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=NORMALIZED_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    filtered = [row for row in rows if row["action"] == "length_filter"]
    status = {
        "status": "pass",
        "minimum_read_length": minimum_read_length,
        "sample_count": len(rows),
        "filtered_sample_count": len(filtered),
        "filtered_samples": [row["sample_id"] for row in filtered],
        "normalized_manifest_sha256": sha256_file(manifest_path),
    }
    (output / "normalization_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(
        output,
        {
            "command": "normalize_inputs.py --deliver-root [external] --audit-dir [safe] --output [safe]",
            "config_hash": "not_applicable",
            "input_hashes": {str(row["sample_id"]): str(row["source_sha256"]) for row in rows},
            "reference_hashes": {},
        },
    )
    return status


def _project_root(path: Path) -> Path:
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise NormalizationError("normalisation output must be inside the project repository")


def _validate_output_root(path: Path) -> None:
    root = _project_root(path)
    safe = (root / "results" / "codex_reanalysis_v1").resolve()
    if path != safe and safe not in path.parents:
        raise NormalizationError("normalisation output must be inside results/codex_reanalysis_v1")


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deliver-root", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--minimum-read-length",
        required=True,
        type=int,
        help="cohort read-length floor; samples below it are filtered to match",
    )
    parser.add_argument(
        "--compressor",
        default=None,
        help="external compressor reading stdin and writing stdout, for example pigz",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    compressor = args.compressor
    if compressor and shutil.which(compressor) is None:
        raise NormalizationError(f"compressor is not available: {compressor}")
    status = run_normalization(
        deliver_root=args.deliver_root,
        audit_dir=args.audit_dir,
        output_dir=args.output,
        minimum_read_length=args.minimum_read_length,
        compressor=compressor,
    )
    sys.stderr.write(
        "normalisation {status}; filtered {filtered_sample_count} of {sample_count} samples\n".format(**status)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
