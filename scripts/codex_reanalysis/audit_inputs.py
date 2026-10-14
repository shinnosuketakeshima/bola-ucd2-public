"""Streaming integrity audit for delivered FASTQ inputs.

No sequence, read identifier, delivery-root path, or private source metadata is
written to stdout or to tracked files.  Results are intended for the ignored
``results/codex_reanalysis_v1/00_input_audit`` directory only.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
from typing import Iterator, TextIO

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.config import load_runtime_config
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.sample_registry import (
    SampleRegistryError,
    build_rows_from_metadata,
    build_sample_registry,
    write_tsv,
)


class FastqValidationError(ValueError):
    """Raised when a FASTQ stream is malformed, truncated, or inconsistent."""


@dataclass(frozen=True)
class FastqSummary:
    read_count: int
    base_count: int
    min_read_length: int
    max_read_length: int
    mean_read_length: float
    mean_base_quality: float
    duplicate_read_ids: int
    cross_sample_read_ids: int = 0


class ReadIdTracker:
    """Exact disk-backed detector for duplicate and cross-sample read IDs.

    Read identifiers are kept only as SHA-256 digests in a temporary SQLite
    database, avoiding an unbounded RAM set and avoiding persistence in audit
    deliverables.
    """

    def __init__(self) -> None:
        self._temporary = tempfile.TemporaryDirectory(prefix="bola-read-id-audit-")
        database = Path(self._temporary.name) / "read_ids.sqlite3"
        self._connection = sqlite3.connect(database)
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute(
            "CREATE TABLE read_ids (read_hash TEXT PRIMARY KEY, sample_id TEXT NOT NULL) WITHOUT ROWID"
        )
        self._connection.execute("BEGIN")

    def observe(self, read_id: str, sample_id: str) -> str | None:
        digest = sha256(read_id.encode("utf-8")).hexdigest()
        try:
            self._connection.execute(
                "INSERT INTO read_ids(read_hash, sample_id) VALUES (?, ?)", (digest, sample_id)
            )
            return None
        except sqlite3.IntegrityError:
            existing = self._connection.execute(
                "SELECT sample_id FROM read_ids WHERE read_hash = ?", (digest,)
            ).fetchone()
            return "same_sample" if existing and existing[0] == sample_id else "cross_sample"

    def close(self) -> None:
        self._connection.commit()
        self._connection.close()
        self._temporary.cleanup()


def _open_fastq(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="strict", newline="")
    return path.open("rt", encoding="ascii", errors="strict", newline="")


def _strip_line(value: str) -> str:
    return value.rstrip("\r\n")


_ASCII_WHITESPACE = bytes(code for code in range(128) if chr(code).isspace())
_MINIMUM_QUALITY_BYTE = 33
_MAXIMUM_QUALITY_BYTE = 126


def _records(path: Path) -> Iterator[tuple[str, int, bytes]]:
    try:
        with _open_fastq(path) as handle:
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
                    raise FastqValidationError(f"truncated FASTQ record {record_number}")
                header = _strip_line(header)
                sequence = _strip_line(sequence)
                separator = _strip_line(separator)
                quality = _strip_line(quality)
                if not header.startswith("@") or len(header) == 1:
                    raise FastqValidationError(f"invalid header in FASTQ record {record_number}")
                if not separator.startswith("+"):
                    raise FastqValidationError(f"invalid separator in FASTQ record {record_number}")
                if len(sequence) != len(quality):
                    raise FastqValidationError(
                        f"quality length mismatch in FASTQ record {record_number}"
                    )
                # The per-character scans below run once per base over ~1 TB of
                # delivered reads, so they use whole-buffer byte operations
                # rather than Python-level character loops. The accepted and
                # rejected character sets are unchanged.
                try:
                    sequence_bytes = sequence.encode("ascii")
                except UnicodeEncodeError:
                    raise FastqValidationError(
                        f"invalid sequence characters in FASTQ record {record_number}"
                    ) from None
                if len(sequence_bytes.translate(None, delete=_ASCII_WHITESPACE)) != len(
                    sequence_bytes
                ):
                    raise FastqValidationError(f"invalid sequence characters in FASTQ record {record_number}")
                try:
                    quality_bytes = quality.encode("ascii")
                except UnicodeEncodeError:
                    raise FastqValidationError(
                        f"invalid quality characters in FASTQ record {record_number}"
                    ) from None
                if quality_bytes and (
                    min(quality_bytes) < _MINIMUM_QUALITY_BYTE
                    or max(quality_bytes) > _MAXIMUM_QUALITY_BYTE
                ):
                    raise FastqValidationError(f"invalid quality characters in FASTQ record {record_number}")
                yield header[1:].split(maxsplit=1)[0], len(sequence_bytes), quality_bytes
    except (gzip.BadGzipFile, EOFError, UnicodeDecodeError) as error:
        raise FastqValidationError("malformed or truncated gzip FASTQ stream") from error


def validate_fastq(path: Path, tracker: ReadIdTracker | None = None, sample_id: str = "") -> FastqSummary:
    """Read a FASTQ stream sequentially and return non-sensitive quality summaries."""
    read_count = base_count = quality_total = duplicates = cross_sample = 0
    min_length: int | None = None
    max_length = 0
    for read_id, length, quality_bytes in _records(Path(path)):
        read_count += 1
        base_count += length
        min_length = length if min_length is None else min(min_length, length)
        max_length = max(max_length, length)
        quality_total += sum(quality_bytes) - 33 * len(quality_bytes)
        if tracker is not None:
            collision = tracker.observe(read_id, sample_id)
            duplicates += collision == "same_sample"
            cross_sample += collision == "cross_sample"
    return FastqSummary(
        read_count=read_count,
        base_count=base_count,
        min_read_length=min_length or 0,
        max_read_length=max_length,
        mean_read_length=(base_count / read_count) if read_count else 0.0,
        mean_base_quality=(quality_total / base_count) if base_count else 0.0,
        duplicate_read_ids=duplicates,
        cross_sample_read_ids=cross_sample,
    )


def sha256_file(path: Path) -> str:
    """Hash the delivered compressed file without retaining any sequence data."""
    digest = sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _safe_fastq_path(deliver_root: Path, relative_path: str) -> Path:
    candidate = (deliver_root / relative_path).resolve()
    if deliver_root not in candidate.parents or not candidate.is_file():
        raise FastqValidationError("expected FASTQ is missing outside or below the supplied delivery root")
    name = candidate.name.lower()
    if name.endswith(".corrupt.bak") or ".rebuild.bak." in name:
        raise FastqValidationError("backup FASTQ files are excluded from analysis")
    return candidate


def _write_table(path: Path, rows: list[dict[str, object]], fieldnames: tuple[str, ...]) -> None:
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_audit(deliver_root: Path, config_path: Path, output: Path) -> dict[str, object]:
    """Run the complete input audit and write only non-sensitive derived summaries."""
    runtime = load_runtime_config(deliver_root, config_path)
    output_dir = Path(output).resolve()
    safe_root = (runtime.project_root / "results" / "codex_reanalysis_v1").resolve()
    if output_dir != safe_root and safe_root not in output_dir.parents:
        raise ValueError("audit output must be inside results/codex_reanalysis_v1")
    output_dir.mkdir(parents=True, exist_ok=True)

    rows, expected_ids = build_rows_from_metadata(
        runtime.project_root / "data/meta/samplesheet.tsv",
        runtime.project_root / "data/meta/deliverable_mapping.tsv",
        runtime.project_root / "data/meta/drb3_alleles.tsv",
        runtime.deliver_root,
    )
    registry = build_sample_registry(rows, expected_ids)
    analytical = [row for row in registry if row["analysis_status"] == "analytical"]
    if len(analytical) != 13:
        raise SampleRegistryError("expected exactly 13 analytical manuscript samples")

    manifest: list[dict[str, object]] = []
    fastq_qc: list[dict[str, object]] = []
    tracker = ReadIdTracker()
    try:
        for row in analytical:
            path = _safe_fastq_path(runtime.deliver_root, row["input_relative_path"])
            summary = validate_fastq(path, tracker, row["sample_id"])
            manifest.append(
                {
                    "deliverable_id": row["deliverable_id"],
                    "sample_id": row["sample_id"],
                    "input_relative_path": row["input_relative_path"],
                    "byte_size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
            fastq_qc.append(
                {
                    "deliverable_id": row["deliverable_id"],
                    "sample_id": row["sample_id"],
                    "fastq_status": "pass",
                    **asdict(summary),
                }
            )
    finally:
        tracker.close()

    _write_table(
        output_dir / "input_manifest.tsv",
        manifest,
        ("deliverable_id", "sample_id", "input_relative_path", "byte_size", "sha256"),
    )
    _write_table(
        output_dir / "fastq_qc.tsv",
        fastq_qc,
        (
            "deliverable_id", "sample_id", "fastq_status", "read_count", "base_count",
            "min_read_length", "max_read_length", "mean_read_length", "mean_base_quality",
            "duplicate_read_ids", "cross_sample_read_ids",
        ),
    )
    write_tsv(output_dir / "sample_registry.tsv", registry)
    input_manifest_path = output_dir / "input_manifest.tsv"
    sample_registry_path = output_dir / "sample_registry.tsv"
    status = {
        "status": "pass",
        "analytical_sample_count": len(analytical),
        "input_manifest_sha256": sha256_file(input_manifest_path),
        "sample_registry_sha256": sha256_file(sample_registry_path),
        "absent_cohort_deliverables": ["DELIVERY_07", "DELIVERY_10"],
        "cross_sample_duplicate_read_ids": sum(int(row["cross_sample_read_ids"]) for row in fastq_qc),
    }
    (output_dir / "audit_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(
        output_dir,
        {
            "command": "audit_inputs.py --deliver-root [external] --config [tracked] --output [safe]",
            "config_hash": runtime.config_sha256,
            "input_hashes": {str(item["input_relative_path"]): str(item["sha256"]) for item in manifest},
            "reference_hashes": {},
        },
    )
    return status


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deliver-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        status = run_audit(args.deliver_root, args.config, args.output)
    except (FastqValidationError, SampleRegistryError, ValueError, OSError) as error:
        print(f"input audit blocked: {error}", file=sys.stderr)
        return 2
    print(f"input audit status: {status['status']}; analytical samples: {status['analytical_sample_count']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
