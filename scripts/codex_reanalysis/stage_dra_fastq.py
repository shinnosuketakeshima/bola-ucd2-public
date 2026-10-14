"""Stage public-named gzip FASTQ for DDBJ DRA upload on the analysis host.

This script does not upload files, contact DDBJ, or mint accessions. It:

- reads a private map of manuscript ID → local analytical FASTQ
- refuses a source inside the immutable delivery root unless
  ``--sanitize-headers`` is set, since only the sanitizing rewrite reads the
  delivery root without ever writing, renaming, copying, or linking inside it
- optionally rewrites FASTQ headers so vendor sample tags never enter INSDC
- writes ``ID01.ont.fastq.gz`` … ``ID13.ont.fastq.gz`` plus filled MD5 tables
  into a stage directory outside the git tree

The private map belongs under ``data/meta/private/`` (gitignored). Never copy
it into the stage directory or into tracked worksheets.

Run on the analysis host that holds (or can read) the mapping-input FASTQ.
Without ``--sanitize-headers``, a source path must never resolve into
``Deliver_ADD``; hardlink/copy staging is for a local analytical copy only.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.build_evidence_matrix import contains_real_sample_id
from scripts.codex_reanalysis.write_ddbj_metadata import (
    COHORT,
    RUN_FIELDS,
    planned_fastq_name,
)

EXPECTED_IDS = tuple(row["sample_id"] for row in COHORT)
DELIVERY_MARKER = "Deliver_ADD"


class StageError(RuntimeError):
    """Raised when a FASTQ cannot be staged for a public DRA upload."""


def _posix(path: Path) -> str:
    return str(path).replace("\\", "/")


def is_delivery_path(path: Path) -> bool:
    """True if *path* is inside the immutable delivery root, even via a symlink."""
    candidates = [_posix(path)]
    try:
        candidates.append(_posix(path.resolve()))
    except OSError:
        pass
    return any(DELIVERY_MARKER in text for text in candidates)


def read_fastq_map(path: Path) -> dict[str, Path]:
    """Read ``manuscript_id`` / ``fastq`` pairs. The map itself is not copied out."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise StageError("FASTQ map is empty")
    missing = [name for name in ("manuscript_id", "fastq") if name not in (rows[0].keys() or ())]
    if missing:
        raise StageError(f"FASTQ map lacks columns: {', '.join(missing)}")
    mapping: dict[str, Path] = {}
    for row in rows:
        sample_id = str(row.get("manuscript_id", "")).strip()
        fastq = str(row.get("fastq", "")).strip()
        if not sample_id and not fastq:
            continue
        if sample_id in mapping:
            raise StageError(f"FASTQ map lists {sample_id} more than once")
        mapping[sample_id] = Path(fastq)
    expected = list(EXPECTED_IDS)
    got = list(mapping)
    if got != expected:
        raise StageError(
            "FASTQ map must list ID01–ID13 once each, in order; "
            f"got {got!r}"
        )
    return mapping


def public_read_id(header_line: str, sample_id: str, index: int) -> str:
    """Keep the first FASTQ token, or mint a public name if that token leaks identity."""
    body = header_line[1:] if header_line.startswith("@") else header_line
    token = body.split()[0] if body.strip() else ""
    if token and not contains_real_sample_id(token):
        return token
    return f"{sample_id}.read{index}"


def header_leaks_identity(header_line: str) -> bool:
    return contains_real_sample_id(header_line)


def scan_headers(path: Path, limit: int | None) -> int:
    """Return how many scanned headers contain a real sample-id token."""
    leaked = 0
    with gzip.open(path, "rt", encoding="ascii", errors="replace", newline="") as handle:
        record = 0
        while True:
            header = handle.readline()
            if header == "":
                break
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise StageError(f"truncated FASTQ in {path.name}")
            record += 1
            if header_leaks_identity(header):
                leaked += 1
            if limit is not None and limit > 0 and record >= limit:
                break
    return leaked


def write_sanitized_gzip(
    source: Path, destination: Path, sample_id: str, compressor: str | None = None
) -> int:
    """Stream a gzip FASTQ, stripping header tags. Returns the number of reads.

    ``compressor`` names an external program (for example ``pigz``) that reads
    the sanitized text on stdin and writes gzip on stdout, used in place of
    Python's single-threaded ``gzip`` module for the write side only; the
    delivery-root source is always read through the standard library.
    """
    partial = destination.with_suffix(destination.suffix + ".partial")
    partial.unlink(missing_ok=True)
    reads = 0
    writer_process: subprocess.Popen | None = None
    try:
        with gzip.open(source, "rt", encoding="ascii", errors="replace", newline="") as incoming:
            if compressor:
                handle = partial.open("wb")
                writer_process = subprocess.Popen(
                    [compressor, "-c"], stdin=subprocess.PIPE, stdout=handle, text=True
                )
                handle.close()
                outgoing = writer_process.stdin
            else:
                outgoing = gzip.open(partial, "wt", encoding="ascii", newline="", compresslevel=6)
            try:
                while True:
                    header = incoming.readline()
                    if header == "":
                        break
                    sequence = incoming.readline()
                    separator = incoming.readline()
                    quality = incoming.readline()
                    if not sequence or not separator or not quality:
                        raise StageError(f"truncated FASTQ while sanitizing {sample_id}")
                    if not header.startswith("@") or not separator.startswith("+"):
                        raise StageError(f"malformed FASTQ record while sanitizing {sample_id}")
                    reads += 1
                    outgoing.write("@" + public_read_id(header.rstrip("\r\n"), sample_id, reads) + "\n")
                    outgoing.write(sequence if sequence.endswith("\n") else sequence + "\n")
                    outgoing.write("+\n")
                    outgoing.write(quality if quality.endswith("\n") else quality + "\n")
            finally:
                outgoing.close()
        if writer_process is not None and writer_process.wait() != 0:
            raise StageError(f"{sample_id}: external compressor failed")
    except BaseException:
        if writer_process is not None and writer_process.poll() is None:
            writer_process.kill()
            writer_process.wait()
        partial.unlink(missing_ok=True)
        raise
    os.replace(partial, destination)
    return reads


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _link_or_copy(source: Path, destination: Path) -> str:
    destination.unlink(missing_ok=True)
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        return "copy"


def stage_one(
    sample_id: str,
    source: Path,
    stage_dir: Path,
    *,
    sanitize_headers: bool,
    header_scan_reads: int,
    compressor: str | None = None,
) -> dict[str, str]:
    if is_delivery_path(source) and not sanitize_headers:
        raise StageError(
            f"{sample_id}: refusing a path inside the delivery root without --sanitize-headers; "
            "a delivery-root source may only be read through the sanitizing rewrite, "
            "which never writes, renames, or links inside the delivery root"
        )
    if not source.is_file():
        raise StageError(f"{sample_id}: FASTQ is absent: {source}")
    destination = stage_dir / planned_fastq_name(sample_id)
    if sanitize_headers:
        reads = write_sanitized_gzip(source, destination, sample_id, compressor=compressor)
        leaked = scan_headers(destination, limit=None)
        if leaked:
            raise StageError(f"{sample_id}: sanitized FASTQ still has identity-bearing headers")
        method = f"sanitize_headers:{reads}"
    else:
        leaked = scan_headers(source, limit=header_scan_reads)
        if leaked:
            raise StageError(
                f"{sample_id}: {leaked} scanned header(s) contain a real sample-id token; "
                "re-run with --sanitize-headers"
            )
        method = _link_or_copy(source, destination)
    return {
        "sample_name": sample_id,
        "file_name": destination.name,
        "file_type": "fastq",
        "md5_checksum": md5_file(destination),
        "stage_method": method,
        "byte_size": str(destination.stat().st_size),
    }


def write_run_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in RUN_FIELDS})


def write_sftp_commands(path: Path, stage_dir: Path) -> None:
    files = " ".join(planned_fastq_name(sample_id) for sample_id in EXPECTED_IDS)
    text = (
        "# Fill ACCOUNT and SUBMISSION_ID from D-way after creating the DRA submission.\n"
        f"# Staged gzip files are in: {stage_dir}\n"
        "sftp -o ServerAliveInterval=60 -o TCPKeepAlive=yes "
        "-i ~/.ssh/id_ed25519 ACCOUNT@ftp-private.ddbj.nig.ac.jp\n"
        "sftp> cd SUBMISSION_ID\n"
        f"sftp> lcd {stage_dir}\n"
        f"sftp> mput {files}\n"
    )
    Path(path).write_text(text, encoding="utf-8")


def stage_dra_fastq(
    fastq_map: Path,
    stage_dir: Path,
    *,
    sanitize_headers: bool,
    header_scan_reads: int,
    compressor: str | None = None,
) -> list[dict[str, str]]:
    mapping = read_fastq_map(fastq_map)
    stage_dir = Path(stage_dir)
    if is_delivery_path(stage_dir):
        raise StageError("stage directory must not sit inside the delivery root")
    stage_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for sample_id, source in mapping.items():
        rows.append(
            stage_one(
                sample_id,
                source,
                stage_dir,
                sanitize_headers=sanitize_headers,
                header_scan_reads=header_scan_reads,
                compressor=compressor,
            )
        )
    write_run_tsv(stage_dir / "dra_run.filled.tsv", rows)
    write_sftp_commands(stage_dir / "sftp_commands.txt", stage_dir)
    summary = stage_dir / "stage_summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_name", "file_name", "md5_checksum", "byte_size", "stage_method"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in writer.fieldnames})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fastq-map",
        type=Path,
        required=True,
        help="private TSV with columns manuscript_id and fastq (gitignored)",
    )
    parser.add_argument(
        "--stage-dir",
        type=Path,
        required=True,
        help="directory for ID01.ont.fastq.gz … ID13.ont.fastq.gz (not inside git)",
    )
    parser.add_argument(
        "--sanitize-headers",
        action="store_true",
        help="rewrite FASTQ headers to drop vendor sample tags before upload",
    )
    parser.add_argument(
        "--header-scan-reads",
        type=int,
        default=10_000,
        help="headers to scan when --sanitize-headers is not set (0 = all)",
    )
    parser.add_argument(
        "--compressor",
        default=None,
        help="external gzip-compatible program (for example pigz) for the write side "
        "of --sanitize-headers; the delivery-root source is always read with gzip",
    )
    args = parser.parse_args()
    if args.compressor and shutil.which(args.compressor) is None:
        raise StageError(f"compressor is not available: {args.compressor}")
    stage_dra_fastq(
        args.fastq_map,
        args.stage_dir,
        sanitize_headers=args.sanitize_headers,
        header_scan_reads=args.header_scan_reads,
        compressor=args.compressor,
    )
    print(f"staged 13 FASTQ under {args.stage_dir}")
    print(f"filled MD5 table: {args.stage_dir / 'dra_run.filled.tsv'}")
    print("copy md5_checksum into docs/codex_reanalysis/ddbj_metadata/dra_run.tsv after upload")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
