"""Verify the delivered FASTQ by content once, and record it durably.

Re-hashing 1.13 TiB takes hours, so it must not be tied to the short lifetime
of a host check. This seal is content-bound rather than time-bound: it stays
valid while the delivered bytes and the audit manifest it was verified against
are unchanged, and mapping records that its input hash was verified live on
this host rather than merely carried over from the audit.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Mapping

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


class InputSealError(RuntimeError):
    """Raised when delivered inputs cannot be verified against the audit."""


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file():
        raise InputSealError(f"input manifest is absent: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if not rows:
        raise InputSealError("input manifest is empty")
    return rows


def build_input_seal(
    deliver_root: Path,
    manifest_path: Path,
    *,
    progress=None,
) -> dict[str, object]:
    """Hash every delivered FASTQ and confirm it against the audit manifest."""
    root = Path(deliver_root).resolve()
    rows = read_manifest(manifest_path)
    samples: dict[str, dict[str, object]] = {}
    for row in rows:
        sample_id = row["sample_id"]
        path = (root / row["input_relative_path"]).resolve()
        if root not in path.parents:
            raise InputSealError(f"delivered path escapes the delivery root: {sample_id}")
        if not path.is_file():
            raise InputSealError(f"delivered FASTQ is absent: {path}")
        size = path.stat().st_size
        if size != int(row["byte_size"]):
            raise InputSealError(f"delivered FASTQ size changed since the audit: {sample_id}")
        digest = sha256_file(path)
        if digest.lower() != row["sha256"].lower():
            raise InputSealError(
                f"delivered FASTQ content changed since the audit: {sample_id}"
            )
        samples[sample_id] = {"sha256": digest, "byte_size": size}
        if progress is not None:
            progress(sample_id, size)
    return {
        "status": "verified",
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "delivery_root": str(root),
        "input_manifest_sha256": sha256_file(Path(manifest_path)),
        "sample_count": len(samples),
        "samples": samples,
    }


def validate_input_seal(seal: Mapping[str, object], manifest_path: Path) -> None:
    """Confirm a seal is complete and still describes the current audit manifest.

    This is cheap on purpose: it is the check a host runs before every long
    job, and it must not re-read a terabyte to do it.
    """
    if seal.get("status") != "verified":
        raise InputSealError(f"input seal is not verified: {seal.get('status')}")
    manifest_hash = sha256_file(Path(manifest_path))
    if seal.get("input_manifest_sha256") != manifest_hash:
        raise InputSealError(
            "input seal was built against a different audit manifest; re-verify the inputs"
        )
    samples = seal.get("samples")
    if not isinstance(samples, dict) or not samples:
        raise InputSealError("input seal records no samples")
    expected = {row["sample_id"] for row in read_manifest(manifest_path)}
    missing = sorted(expected - set(samples))
    if missing:
        raise InputSealError(f"input seal is missing samples: {', '.join(missing)}")


def verified_sample(seal: Mapping[str, object], sample_id: str) -> dict[str, object]:
    samples = seal.get("samples")
    if not isinstance(samples, dict) or sample_id not in samples:
        raise InputSealError(f"input seal does not cover {sample_id}")
    record = samples[sample_id]
    if not isinstance(record, dict) or not record.get("sha256"):
        raise InputSealError(f"input seal record is incomplete for {sample_id}")
    return dict(record)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deliver-root", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    def report(sample_id: str, size: int) -> None:
        sys.stderr.write(f"  verified {sample_id:14s} {size / 1e9:6.1f} GB\n")
        sys.stderr.flush()

    try:
        seal = build_input_seal(args.deliver_root, args.input_manifest, progress=report)
    except InputSealError as error:
        sys.stderr.write(f"INPUT VERIFICATION FAILED: {error}\n")
        return 1
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(f"{output.suffix}.partial")
    temporary.write_text(json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    sys.stderr.write(f"\nverified {seal['sample_count']} delivered FASTQ against the audit\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
