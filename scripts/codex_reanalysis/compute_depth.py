"""Extract per-locus depth from one sample's target BAM at MAPQ 0/10/20.

Zero-coverage positions are retained deliberately: dropping them (as
``samtools depth`` does without ``-a``) would silently overstate a locus's
median depth for any paralogue the reference under-represents, which is
exactly the failure mode ``docs/HANDOFF_PLAN_PROGRESS.md`` §2b documents for
the legacy single-reference depth-CNV pipeline.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

DEFAULT_MAPQ_PASSES = (0, 10, 20)


def build_depth_command(
    bam: str, contig: str, start: int, end: int, mapq: int, samtools: str = "samtools"
) -> list[str]:
    """Build a ``samtools depth`` command over one 0-based half-open BED interval."""
    if end <= start:
        raise ValueError(f"invalid interval: {start}-{end}")
    region = f"{contig}:{start + 1}-{end}"
    return [samtools, "depth", "-a", "-Q", str(mapq), "-r", region, bam]


def parse_depth_output(text: str) -> list[tuple[int, int]]:
    """Parse ``samtools depth`` lines (contig, pos, depth) into (pos, depth) pairs."""
    positions: list[tuple[int, int]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        _, pos, depth = line.split("\t")
        positions.append((int(pos), int(depth)))
    return positions


def run_depth(
    bam: str, contig: str, start: int, end: int, mapq: int, samtools: str = "samtools"
) -> list[tuple[int, int]]:
    command = build_depth_command(bam, contig, start, end, mapq, samtools)
    result = subprocess.run(command, capture_output=True, check=True, text=True)
    return parse_depth_output(result.stdout)


def sha256_paths(paths: Mapping[str, Path]) -> dict[str, str]:
    """Return SHA-256 hex digests keyed by a stable artefact name."""
    return {name: sha256_file(path) for name, path in paths.items()}


def require_complete_depth(
    positions: Sequence[tuple[int, int]], expected_length: int, locus: str, mapq: int
) -> None:
    """Fail if ``samtools depth -a`` did not return every base in the BED interval."""
    if len(positions) != expected_length:
        raise ValueError(
            f"{locus} MAPQ {mapq}: expected {expected_length} positions, got {len(positions)}"
        )


def summarize_depth(positions: Sequence[tuple[int, int]]) -> dict:
    """Minimal, zero-retaining depth summary: median depth and zero-coverage count."""
    depths = [depth for _, depth in positions]
    if not depths:
        return {"median": 0.0, "zero_positions": 0}
    return {"median": float(statistics.median(depths)), "zero_positions": sum(1 for d in depths if d == 0)}


def depth_row(positions: Sequence[tuple[int, int]], expected_length: int) -> dict:
    """Full per-locus depth summary used by the normalization stage."""
    summary = summarize_depth(positions)
    depths = [depth for _, depth in positions]
    mean = float(statistics.mean(depths)) if depths else 0.0
    mad = float(statistics.median(abs(d - summary["median"]) for d in depths)) if depths else 0.0
    breadth = 1.0 - (summary["zero_positions"] / expected_length) if expected_length else 0.0
    return {
        "median": summary["median"],
        "mean": round(mean, 3),
        "mad": round(mad, 3),
        "zero_positions": summary["zero_positions"],
        "n_positions": len(depths),
        "breadth": round(breadth, 4),
    }


def read_bed(bed_path: Path) -> list[tuple[str, int, int, str]]:
    targets = []
    with Path(bed_path).open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            contig, start, end, name = row[0], int(row[1]), int(row[2]), row[3]
            targets.append((contig, start, end, name))
    return targets


def compute_sample_depth(
    bam: Path,
    bed_path: Path,
    mapq_passes: Sequence[int] = DEFAULT_MAPQ_PASSES,
    samtools: str = "samtools",
) -> list[dict]:
    rows = []
    for contig, start, end, name in read_bed(bed_path):
        length = end - start
        for mapq in mapq_passes:
            positions = run_depth(str(bam), contig, start, end, mapq, samtools)
            require_complete_depth(positions, length, name, mapq)
            row = depth_row(positions, length)
            row.update({"locus": name, "contig": contig, "start": start, "end": end, "mapq": mapq})
            rows.append(row)
    return rows


def write_depth_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "locus",
        "contig",
        "start",
        "end",
        "mapq",
        "median",
        "mean",
        "mad",
        "breadth",
        "zero_positions",
        "n_positions",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path, help="Analytical per-locus BED (targets.bed)")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--output", required=True, type=Path, help="Per-sample output directory")
    parser.add_argument("--mapq", nargs="+", type=int, default=list(DEFAULT_MAPQ_PASSES))
    args = parser.parse_args(argv)

    if not args.bam.is_file():
        raise FileNotFoundError(f"target BAM does not exist: {args.bam}")

    rows = compute_sample_depth(args.bam, args.bed, args.mapq, args.samtools)
    args.output.mkdir(parents=True, exist_ok=True)
    write_depth_tsv(args.output / "depth_raw.tsv", rows)

    status_payload = {
        "sample_id": args.sample_id,
        "locus_count": len({row["locus"] for row in rows}),
        "mapq_passes": list(args.mapq),
        "row_count": len(rows),
        "status": "passed",
    }
    (args.output / "depth_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        args.output,
        {
            "input_hashes": sha256_paths({"bam": args.bam, "bed": args.bed}),
            "reference_hashes": sha256_paths({"bed": args.bed}),
        },
    )
    print(f"depth computed for {args.sample_id}: {len(rows)} rows -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
