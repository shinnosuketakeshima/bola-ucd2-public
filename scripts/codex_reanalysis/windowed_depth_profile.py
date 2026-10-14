"""1-kb-bin windowed depth profiling of BTNL2 and BTN2A1-like on ARS-UCD2.0.

Ports the sub-locus breakpoint-shape observation behind manuscript Fig. 3B/C
(BTNL2 and BTN2A1-like/LOC104969893) from the ARS-UCD1.2-generation ad-hoc
analysis recorded in docs/HANDOFF_PLAN_PROGRESS.md SS2c-2d
(results/cnv_corrected/) to the ARS-UCD2.0 codex_reanalysis pipeline: same
gene-level depth signal as 03_depth_cnv (samtools depth -a, so zero-coverage
positions are retained rather than silently dropped), at finer 1-kb-bin
resolution, over the full BTNL2 and BTN2A1-like gene spans, normalized by
each sample's already-computed MAPQ>=20 control-region median depth
(03_depth_cnv/normalized_depth_matrix.tsv) so bin-level and gene-level ratios
are on the same scale.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
from pathlib import Path
from typing import Sequence

MAPQ = 20
BIN_SIZE = 1000

# (locus label, contig, 0-based start, 0-based end) -- coordinate_registry.tsv
LOCI = [
    ("BTNL2", "NC_037350.1", 25_858_353, 25_883_555),
    ("BTN2A1L_LOC104969893", "NC_037350.1", 26_037_136, 26_114_879),
]

SAMPLE_ORDER = [
    "ID01",
    "ID02",
    "ID03",
    "ID04",
    "ID05",
    "ID06",
    "ID07",
    "ID08",
    "ID09",
    "ID10",
    "ID11",
    "ID12",
    "ID13",
]


def run_depth(bam: Path, contig: str, start: int, end: int, samtools: str) -> list[int]:
    """samtools depth -a -Q 20 over a 0-based half-open interval; returns per-base depth."""
    region = f"{contig}:{start + 1}-{end}"
    command = [samtools, "depth", "-a", "-Q", str(MAPQ), "-r", region, str(bam)]
    result = subprocess.run(command, capture_output=True, check=True, text=True)
    depths: list[int] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        _, _, depth = line.split("\t")
        depths.append(int(depth))
    expected = end - start
    if len(depths) != expected:
        raise ValueError(f"{contig}:{start}-{end}: expected {expected} positions, got {len(depths)}")
    return depths


def bin_depths(depths: Sequence[int], start: int, bin_size: int = BIN_SIZE) -> list[tuple[int, int, float]]:
    """Mean depth per fixed-size bin from a gene's start; the final bin may be shorter."""
    rows = []
    for offset in range(0, len(depths), bin_size):
        chunk = depths[offset : offset + bin_size]
        bin_start = start + offset
        bin_end = bin_start + len(chunk)
        rows.append((bin_start, bin_end, statistics.mean(chunk)))
    return rows


def load_control_median(matrix_path: Path) -> dict[str, float]:
    """sample_id -> control_median_depth at MAPQ 20, read from the existing 03_depth_cnv matrix."""
    control_median: dict[str, float] = {}
    with matrix_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["mapq"]) != MAPQ:
                continue
            sample_id = row["sample_id"]
            value = float(row["control_median_depth"])
            if sample_id in control_median and control_median[sample_id] != value:
                raise ValueError(f"inconsistent control_median_depth for {sample_id}")
            control_median[sample_id] = value
    return control_median


def profile_sample(
    sample_id: str, bam: Path, samtools: str
) -> list[dict]:
    rows = []
    for locus, contig, start, end in LOCI:
        depths = run_depth(bam, contig, start, end, samtools)
        for bin_start, bin_end, mean_depth in bin_depths(depths, start):
            rows.append(
                {
                    "locus": locus,
                    "sample_id": sample_id,
                    "contig": contig,
                    "bin_start": bin_start,
                    "bin_end": bin_end,
                    "mean_depth": round(mean_depth, 3),
                }
            )
    return rows


def write_tsv(path: Path, rows: Sequence[dict], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, default=None, help="defaults to <root>/results/codex_reanalysis_v1/02_mapping")
    parser.add_argument("--depth-matrix", type=Path, default=None, help="defaults to <root>/results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv")
    parser.add_argument("--output", type=Path, default=None, help="defaults to <root>/results/codex_reanalysis_v1/03_depth_cnv/windowed")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args()

    root = args.project_root.resolve()
    mapping_dir = args.mapping_dir or root / "results/codex_reanalysis_v1/02_mapping"
    depth_matrix = args.depth_matrix or root / "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv"
    output = args.output or root / "results/codex_reanalysis_v1/03_depth_cnv/windowed"
    output.mkdir(parents=True, exist_ok=True)

    control_median = load_control_median(depth_matrix)
    missing = [sample for sample in SAMPLE_ORDER if sample not in control_median]
    if missing:
        raise SystemExit(f"normalized_depth_matrix.tsv missing control_median_depth for: {missing}")

    all_rows: list[dict] = []
    for sample_id in SAMPLE_ORDER:
        bam = mapping_dir / f"{sample_id}.targets.bam"
        if not bam.is_file():
            raise SystemExit(f"missing target BAM: {bam}")
        rows = profile_sample(sample_id, bam, args.samtools)
        sample_control_median = control_median[sample_id]
        for row in rows:
            row["control_median_depth"] = sample_control_median
            row["normalized_ratio"] = round(row["mean_depth"] / sample_control_median, 4) if sample_control_median > 0 else None
        all_rows.extend(rows)
        print(f"{sample_id}: {len(rows)} bins")

    fieldnames = ["locus", "sample_id", "contig", "bin_start", "bin_end", "mean_depth", "control_median_depth", "normalized_ratio"]
    write_tsv(output / "windowed_depth_bins.tsv", all_rows, fieldnames)

    samtools_version = subprocess.run(
        [args.samtools, "--version"], capture_output=True, check=True, text=True
    ).stdout.splitlines()[0]

    status = {
        "status": "passed",
        "loci": [locus for locus, *_ in LOCI],
        "sample_count": len(SAMPLE_ORDER),
        "bin_size": BIN_SIZE,
        "mapq": MAPQ,
        "samtools_version": samtools_version,
        "row_count": len(all_rows),
    }
    (output / "windowed_depth_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(all_rows)} bins for {len(SAMPLE_ORDER)} samples -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
