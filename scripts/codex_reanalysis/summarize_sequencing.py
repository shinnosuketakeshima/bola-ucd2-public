"""Per-animal sequencing and target-mapping summary (manuscript Table S1).

Joins three existing outputs and one measurement made here:

* read count, yield, mean read length and mean base quality of the delivered
  FASTQ (`00_input_audit/fastq_qc.tsv`);
* primary alignments retained in the target-interval BAM, and the median and
  N50 length of those primary reads (measured here with `samtools view
  -F 0x904`; these describe reads overlapping the analysis intervals, not the
  whole run);
* the control median depth at MAPQ >=20 used for normalization
  (`03_depth_cnv/normalized_depth_matrix.tsv`).

Whole-run read N50 is not reported: it needs the delivered FASTQ, which is
outside this workspace.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance


def n50(lengths: Sequence[int]) -> int:
    if not lengths:
        raise ValueError("n50 of an empty length list")
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def primary_read_lengths(samtools: str, bam: Path) -> list[int]:
    output = subprocess.run(
        [samtools, "view", "-F", "0x904", str(bam)], check=True, capture_output=True, text=True
    ).stdout
    return [len(line.split("\t", 10)[9]) for line in output.splitlines() if line]


def control_median_depths(matrix: Path, mapq: str = "20") -> dict[str, float]:
    depths: dict[str, float] = {}
    with matrix.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["mapq"] == mapq:
                depths.setdefault(row["sample_id"], float(row["control_median_depth"]))
    return depths


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq-qc", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--depth-matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="directory for sequencing_summary.tsv")
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    controls = control_median_depths(args.depth_matrix)
    rows = []
    with args.fastq_qc.open(encoding="utf-8", newline="") as handle:
        for qc in csv.DictReader(handle, delimiter="\t"):
            sample = qc["sample_id"]
            lengths = primary_read_lengths(args.samtools, args.bam_dir / f"{sample}.targets.bam")
            rows.append(
                {
                    "sample_id": sample,
                    "reads": int(qc["read_count"]),
                    "yield_gb": f"{int(qc['base_count']) / 1e9:.1f}",
                    "mean_read_length_kb": f"{float(qc['mean_read_length']) / 1000:.1f}",
                    "mean_base_quality": f"{float(qc['mean_base_quality']):.1f}",
                    "target_primary_alignments": len(lengths),
                    "target_read_median_kb": f"{statistics.median(lengths) / 1000:.1f}",
                    "target_read_n50_kb": f"{n50(lengths) / 1000:.1f}",
                    "control_median_depth_mapq20": f"{controls[sample]:g}",
                }
            )
    if len(rows) != 13:
        raise ValueError(f"expected 13 samples, found {len(rows)}")
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "sequencing_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_provenance(args.output, {"inputs": [str(args.fastq_qc), str(args.depth_matrix)], "bam_filter": "-F 0x904"})
    print(f"Sequencing summary written: {args.output / 'sequencing_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
