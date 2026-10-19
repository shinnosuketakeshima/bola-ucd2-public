"""Mapping checks behind the BTNL2 depth elevation (supplementary-alignment partners).

Depth alone cannot say whether extra reads at a locus are a tandem copy of
reference sequence or come from a rearranged or paralogous segment. For one
interval (default: the BTNL2 3' segment that is elevated in ID08), this stage
counts, in every animal's target BAM:

* primary reads (`-F 0x904`) at MAPQ >= `--min-mapq`, and the share with MAPQ
  < `--min-mapq` among all primary reads;
* secondary alignments retained in the BAM;
* primary reads whose SA tag places a supplementary segment with MAPQ >=
  `--min-mapq` in each `--partner` window.

A read is counted once per partner window. These are descriptive counts on
target-only BAMs: segments landing outside the retained intervals are not
visible, so absence of a partner is not evidence of absence.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance

DEFAULT_INTERVAL = "NC_037350.1:25867555-25883555"
DEFAULT_PARTNERS = (
    "NC_037350.1:25940001-25950000",
    "NC_037350.1:26150001-26170000",
    "NC_037350.1:26240001-26270000",
)


def parse_region(region: str) -> tuple[str, int, int]:
    contig, _, span = region.partition(":")
    start, _, end = span.partition("-")
    if not contig or not start or not end or int(start) > int(end):
        raise ValueError(f"region must be CONTIG:START-END, got {region!r}")
    return contig, int(start), int(end)


def partner_hits(sa_tag: str, partners: Sequence[tuple[str, int, int]], min_mapq: int) -> set[int]:
    """Indices of partner windows holding a supplementary segment of this read."""
    hits = set()
    for entry in sa_tag.split(";"):
        if not entry:
            continue
        contig, position, _strand, _cigar, mapq, _nm = entry.split(",")
        if int(mapq) < min_mapq:
            continue
        for index, (p_contig, p_start, p_end) in enumerate(partners):
            if contig == p_contig and p_start <= int(position) <= p_end:
                hits.add(index)
    return hits


def summarize_records(lines: Sequence[str], partners: Sequence[tuple[str, int, int]], min_mapq: int) -> dict:
    primary = low_mapq = secondary = kept = 0
    partner_counts = [0] * len(partners)
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        if flag & 0x100:
            secondary += 1
            continue
        if flag & 0x804:
            continue
        primary += 1
        if int(fields[4]) < min_mapq:
            low_mapq += 1
            continue
        kept += 1
        tag = next((t[5:] for t in fields[11:] if t.startswith("SA:Z:")), "")
        for index in partner_hits(tag, partners, min_mapq):
            partner_counts[index] += 1
    return {
        "primary_reads": primary,
        "low_mapq_fraction": f"{low_mapq / primary:.2f}" if primary else "",
        "secondary_alignments": secondary,
        "primary_mapq_reads": kept,
        "partner_counts": partner_counts,
    }


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", default=DEFAULT_INTERVAL)
    parser.add_argument("--partner", action="append", default=None)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    parse_region(args.interval)
    partner_regions = args.partner or list(DEFAULT_PARTNERS)
    partners = [parse_region(p) for p in partner_regions]
    with args.registry.open(encoding="utf-8", newline="") as handle:
        samples = [r["sample_id"] for r in csv.DictReader(handle, delimiter="\t") if r["analysis_status"] == "analytical"]
    if len(samples) != 13:
        raise ValueError(f"expected 13 analytical animals, found {len(samples)}")
    rows = []
    for sample in samples:
        output = subprocess.run(
            [args.samtools, "view", str(args.bam_dir / f"{sample}.targets.bam"), args.interval],
            check=True, capture_output=True, text=True,
        ).stdout
        summary = summarize_records(output.splitlines(), partners, args.min_mapq)
        row = {"sample_id": sample, **{k: v for k, v in summary.items() if k != "partner_counts"}}
        for region, count in zip(partner_regions, summary["partner_counts"]):
            row[f"SA_{region}"] = count
            row[f"SA_{region}_fraction"] = f"{count / summary['primary_mapq_reads']:.2f}" if summary["primary_mapq_reads"] else ""
        rows.append(row)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "btn_split_reads.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    status = {"status": "passed", "samples": len(rows), "interval": args.interval, "partners": partner_regions,
              "min_mapq": args.min_mapq}
    (args.output / "btn_split_reads_status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n",
                                                             encoding="utf-8")
    write_provenance(
        args.output,
        {"inputs": {"registry": str(args.registry)},
         "parameters": {"interval": args.interval, "partners": partner_regions, "min_mapq": args.min_mapq,
                        "primary_filter": "-F 0x904"}},
    )
    print(f"BTN split-read summary written: {args.output / 'btn_split_reads.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
