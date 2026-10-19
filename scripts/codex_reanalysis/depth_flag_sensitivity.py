"""Recompute cohort-relative depth flags with one animal left out.

The cohort median and MAD behind every depth flag (`normalize_depth.py`) are
computed from all 13 animals, so an animal whose identity or sex is in doubt
also shapes the flags of the others. This stage removes `--exclude`
(e.g. ID02, whose depth-based sex disagrees with its record), recomputes the
median, MAD and flag of every remaining animal with the same
`interpretation_flag` rule, and lists each flag that appears or disappears.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.normalize_depth import CONSISTENT, UNREPRESENTED_EXCEPTION, interpretation_flag
from scripts.codex_reanalysis.provenance import write_provenance

NOT_FLAGGED = {CONSISTENT, UNREPRESENTED_EXCEPTION}


def recompute_flags(ratios: dict[str, float], controls: dict[str, float], minimum_control_depth: float) -> dict[str, str]:
    values = list(ratios.values())
    median = statistics.median(values)
    mad = statistics.median(abs(v - median) for v in values)
    return {
        sample: interpretation_flag(ratio, median, mad, controls[sample], minimum_control_depth)
        for sample, ratio in ratios.items()
    }


def compare(matrix_rows: Sequence[dict], exclude: str, mapq: str, minimum_control_depth: float) -> list[dict]:
    by_locus: dict[str, list[dict]] = {}
    for row in matrix_rows:
        if row["mapq"] == mapq:
            by_locus.setdefault(row["locus"], []).append(row)
    if not any(r["sample_id"] == exclude for rows in by_locus.values() for r in rows):
        raise ValueError(f"{exclude} is not in the depth matrix")
    changes = []
    for locus, rows in sorted(by_locus.items()):
        kept = [r for r in rows if r["sample_id"] != exclude]
        ratios = {r["sample_id"]: float(r["normalized_ratio"]) for r in kept}
        controls = {r["sample_id"]: float(r["control_median_depth"]) for r in kept}
        reduced = recompute_flags(ratios, controls, minimum_control_depth)
        for row in kept:
            before, after = row["interpretation"], reduced[row["sample_id"]]
            if (before in NOT_FLAGGED) != (after in NOT_FLAGGED) or (before not in NOT_FLAGGED and before != after):
                changes.append(
                    {"locus": locus, "sample_id": row["sample_id"], "ratio": row["normalized_ratio"],
                     "flag_all_animals": before, f"flag_without_{exclude}": after}
                )
    return changes


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="03_depth_cnv/normalized_depth_matrix.tsv")
    parser.add_argument("--exclude", required=True, help="sample_id to leave out")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mapq", default="20")
    parser.add_argument("--minimum-control-depth", type=float, default=10.0)
    args = parser.parse_args(argv)
    with args.matrix.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    changes = compare(rows, args.exclude, args.mapq, args.minimum_control_depth)
    args.output.mkdir(parents=True, exist_ok=True)
    fields = ["locus", "sample_id", "ratio", "flag_all_animals", f"flag_without_{args.exclude}"]
    with (args.output / f"flag_changes_without_{args.exclude}.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(changes)
    status = {"status": "passed", "excluded": args.exclude, "mapq": args.mapq, "flag_changes": len(changes)}
    (args.output / "flag_sensitivity_status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n",
                                                              encoding="utf-8")
    write_provenance(
        args.output,
        {"inputs": {"matrix": str(args.matrix)},
         "parameters": {"exclude": args.exclude, "mapq": args.mapq, "minimum_control_depth": args.minimum_control_depth}},
    )
    print(f"{len(changes)} flag changes without {args.exclude}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
