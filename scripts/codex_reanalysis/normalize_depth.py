"""Normalize per-sample depth against distributed single-copy controls.

Aggregates the per-sample ``depth_raw.tsv`` tables written by
``compute_depth.py`` into one cohort-level normalized matrix and report. A
ratio is reported for every locus/sample/MAPQ combination; it is never
collapsed into an integer copy-number call here, because a single linear
reference cannot distinguish "diploid at reduced identity" from "haploid
loss" for a paralogous BoLA locus without independent structural evidence
(assembly, read-multihit, or wet typing) -- see
``docs/HANDOFF_PLAN_PROGRESS.md`` §2b-§2h for the concrete cases where that
distinction mattered.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

CONTROL_PREFIX = "control_"
OUTLIER_THRESHOLD = 3.5
UNREPRESENTED_EXCEPTION = "present_at_otherwise_unrepresented_locus"
MODIFIED_Z_ELEVATED = "elevated_relative_to_cohort"
MODIFIED_Z_REDUCED = "reduced_relative_to_cohort"
INSUFFICIENT_CONTROL = "insufficient_control_signal"
CONSISTENT = "consistent_with_cohort"


def normalized_ratio(locus_median: float, control_median: float) -> float:
    """Return locus depth relative to the sample's control-locus median depth."""
    if control_median <= 0:
        raise ValueError("control median depth must be positive")
    return locus_median / control_median


def control_median_depth(control_locus_medians: Mapping[str, float]) -> float:
    """Summarize the sample's distributed single-copy controls into one median.

    Using the median across control loci (rather than any single control
    gene) means one poorly mapping control does not silently shift every
    normalized ratio in the sample.
    """
    values = list(control_locus_medians.values())
    if not values:
        raise ValueError("at least one control locus is required")
    return float(statistics.median(values))


def interpretation_flag(ratio: float, cohort_median: float, cohort_mad: float, control_median: float, minimum_control_depth: float) -> str:
    """Classify a ratio without asserting an integer copy number.

    ``insufficient_control_signal`` takes priority: a low-depth control
    makes the ratio itself unreliable regardless of where it falls.

    When the cohort MAD is zero the modified z-score is undefined. A
    zero-median locus with a non-zero sample is reference representation
    (or haplotype match), not copy-number elevation.
    """
    if control_median < minimum_control_depth:
        return INSUFFICIENT_CONTROL
    if cohort_mad == 0:
        if ratio == cohort_median:
            return CONSISTENT
        if cohort_median == 0 and ratio > 0:
            return UNREPRESENTED_EXCEPTION
        return MODIFIED_Z_ELEVATED if ratio > cohort_median else MODIFIED_Z_REDUCED
    modified_z = 0.6745 * (ratio - cohort_median) / cohort_mad
    if abs(modified_z) <= OUTLIER_THRESHOLD:
        return CONSISTENT
    return MODIFIED_Z_ELEVATED if modified_z > 0 else MODIFIED_Z_REDUCED


def read_depth_raw(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["mapq"] = int(row["mapq"])
        row["median"] = float(row["median"])
        row["breadth"] = float(row["breadth"])
    return rows


def load_cohort_depth(root: Path, sample_ids: Sequence[str]) -> dict[str, list[dict]]:
    cohort = {}
    for sample_id in sample_ids:
        path = Path(root) / sample_id / "depth_raw.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"missing depth_raw.tsv for sample: {sample_id}")
        cohort[sample_id] = read_depth_raw(path)
    return cohort


def build_normalized_matrix(cohort_depth: Mapping[str, list[dict]], minimum_control_depth: float = 10.0) -> list[dict]:
    """Build one ratio row per sample/locus/MAPQ, then flag each against the cohort."""
    by_sample_mapq_locus: dict[tuple[str, int], dict[str, dict]] = {}
    for sample_id, rows in cohort_depth.items():
        for row in rows:
            key = (sample_id, row["mapq"])
            by_sample_mapq_locus.setdefault(key, {})[row["locus"]] = row

    sample_control_median: dict[tuple[str, int], float] = {}
    for (sample_id, mapq), loci in by_sample_mapq_locus.items():
        controls = {name: row["median"] for name, row in loci.items() if name.startswith(CONTROL_PREFIX)}
        sample_control_median[(sample_id, mapq)] = control_median_depth(controls)

    unnormalized: list[dict] = []
    for (sample_id, mapq), loci in by_sample_mapq_locus.items():
        control_median = sample_control_median[(sample_id, mapq)]
        for locus_name, row in loci.items():
            if locus_name.startswith(CONTROL_PREFIX):
                continue
            ratio = normalized_ratio(row["median"], control_median) if control_median > 0 else float("inf")
            unnormalized.append(
                {
                    "sample_id": sample_id,
                    "locus": locus_name,
                    "mapq": mapq,
                    "locus_median_depth": row["median"],
                    "control_median_depth": round(control_median, 3),
                    "breadth": row["breadth"],
                    "normalized_ratio": round(ratio, 4),
                }
            )

    ratios_by_locus_mapq: dict[tuple[str, int], list[float]] = {}
    for entry in unnormalized:
        ratios_by_locus_mapq.setdefault((entry["locus"], entry["mapq"]), []).append(entry["normalized_ratio"])

    cohort_stats = {
        key: (statistics.median(values), statistics.median(abs(v - statistics.median(values)) for v in values))
        for key, values in ratios_by_locus_mapq.items()
    }

    for entry in unnormalized:
        cohort_median, cohort_mad = cohort_stats[(entry["locus"], entry["mapq"])]
        entry["cohort_median_ratio"] = round(cohort_median, 4)
        entry["cohort_mad_ratio"] = round(cohort_mad, 4)
        entry["interpretation"] = interpretation_flag(
            entry["normalized_ratio"], cohort_median, cohort_mad, entry["control_median_depth"], minimum_control_depth
        )

    return sorted(unnormalized, key=lambda item: (item["locus"], item["mapq"], item["sample_id"]))


def locus_range_rows(rows: Sequence[dict], mapq: int = 20) -> list[dict]:
    """Per-locus min/median/max at one MAPQ pass, independent of z-score flags."""
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        if row["mapq"] != mapq:
            continue
        grouped.setdefault(row["locus"], []).append(row)
    summary = []
    for locus_name in sorted(grouped):
        values = grouped[locus_name]
        ratios = [item["normalized_ratio"] for item in values]
        interpretations = [item["interpretation"] for item in values]
        summary.append(
            {
                "locus": locus_name,
                "minimum": min(ratios),
                "median": statistics.median(ratios),
                "maximum": max(ratios),
                "n_zero": sum(1 for ratio in ratios if ratio == 0),
                "n_modified_z": sum(
                    1 for flag in interpretations if flag in {MODIFIED_Z_ELEVATED, MODIFIED_Z_REDUCED}
                ),
                "n_unrepresented_exception": sum(1 for flag in interpretations if flag == UNREPRESENTED_EXCEPTION),
            }
        )
    return summary


def write_matrix_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "locus",
        "mapq",
        "sample_id",
        "locus_median_depth",
        "control_median_depth",
        "breadth",
        "normalized_ratio",
        "cohort_median_ratio",
        "cohort_mad_ratio",
        "interpretation",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_report(rows: Sequence[dict], sample_ids: Sequence[str]) -> str:
    mapq20 = [row for row in rows if row["mapq"] == 20]
    modified_z = [row for row in mapq20 if row["interpretation"] in {MODIFIED_Z_ELEVATED, MODIFIED_Z_REDUCED}]
    unrepresented = [row for row in mapq20 if row["interpretation"] == UNREPRESENTED_EXCEPTION]
    insufficient = [row for row in mapq20 if row["interpretation"] == INSUFFICIENT_CONTROL]
    ranges = locus_range_rows(rows, mapq=20)
    loci = sorted({row["locus"] for row in mapq20})
    lines = [
        "# ARS-UCD2.0 Normalized-Depth / CNV-Signal Report",
        "",
        f"{len(sample_ids)} samples x {len(loci)} loci, MAPQ>=20 pass shown below "
        "(MAPQ 0/10/20 are all recorded in `normalized_depth_matrix.tsv`).",
        "",
        "**Ratios are depth relative to each sample's median across six distributed single-copy "
        "control genes (ALB, MSTN, PRNP, RPP30, TBP, TERT). A ratio is not a validated integer "
        "copy number: for a paralogous BoLA locus, reduced ratio may reflect real copy loss or may "
        "reflect reference-allele divergence from this reference haplotype (see "
        "`docs/HANDOFF_PLAN_PROGRESS.md` §2b for a concrete case). Locus-specific claims require "
        "independent support (assembly, read-multihit, or wet typing) before publication.**",
        "",
        "Interpretation split at MAPQ>=20:",
        "",
        f"- modified z-score > {OUTLIER_THRESHOLD}: {len(modified_z)} of {len(mapq20)} rows",
        f"- non-zero depth at a locus whose cohort median ratio is 0: {len(unrepresented)} of {len(mapq20)} rows",
        f"- insufficient control depth: {len(insufficient)} of {len(mapq20)} rows",
        "",
        "Locus ranges are the source of truth for heterogeneity (for example JSP.1 near-zero "
        "depth in some animals). A row can sit at the extreme of a wide range without crossing "
        "the modified-z threshold when the cohort MAD is large. Class III loci with a very small "
        "MAD can flag modest ratio differences that are not deletions.",
        "",
        "## Locus ranges (MAPQ>=20)",
        "",
        "| locus | minimum | median | maximum | n_zero | n_modified_z | n_unrepresented_exception |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in ranges:
        lines.append(
            f"| {row['locus']} | {row['minimum']} | {row['median']} | {row['maximum']} | "
            f"{row['n_zero']} | {row['n_modified_z']} | {row['n_unrepresented_exception']} |"
        )
    lines.extend(
        [
            "",
            "## Modified-z flagged rows (MAPQ>=20)",
            "",
            "| locus | sample_id | normalized_ratio | cohort_median_ratio | cohort_mad_ratio | control_median_depth | interpretation |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    if modified_z:
        for row in sorted(modified_z, key=lambda item: (item["locus"], item["sample_id"])):
            lines.append(
                f"| {row['locus']} | {row['sample_id']} | {row['normalized_ratio']} | "
                f"{row['cohort_median_ratio']} | {row['cohort_mad_ratio']} | "
                f"{row['control_median_depth']} | {row['interpretation']} |"
            )
    else:
        lines.append("| - | - | - | - | - | - | none flagged |")
    lines.extend(
        [
            "",
            "## Non-zero at otherwise unrepresented loci (MAPQ>=20)",
            "",
            "These rows are not copy-number gains. The cohort median ratio is 0, so the locus "
            "is unrepresented (or highly divergent) in this population on the linear reference. "
            "A non-zero sample is an exception to that absence.",
            "",
            "| locus | sample_id | normalized_ratio | locus_median_depth | breadth | interpretation |",
            "|---|---|---|---|---|---|",
        ]
    )
    if unrepresented:
        for row in sorted(unrepresented, key=lambda item: (item["locus"], item["sample_id"])):
            lines.append(
                f"| {row['locus']} | {row['sample_id']} | {row['normalized_ratio']} | "
                f"{row['locus_median_depth']} | {row['breadth']} | {row['interpretation']} |"
            )
    else:
        lines.append("| - | - | - | - | - | none |")
    lines.append("")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="03_depth_cnv stage directory")
    parser.add_argument("--registry", required=True, type=Path, help="00_input_audit/sample_registry.tsv")
    parser.add_argument("--output", required=True, type=Path, help="Markdown report destination")
    parser.add_argument("--minimum-control-depth", type=float, default=10.0)
    args = parser.parse_args(argv)

    with args.registry.open(encoding="utf-8", newline="") as handle:
        sample_ids = sorted(
            row["sample_id"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("analysis_status", "").strip() == "analytical" and row.get("sample_id", "").strip()
        )

    cohort_depth = load_cohort_depth(args.root, sample_ids)
    rows = build_normalized_matrix(cohort_depth, args.minimum_control_depth)

    write_matrix_tsv(args.root / "normalized_depth_matrix.tsv", rows)
    mapq20 = [row for row in rows if row["mapq"] == 20]
    modified_z = [row for row in mapq20 if row["interpretation"] in {MODIFIED_Z_ELEVATED, MODIFIED_Z_REDUCED}]
    unrepresented = [row for row in mapq20 if row["interpretation"] == UNREPRESENTED_EXCEPTION]
    insufficient = [row for row in mapq20 if row["interpretation"] == INSUFFICIENT_CONTROL]
    status_payload = {
        "sample_count": len(sample_ids),
        "locus_count": len({row["locus"] for row in rows}),
        "flagged_row_count_mapq20": len(modified_z) + len(unrepresented) + len(insufficient),
        "modified_z_flagged_row_count_mapq20": len(modified_z),
        "unrepresented_exception_count_mapq20": len(unrepresented),
        "insufficient_control_count_mapq20": len(insufficient),
        "status": "passed",
    }
    (args.root / "normalized_depth_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        args.root,
        {
            "input_hashes": {
                sample_id: sha256_file(args.root / sample_id / "depth_raw.tsv") for sample_id in sample_ids
            }
        },
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(rows, sample_ids), encoding="utf-8")
    print(f"normalized depth report written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
