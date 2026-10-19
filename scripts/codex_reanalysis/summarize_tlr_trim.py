"""Sex-aware TLR/TRIM depth summary on the ARS-UCD2.0 codex_reanalysis pipeline.

TLR7 and TLR8 are X-linked (`data/meta/tlr_gene_coordinates.bed`); a male animal is
hemizygous there, so its depth ratio against the cohort's autosomal-normalized control
set is expected near 0.5 of a female's, not a CNV. `normalize_depth.py`'s generic
cohort-median/MAD classifier does not know this -- it treats the whole TLR7/TLR8
column as one mixed-sex distribution and calls everything "consistent_with_cohort" by
construction (the mixed distribution's own spread absorbs the sex effect). This script
adds the layer the legacy `scripts/summarize_tlr_sexaware.py` provided: infer each
animal's sex from mean(TLR7, TLR8) depth ratio, then report autosomal (non-TLR7/TLR8)
TLR/TRIM loci flagged by the generic classifier as the real candidate list.

The sex-inference threshold is derived from this cohort's own data, not reused from
the legacy pipeline's fixed 0.70/0.85 cutoffs -- the two pipelines normalize against a
different control set (six distributed single-copy genes here vs. three fixed legacy
regions), so the absolute ratio scale is not directly comparable between them.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path
from typing import Sequence

MALE = "male_inferred"
FEMALE = "female_inferred"
UNCERTAIN = "uncertain"

X_LINKED_GENES = ("TLR7", "TLR8")
MAPQ = "20"


class TlrTrimSummaryError(ValueError):
    """Raised when the depth matrix cannot be summarized safely."""


def read_matrix_rows(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compute_sex_inference_ratios(rows: Sequence[dict[str, str]]) -> dict[str, float]:
    """Return {sample_id: mean(TLR7, TLR8) normalized_ratio at MAPQ 20}."""
    by_sample: dict[str, dict[str, float]] = {}
    for row in rows:
        if row["mapq"] != MAPQ or row["locus"] not in X_LINKED_GENES:
            continue
        by_sample.setdefault(row["sample_id"], {})[row["locus"]] = float(row["normalized_ratio"])
    ratios = {}
    for sample_id, values in by_sample.items():
        if set(values) != set(X_LINKED_GENES):
            raise TlrTrimSummaryError(f"{sample_id} is missing TLR7 or TLR8 at MAPQ {MAPQ}")
        ratios[sample_id] = statistics.mean(values.values())
    return ratios


def find_sex_inference_gap(ratios: Sequence[float]) -> tuple[float, float]:
    """Return (threshold, gap_width): the midpoint of the widest gap in sorted ratios.

    A wide, unambiguous gap is what justifies calling every animal confidently rather
    than carrying a fixed external threshold; a narrow gap is reported so the caller
    can decide whether "uncertain" is warranted (see `classify_sex`).
    """
    ordered = sorted(ratios)
    if len(ordered) < 2:
        raise TlrTrimSummaryError("need at least 2 samples to find a sex-inference gap")
    gaps = [(ordered[i + 1] - ordered[i], (ordered[i + 1] + ordered[i]) / 2) for i in range(len(ordered) - 1)]
    width, threshold = max(gaps)
    return threshold, width


def classify_sex(ratio: float, threshold: float, *, uncertainty_band: float) -> str:
    if ratio < threshold - uncertainty_band:
        return MALE
    if ratio > threshold + uncertainty_band:
        return FEMALE
    return UNCERTAIN


def autosomal_outliers(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """Rows for TLR/TRIM loci (excluding TLR7/TLR8) not classified consistent_with_cohort, at MAPQ 20."""
    outliers = []
    for row in rows:
        if row["mapq"] != MAPQ:
            continue
        locus = row["locus"]
        if not (locus.startswith("TLR") or locus.startswith("TRIM")):
            continue
        if locus in X_LINKED_GENES:
            continue
        if row["interpretation"] != "consistent_with_cohort":
            outliers.append(row)
    return outliers


def render_report(
    sex_ratios: dict[str, float],
    threshold: float,
    gap_width: float,
    uncertainty_band: float,
    outliers: Sequence[dict[str, str]],
) -> str:
    lines = [
        "# TLR/TRIM depth summary, sex-aware (ARS-UCD2.0)",
        "",
        "Source: `results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv` "
        "(MAPQ 20 rows), the same cohort-internal six-control-gene normalization used "
        "throughout `DEPTH_CNV_REPORT.md`. TLR7/TLR8 are X-linked; the generic "
        "cohort-median/MAD classifier does not know this and folds male hemizygosity "
        "into the column's own spread, so it never flags TLR7/TLR8 on its own -- sex "
        "must be inferred separately before those two loci can be interpreted.",
        "",
        "## Sex inference",
        "",
        f"Threshold: mean(TLR7, TLR8) normalized ratio at MAPQ 20, split at the widest "
        f"gap in this cohort's own sorted values ({threshold:.3f}, gap width "
        f"{gap_width:.3f}); a ±{uncertainty_band:.2f} band around the threshold is "
        "reported as uncertain rather than forced to a side. This threshold is derived "
        "from this pipeline's own control-normalized scale, not reused from the legacy "
        "TLR_CNV_REPORT.md's fixed 0.70/0.85 cutoffs, which were computed against a "
        "different (three-fixed-region) control set.",
        "",
        "**This table is not a substitute for `RESULTS_EN.md` Table 1's inferred_sex "
        "column, and does not override it.** Table 1's sex calls were validated "
        "against primary blood-sampling records (9/9 concordance for the animals with "
        "a recorded sex field) and that validation is untouched by this pass. One "
        "animal (ID02) disagrees between the two pipelines: Table 1 "
        "calls it male (legacy TLR7/TLR8 = 0.80/0.47, mean 0.63, primary-record-"
        "confirmed), this pipeline's gap-based threshold calls it female (TLR7/TLR8 = "
        "0.96/0.65, mean 0.81) -- its TLR7 ratio is elevated under both "
        "normalizations (0.80 under legacy, right at their own 0.85 female cutoff; "
        "0.96 here), so this is a pre-existing borderline signal common to both "
        "methods, not a new artifact. Given Table 1's call is anchored to a real "
        "record and this one is not, Table 1's male call stands; this report's sex "
        "column exists only to normalize the autosomal TLR/TRIM candidate list below, "
        "and this one disagreement is flagged rather than resolved.",
        "",
        "| sample_id | mean(TLR7,TLR8) ratio | inferred sex |",
        "|---|---|---|",
    ]
    for sample_id in sorted(sex_ratios, key=lambda s: sex_ratios[s]):
        ratio = sex_ratios[sample_id]
        sex = classify_sex(ratio, threshold, uncertainty_band=uncertainty_band)
        lines.append(f"| {sample_id} | {ratio:.4f} | {sex} |")

    male_count = sum(1 for r in sex_ratios.values() if classify_sex(r, threshold, uncertainty_band=uncertainty_band) == MALE)
    female_count = sum(1 for r in sex_ratios.values() if classify_sex(r, threshold, uncertainty_band=uncertainty_band) == FEMALE)
    uncertain_count = len(sex_ratios) - male_count - female_count
    lines += [
        "",
        f"male_inferred: {male_count}; female_inferred: {female_count}; uncertain: {uncertain_count}",
        "",
        "## Autosomal TLR/TRIM loci flagged relative to cohort (MAPQ 20, TLR7/TLR8 excluded)",
        "",
        "These are candidates flagged by the generic cohort-median/MAD classifier "
        "(`normalize_depth.py`), the same threshold logic applied to every locus in "
        "the panel -- a statistical outlier relative to this cohort's own spread, not "
        "a validated integer copy-number call (`docs/codex_reanalysis/DEPTH_CNV_REPORT.md` "
        "states this explicitly). This is a more sensitive, cohort-relative flag than "
        "the legacy pipeline's fixed 0.15/0.60/1.60 absolute thresholds, so the "
        "candidate counts below are not directly comparable to "
        "`results/tlr/TLR_CNV_REPORT.md` / `results/trim/TRIM_CNV_REPORT.md`.",
        "",
    ]
    if outliers:
        lines.append("| locus | sample_id | normalized_ratio | cohort_median_ratio | interpretation |")
        lines.append("|---|---|---|---|---|")
        for row in sorted(outliers, key=lambda r: (r["locus"], r["sample_id"])):
            lines.append(
                f"| {row['locus']} | {row['sample_id']} | {row['normalized_ratio']} | "
                f"{row['cohort_median_ratio']} | {row['interpretation']} |"
            )
    else:
        lines.append("None.")
    lines.append("")
    return "\n".join(lines)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--uncertainty-band", type=float, default=0.05)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    rows = read_matrix_rows(args.matrix)
    sex_ratios = compute_sex_inference_ratios(rows)
    threshold, gap_width = find_sex_inference_gap(list(sex_ratios.values()))
    outliers = autosomal_outliers(rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        render_report(sex_ratios, threshold, gap_width, args.uncertainty_band, outliers),
        encoding="utf-8",
    )
    print(f"TLR/TRIM sex-aware summary written: {args.output} ({len(outliers)} autosomal outliers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
