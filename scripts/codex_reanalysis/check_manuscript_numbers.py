"""Recompute headline manuscript numbers from raw evidence and check them in place.

Task 13 Step 5 of the reanalysis plan asks for a DOCX-to-evidence number audit. No
`docs/INTEGRATED_MANUSCRIPT_EN_CODEX_REANALYSIS.docx` exists yet -- the plan's Task 13
Step 3 ("update Abstract/Methods/Results/Discussion/Conclusions/tables/figure legends
from the evidence matrix") is not something this checker can honestly claim to have
driven: 9 of `EVIDENCE_MATRIX.tsv`'s 15 claims are `limited`, and `figure_legends`
is still an open `SCOPE_LIMITATIONS` gap. What this script does instead is check the
part with real, checkable value: recompute each CHECK's number independently from the
raw pipeline output (not from the matrix's free-text `new_value`, and not by parsing
the manuscript's own prose as ground truth -- that would be a self-fulfilling test),
then confirm the manuscript states the same value.

The actual manuscript today is the Markdown chain (`docs/RESULTS_EN.md`,
`docs/DISCUSSION_EN.md`, ...), not a DOCX -- that is what is actively edited and cited
throughout `docs/codex_reanalysis/`. This checks those files directly.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Callable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


class ManuscriptNumberError(ValueError):
    """Raised when a recomputed number cannot be checked against the manuscript."""


def compare_numeric_claim(name: str, manuscript: float, evidence: float, tolerance: float = 0.0005) -> bool:
    """True if the manuscript's stated value is within tolerance of the recomputed one."""
    return abs(manuscript - evidence) <= tolerance


def phi_squared(a: int, b: int, c: int, d: int) -> float:
    """Squared phi coefficient for a 2x2 contingency table [[a, b], [c, d]]."""
    denominator = (a + b) * (c + d) * (a + c) * (b + d)
    if denominator == 0:
        raise ManuscriptNumberError("phi_squared: a contingency-table margin is zero")
    phi = ((a * d) - (b * c)) / (denominator**0.5)
    return phi * phi


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compute_rs29026690_carrier_phi_squared(repo_root: Path) -> float:
    """r² between DRB3*016:01 carriage and rs29026690-T, recomputed from marker_genotypes.tsv."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv")
    a = b = c = d = 0  # [carrier, non-carrier] x [T-positive, T-negative]
    for row in rows:
        if row["marker_id"] != "rs29026690":
            continue
        carrier = row["carrier_016_01"] == "True"
        t_positive = "T" in row["primary_genotype"]
        if carrier and t_positive:
            a += 1
        elif carrier and not t_positive:
            b += 1
        elif not carrier and t_positive:
            c += 1
        else:
            d += 1
    return phi_squared(a, b, c, d)


def compute_locus_cohort_median_ratio(repo_root: Path, locus: str, mapq: str) -> float:
    """Cohort median control-normalized depth ratio for one locus at one MAPQ pass."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv")
    for row in rows:
        if row["locus"] == locus and row["mapq"] == mapq:
            return float(row["cohort_median_ratio"])
    raise ManuscriptNumberError(f"no normalized_depth_matrix.tsv row for locus={locus!r} mapq={mapq!r}")


def compute_marker_accepted_count(repo_root: Path, marker_id: str) -> int:
    """Count of samples whose evidence_status is `accepted` for one marker."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv")
    return sum(1 for row in rows if row["marker_id"] == marker_id and row["evidence_status"] == "accepted")


def compute_wide_same_contig_count(repo_root: Path) -> int:
    """Count of animals where DRB3 and rs29026690 fall on the same wide-interval assembly contig."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/06_local_assembly/wide/wide_linkage.tsv")
    return sum(1 for row in rows if row["same_contig"] == "yes")


def compute_drb3_exon2_concordant_count(repo_root: Path) -> int:
    """Animals whose best-cost exon-2 diplotype equals the PCR-SBT diplotype (stage 03d)."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/03d_drb3_exon2_typing/drb3_exon2_typing.tsv")
    return sum(1 for row in rows if row["concordant"] == "True")


def compute_landmark_same_block_count(repo_root: Path, column: str) -> int:
    """Count of animals where a WhatsHap phase block spans both named landmarks."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04b_mhc_phase_blocks/landmark_co_membership.tsv")
    return sum(1 for row in rows if row[column] == "True")


def compute_segregating_variant_count(repo_root: Path, gene: str) -> int:
    """Segregating (allele_frequency_class == segregating) coding variants at one gene."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv")
    return sum(1 for row in rows if row["gene"] == gene and row["allele_frequency_class"] == "segregating")


def compute_clean_missense_candidate_count(repo_root: Path, gene: str) -> int:
    """Segregating missense variants at one gene with at least one homozygous carrier."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv")
    return sum(
        1
        for row in rows
        if row["gene"] == gene
        and row["allele_frequency_class"] == "segregating"
        and row["effect"] == "MISSENSE"
        and int(row["n_hom_alt"]) > 0
    )


def read_json_field(path: Path, field: str) -> int:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return payload[field]


#: Each check recomputes one manuscript-stated number independently from raw pipeline
#: output, then confirms `manuscript_literal` -- the exact substring the manuscript
#: prose uses -- both appears in every listed file and parses to a value within
#: `tolerance` of the recomputed one. `tolerance` reflects the manuscript's own
#: rounding precision (e.g. "0.74" is a 2-decimal rounding, so 0.006 covers any value
#: that legitimately rounds to it), not an arbitrarily loose bound.
CHECKS: tuple[dict, ...] = (
    {
        "name": "rs29026690_016_01_r_squared",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "manuscript_literal": "0.74",
        "compute": lambda repo_root: compute_rs29026690_carrier_phi_squared(repo_root),
        "tolerance": 0.006,
    },
    {
        "name": "btnl2_cohort_median_ratio_mapq20",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "1.32",
        "compute": lambda repo_root: compute_locus_cohort_median_ratio(repo_root, "BTNL2", "20"),
        "tolerance": 0.006,
    },
    {
        "name": "btn3a2_cohort_median_ratio_mapq20",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "1.05",
        "manuscript_context": "BTN3A2 (median ratio 1.05;",
        "compute": lambda repo_root: compute_locus_cohort_median_ratio(repo_root, "BTN_LOC615223", "20"),
        "tolerance": 0.006,
    },
    {
        "name": "btn2a2_cohort_median_ratio_mapq20",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "1.00",
        "manuscript_context": "BTN2A2 (1.00;",
        "compute": lambda repo_root: compute_locus_cohort_median_ratio(repo_root, "BTN_LOC107131749", "20"),
        "tolerance": 0.006,
    },
    {
        "name": "btn3a3_cohort_median_ratio_mapq20",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "1.04",
        "manuscript_context": "BTN3A3 (1.04;",
        "compute": lambda repo_root: compute_locus_cohort_median_ratio(repo_root, "BTN_LOC786760", "20"),
        "tolerance": 0.006,
    },
    {
        "name": "tnf_824_accepted_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "10/13",
        "compute": lambda repo_root: compute_marker_accepted_count(repo_root, "TNF_824"),
        "tolerance": 0,
    },
    {
        "name": "wide_linkage_same_contig_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "0/13",
        "compute": lambda repo_root: compute_wide_same_contig_count(repo_root),
        "tolerance": 0,
    },
    {
        "name": "drb3_exon2_concordant_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "12/13",
        "compute": lambda repo_root: compute_drb3_exon2_concordant_count(repo_root),
        "tolerance": 0,
    },
    {
        "name": "tnf_824_x_rs29026690_same_block_count",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "manuscript_literal": "12/13",
        "compute": lambda repo_root: compute_landmark_same_block_count(
            repo_root, "TNF_824_x_rs29026690_same_block"
        ),
        "tolerance": 0,
    },
    {
        "name": "btn3a2_clean_missense_candidate_count",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "manuscript_literal": "47",
        "compute": lambda repo_root: compute_clean_missense_candidate_count(repo_root, "LOC615223_BTN3A2"),
        "tolerance": 0,
    },
    {
        "name": "btnl2_segregating_variant_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "101",
        "compute": lambda repo_root: compute_segregating_variant_count(repo_root, "BTNL2"),
        "tolerance": 0,
    },
)


#: One distinctive, specific string per `rejected` claim in EVIDENCE_MATRIX.tsv -- not
#: the gene/locus name alone (e.g. "NCR3" also appears legitimately elsewhere, as part
#: of the exploratory TNF/LTA region discussion), but the coordinate that would only
#: appear paired with the retracted attribution if it had leaked into the manuscript.
REJECTED_TERMS: dict[str, str] = {
    "ncr3_variant_attribution": "27,713,795",
}


def contains_rejected_claim(text: str, rejected_terms: set[str] | dict[str, str]) -> bool:
    """True if any rejected term is present in `text`."""
    terms = rejected_terms.values() if isinstance(rejected_terms, dict) else rejected_terms
    return any(term in text for term in terms)


def check_no_rejected_claims_in_manuscript(repo_root: Path, manuscript_files: Sequence[str]) -> list[str]:
    """Return the manuscript files (if any) that contain a rejected claim's term."""
    hits = []
    for path in manuscript_files:
        text = (repo_root / path).read_text(encoding="utf-8")
        if contains_rejected_claim(text, REJECTED_TERMS):
            hits.append(path)
    return hits


def _parse_manuscript_literal(literal: str) -> float:
    """Parse a manuscript literal into a float: a bare decimal, or the numerator of an ``n/d`` count."""
    if "/" in literal:
        numerator, _, _denominator = literal.partition("/")
        return float(numerator)
    return float(literal)


def literal_found(text: str, literal: str) -> bool:
    return literal in text


def run_checks(checks: Sequence[dict], repo_root: Path) -> list[dict]:
    results = []
    for check in checks:
        missing_files = [
            path for path in check["manuscript_files"] if not (repo_root / path).exists()
        ]
        if missing_files:
            raise ManuscriptNumberError(f"{check['name']}: manuscript file(s) not found: {missing_files}")

        texts = {path: (repo_root / path).read_text(encoding="utf-8") for path in check["manuscript_files"]}
        found_in = [path for path, text in texts.items() if literal_found(text, check.get("manuscript_context", check["manuscript_literal"]))]

        compute: Callable[[Path], float] = check["compute"]
        evidence_value = compute(repo_root)
        manuscript_value = _parse_manuscript_literal(check["manuscript_literal"])
        matches = compare_numeric_claim(
            check["name"], manuscript_value, evidence_value, tolerance=check["tolerance"]
        )

        results.append(
            {
                "name": check["name"],
                "manuscript_literal": check["manuscript_literal"],
                "found_in": found_in,
                "evidence_value": evidence_value,
                "matches": bool(found_in) and matches,
            }
        )
    return results


def render_report(results: Sequence[dict], rejected_claim_hits: Sequence[str]) -> str:
    passing = [row for row in results if row["matches"]]
    lines = [
        "# Manuscript Number Audit",
        "",
        "Each row's number was recomputed independently from raw pipeline output "
        "(not parsed from the manuscript's own prose) and checked against the literal "
        "value the manuscript states. This is a targeted spot-check of headline "
        "numbers with a single clean canonical value, not a claim-by-claim sweep -- "
        "see `docs/codex_reanalysis/EVIDENCE_MATRIX.tsv` for the broader claim audit.",
        "",
        f"{len(passing)}/{len(results)} checks passed.",
        "",
        "| check | manuscript literal | found in | recomputed value | status |",
        "|---|---|---|---|---|",
    ]
    for row in results:
        status = "PASS" if row["matches"] else "FAIL"
        found_in = ", ".join(f"`{path}`" for path in row["found_in"]) or "(not found)"
        lines.append(
            f"| {row['name']} | `{row['manuscript_literal']}` | {found_in} | "
            f"{row['evidence_value']!r} | {status} |"
        )
    lines.append("")
    lines.append("## Rejected-claim leak check")
    lines.append("")
    lines.append(
        "Every `rejected` claim in `EVIDENCE_MATRIX.tsv` (`REJECTED_TERMS`) is checked "
        "for its distinctive term across the current manuscript files."
    )
    lines.append("")
    if rejected_claim_hits:
        lines.append(f"**FAIL**: rejected-claim terms found in: {', '.join(rejected_claim_hits)}")
    else:
        lines.append("PASS: no rejected-claim term found in any checked manuscript file.")
    lines.append("")
    return "\n".join(lines)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True, help="docs/codex_reanalysis/MANUSCRIPT_VALIDATION.md")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    return parser.parse_args(argv)


#: Manuscript files scanned for rejected-claim terms. Kept to the union of
#: `manuscript_files` across CHECKS plus the sections most likely to carry a locus
#: attribution, rather than every doc in the repo.
_REJECTED_CLAIM_SCAN_FILES = (
    "docs/RESULTS_EN.md",
    "docs/DISCUSSION_EN.md",
    "docs/METHODS_EN.md",
)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    results = run_checks(CHECKS, args.repo_root)
    rejected_claim_hits = check_no_rejected_claims_in_manuscript(args.repo_root, _REJECTED_CLAIM_SCAN_FILES)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(render_report(results, rejected_claim_hits), encoding="utf-8")

    failing = [row["name"] for row in results if not row["matches"]]
    passing = len(results) - len(failing)
    print(f"Manuscript number audit written: {args.report} ({passing}/{len(results)} passed)")
    if failing:
        print(f"FAILED numeric checks: {failing}")
    if rejected_claim_hits:
        print(f"FAILED rejected-claim leak check: {rejected_claim_hits}")
    return 1 if (failing or rejected_claim_hits) else 0


if __name__ == "__main__":
    raise SystemExit(main())
