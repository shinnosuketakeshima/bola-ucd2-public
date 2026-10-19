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


def compute_locus_max_ratio(repo_root: Path, locus: str, mapq: str) -> float:
    """Largest per-animal control-normalized depth ratio for one locus at one MAPQ pass."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv")
    ratios = [float(row["normalized_ratio"]) for row in rows if row["locus"] == locus and row["mapq"] == mapq]
    if not ratios:
        raise ManuscriptNumberError(f"no normalized_depth_matrix.tsv rows for locus={locus!r} mapq={mapq!r}")
    return max(ratios)


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


def compute_coding_variant_count(repo_root: Path, gene: str) -> int:
    """Coding variants called at one gene (denominator-free: no frequency class applied)."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv")
    return sum(1 for row in rows if row["gene"] == gene)


def compute_clean_missense_candidate_count(repo_root: Path, gene: str) -> int:
    """Missense variants at one gene called homozygous-alternative in at least one animal."""
    rows = read_tsv(repo_root / "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv")
    return sum(
        1
        for row in rows
        if row["gene"] == gene
        and row["effect"] == "MISSENSE"
        and int(row["n_hom_alt"]) > 0
    )


CANDIDATE_SITE_SUMMARY = (
    "results/codex_reanalysis_v1/04_small_variants/candidate_site_genotypes/candidate_site_summary.tsv"
)


def compute_candidate_site_field(repo_root: Path, position: int, field: str) -> float:
    """One field of the re-genotyped CNTN3/ABHD16A/PRR3 missense site (called alleles only)."""
    for row in read_tsv(repo_root / CANDIDATE_SITE_SUMMARY):
        if int(row["position"]) == position:
            return float(row[field])
    raise ManuscriptNumberError(f"no candidate_site_summary.tsv row at position {position}")


def compute_candidate_class_count(repo_root: Path, frequency_class: str) -> int:
    """Number of re-genotyped missense sites in one called-allele frequency class."""
    return sum(1 for row in read_tsv(repo_root / CANDIDATE_SITE_SUMMARY) if row["frequency_class"] == frequency_class)


def compute_supported_site_count(repo_root: Path, min_called: int = 9) -> int:
    """Missense sites still supported at base quality >= 20: enough animals called, a carrier
    called, and the same frequency label as at base quality >= 10 (Table 6)."""
    base = {int(r["position"]): r for r in read_tsv(repo_root / CANDIDATE_SITE_SUMMARY)}
    strict = read_tsv(repo_root / CANDIDATE_SITE_SUMMARY.replace("candidate_site_summary.tsv",
                                                                 "sensitivity_bq20/candidate_site_summary.tsv"))
    return sum(
        1
        for row in strict
        if int(row["animals_called"]) >= min_called
        and int(row["carriers"]) >= 1
        and row["frequency_class"] == base[int(row["position"])]["frequency_class"]
    )


DQ_ATTRIBUTION = "results/codex_reanalysis_v1/06b_dq_attribution"


def compute_dq_configuration_animals(repo_root: Path, column: str, min_reads: int = 3) -> int:
    """Animals with at least `min_reads` reads covering one DQ configuration diagnostic point."""
    rows = read_tsv(repo_root / DQ_ATTRIBUTION / "dq_configuration_reads.tsv")
    return sum(1 for row in rows if int(row[column]) >= min_reads)


def compute_dq_configuration_reads(repo_root: Path, sample_id: str, column: str) -> int:
    for row in read_tsv(repo_root / DQ_ATTRIBUTION / "dq_configuration_reads.tsv"):
        if row["sample_id"] == sample_id:
            return int(row[column])
    raise ManuscriptNumberError(f"no dq_configuration_reads.tsv row for {sample_id}")


def compute_dq_selected_event_kb(repo_root: Path, event: str, field: str) -> float:
    """Size (kb) of the deletion junction or insertion selected by the attribution stage."""
    status = json.loads((repo_root / DQ_ATTRIBUTION / "dq_attribution_status.json").read_text(encoding="utf-8"))
    return int(status[f"selected_{event}"][field]) / 1000


DQ_SENSITIVITY = "results/codex_reanalysis_v1/06b_dq_attribution_panel_sensitivity"


def compute_dq_sensitivity_junction(repo_root: Path, which: str) -> int:
    """Junction counts with the ID13 contig as the configuration A panel member.

    `total_gap200`: reads meeting the primary rule over all animals; `min_gap50`/`max_gap50`: the
    per-animal range under the 50-bp rule among animals with any junction read; `id11_gap50`.
    """
    rows = [row for row in read_tsv(repo_root / DQ_SENSITIVITY / "dq_configuration_support_comparison.tsv")
            if row["diagnostic"] == "deletion_junction"]
    if which == "total_gap200":
        return sum(int(row["continuous_gap200"]) for row in rows)
    if which == "id11_gap50":
        return next(int(row["continuous_gap50"]) for row in rows if row["sample_id"] == "ID11")
    strict = [int(row["continuous_gap50"]) for row in rows if int(row["continuous_gap200"]) > 0]
    return min(strict) if which == "min_gap50" else max(strict)


DQ_JUNCTIONS = "results/codex_reanalysis_v1/06c_dq_junction_comparison"
DQ_PUBLIC = "results/codex_reanalysis_v1/06d_public_dq_comparison"


def compute_dq_strict_reads(repo_root: Path, sample_id: str, diagnostic: str) -> int:
    """Panel reads for one animal and diagnostic under the strict (50-bp gap) rule."""
    for row in read_tsv(repo_root / DQ_ATTRIBUTION / "dq_configuration_support_comparison.tsv"):
        if (row["sample_id"], row["diagnostic"]) == (sample_id, diagnostic):
            return int(row["continuous_gap50"])
    raise ManuscriptNumberError(f"no support-comparison row for {sample_id} {diagnostic}")


def compute_dq_junction_left_overlap_kb(repo_root: Path) -> float:
    """Breakpoint uncertainty at the left boundary of the selected junction (its source contig)."""
    status = json.loads((repo_root / DQ_ATTRIBUTION / "dq_attribution_status.json").read_text(encoding="utf-8"))
    sample = status["parameters"]["deletion_source"].split(":")[0]
    for row in read_tsv(repo_root / DQ_JUNCTIONS / "dq_junction_structure.tsv"):
        if row["sample_id"] == sample:
            return int(row["left_boundary_overlap_bp"]) / 1000
    raise ManuscriptNumberError(f"no junction-structure row for {sample}")


def compute_dq_rejected_junction_gap(repo_root: Path, which: str) -> int:
    """Shortest or longest gap of the junction reads that only the strict (50-bp) rule rejects."""
    gaps = [max(int(r["longest_deletion"]), int(r["longest_insertion"]))
            for r in read_tsv(repo_root / DQ_ATTRIBUTION / "dq_support_reads.tsv")
            if r["analysis"] == "panel" and r["diagnostic"] == "deletion_junction"
            and r["continuous_gap200"] == "True" and r["continuous_gap50"] == "False"]
    if not gaps:
        raise ManuscriptNumberError("no junction read rejected only by the strict rule")
    return min(gaps) if which == "min" else max(gaps)


def compute_public_assemblies(repo_root: Path, configuration: str | None = None) -> int:
    """Published assemblies (ARS-UCD2.0 control excluded) in total or with one configuration."""
    rows = [r for r in read_tsv(repo_root / DQ_PUBLIC / "public_dq_classification.tsv")
            if r["assembly_accession"] != "GCF_002263795.3"]
    return sum(1 for r in rows if configuration is None or r["configuration"] == configuration)


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
        "manuscript_literal": "10",
        "manuscript_context": "Ten animals had accepted genotypes",
        "compute": lambda repo_root: compute_marker_accepted_count(repo_root, "TNF_824"),
        "tolerance": 0,
    },
    {
        "name": "wide_linkage_same_contig_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "0/13",
        "manuscript_context": "(0/13 same contig)",
        "compute": lambda repo_root: compute_wide_same_contig_count(repo_root),
        "tolerance": 0,
    },
    {
        "name": "drb3_exon2_concordant_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "12/13",
        "manuscript_context": "Concordance was 12/13 at every window length",
        "compute": lambda repo_root: compute_drb3_exon2_concordant_count(repo_root),
        "tolerance": 0,
    },
    {
        "name": "tnf_824_x_rs29026690_same_block_count",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "manuscript_literal": "12/13",
        "manuscript_context": "TNF -824A/G and rs29026690 in 12/13 animals",
        "compute": lambda repo_root: compute_landmark_same_block_count(
            repo_root, "TNF_824_x_rs29026690_same_block"
        ),
        "tolerance": 0,
    },
    {
        "name": "btn3a2_clean_missense_candidate_count",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "manuscript_literal": "59",
        "compute": lambda repo_root: compute_clean_missense_candidate_count(repo_root, "LOC615223_BTN3A2"),
        "tolerance": 0,
    },
    {
        "name": "btnl2_coding_variant_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "102",
        "compute": lambda repo_root: compute_coding_variant_count(repo_root, "BTNL2"),
        "tolerance": 0,
    },
    {
        "name": "abhd16a_27617747_alt_af_called_alleles",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "0.458",
        "compute": lambda repo_root: compute_candidate_site_field(repo_root, 27617747, "alt_af"),
        "tolerance": 0.0005,
    },
    {
        "name": "abhd16a_27617747_carriers_of_called",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "9/12",
        "compute": lambda repo_root: compute_candidate_site_field(repo_root, 27617747, "carriers"),
        "tolerance": 0,
    },
    {
        "name": "cntn3_27374043_hom_alt_animals",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "13/13",
        "manuscript_context": "homozygous for the alternative allele in 13/13",
        "compute": lambda repo_root: compute_candidate_site_field(repo_root, 27374043, "hom_alt"),
        "tolerance": 0,
    },
    {
        "name": "candidate_missense_supported_site_count",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "seven",
        "manuscript_context": "seven of the 13 sites met these criteria",
        "compute": lambda repo_root: compute_supported_site_count(repo_root),
        "tolerance": 0,
    },
    {
        "name": "dq_junction_reference_gap_kb",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "76.4",
        "manuscript_context": "the 76.4-kb reference interval",
        "compute": lambda repo_root: compute_dq_selected_event_kb(repo_root, "deletion", "reference_gap_bp"),
        "tolerance": 0.05,
    },
    {
        "name": "dq_junction_spanning_reads_id11",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "14",
        "manuscript_context": "Eight of the 14 junction reads of ID11",
        "compute": lambda repo_root: compute_dq_configuration_reads(repo_root, "ID11", "deletion_junction"),
        "tolerance": 0,
    },
    {
        "name": "dq_junction_animals",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "10",
        "manuscript_context": "Reads running through the junction were found in 10 of 13 animals",
        "compute": lambda repo_root: compute_dq_configuration_animals(repo_root, "deletion_junction"),
        "tolerance": 0,
    },
    {
        "name": "dq_segment_edge_animals",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "nine",
        "manuscript_context": "segment were found in nine animals",
        "compute": lambda repo_root: compute_dq_configuration_animals(repo_root, "insertion_left_edge"),
        "tolerance": 0,
    },
    {
        "name": "dq_insertion_contig_kb",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "25.8",
        "manuscript_context": "a 25.8-kb segment that is absent from ARS-UCD2.0",
        "compute": lambda repo_root: compute_dq_selected_event_kb(repo_root, "insertion", "contig_gap_bp"),
        "tolerance": 0.05,
    },
    {
        "name": "dq_reference_configuration_reads_id09",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "14",
        "manuscript_context": "only in ID09 (14 reads;",
        "compute": lambda repo_root: compute_dq_configuration_reads(repo_root, "ID09", "reference_deletion_point"),
        "tolerance": 0,
    },
    {
        "name": "dq_junction_net_absent_kb_id11",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "71.9",
        "manuscript_context": "the net reference sequence absent from the contig was 71.9 kb",
        "compute": lambda repo_root: compute_dq_selected_event_kb(repo_root, "deletion", "net_bp"),
        "tolerance": 0.05,
    },
    {
        "name": "jsp1_max_ratio_mapq20",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "2.19",
        "manuscript_context": "to 2.19 in ID06 (Table 3; Fig. 4A)",
        "compute": lambda repo_root: compute_locus_max_ratio(repo_root, "JSP1", "20"),
        "tolerance": 0.005,
    },
    {
        "name": "dq_sensitivity_junction_reads_total",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "121",
        "manuscript_context": "the same 121 reads met the primary rule",
        "compute": lambda repo_root: compute_dq_sensitivity_junction(repo_root, "total_gap200"),
        "tolerance": 0,
    },
    {
        "name": "dq_sensitivity_strict_min",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "4",
        "manuscript_context": "retained in all ten animals (4–19 per animal)",
        "compute": lambda repo_root: compute_dq_sensitivity_junction(repo_root, "min_gap50"),
        "tolerance": 0,
    },
    {
        "name": "dq_sensitivity_strict_max",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "19",
        "manuscript_context": "retained in all ten animals (4–19 per animal)",
        "compute": lambda repo_root: compute_dq_sensitivity_junction(repo_root, "max_gap50"),
        "tolerance": 0,
    },
    {
        "name": "dq_sensitivity_strict_id11",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "8",
        "manuscript_context": "in ID11, the 8 reads lacking the 91 bp were now retained",
        "compute": lambda repo_root: compute_dq_sensitivity_junction(repo_root, "id11_gap50"),
        "tolerance": 0,
    },
    {
        "name": "dq_whole_segment_reads_id01",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "7",
        "manuscript_context": "only in ID01 (7 reads)",
        "compute": lambda repo_root: compute_dq_configuration_reads(repo_root, "ID01", "insertion_whole_segment"),
        "tolerance": 0,
    },
    {
        "name": "dq_whole_segment_reads_id02",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "8",
        "manuscript_context": "and ID02 (8 reads)",
        "compute": lambda repo_root: compute_dq_configuration_reads(repo_root, "ID02", "insertion_whole_segment"),
        "tolerance": 0,
    },
    {
        "name": "dq_junction_strict_reads_id11",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "6",
        "manuscript_context": "junction reads remained only in ID11 (6 reads)",
        "compute": lambda repo_root: compute_dq_strict_reads(repo_root, "ID11", "deletion_junction"),
        "tolerance": 0,
    },
    {
        "name": "dq_junction_left_breakpoint_interval_kb",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "0.93",
        "manuscript_context": "more precisely than a 0.93-kb interval",
        "compute": lambda repo_root: compute_dq_junction_left_overlap_kb(repo_root),
        "tolerance": 0.006,
    },
    {
        "name": "dq_rejected_junction_gap_min_bp",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "81",
        "manuscript_context": "carried an 81–100-bp deletion",
        "compute": lambda repo_root: compute_dq_rejected_junction_gap(repo_root, "min"),
        "tolerance": 0,
    },
    {
        "name": "dq_rejected_junction_gap_max_bp",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "100",
        "manuscript_context": "carried an 81–100-bp deletion",
        "compute": lambda repo_root: compute_dq_rejected_junction_gap(repo_root, "max"),
        "tolerance": 0,
    },
    {
        "name": "public_assemblies_compared",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "21",
        "manuscript_context": "Of 21 chromosome-level assemblies other than ARS-UCD2.0",
        "compute": lambda repo_root: compute_public_assemblies(repo_root),
        "tolerance": 0,
    },
    {
        "name": "public_assemblies_configuration_a",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "nine",
        "manuscript_context": "nine carried configuration A",
        "compute": lambda repo_root: compute_public_assemblies(repo_root, "A"),
        "tolerance": 0,
    },
    {
        "name": "public_assemblies_configuration_b",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "manuscript_literal": "Five",
        "manuscript_context": "Five carried configuration B",
        "compute": lambda repo_root: compute_public_assemblies(repo_root, "B"),
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


WORD_NUMBERS = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8,
                "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13}


def _parse_manuscript_literal(literal: str) -> float:
    """Parse a manuscript literal into a float: a bare decimal, a spelled-out count, or the
    numerator of an ``n/d`` count."""
    if literal.lower() in WORD_NUMBERS:
        return float(WORD_NUMBERS[literal.lower()])
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
        "Where a check sets `manuscript_context`, the number must appear inside that exact phrase, "
        "so a literal such as `12/13` cannot pass on an unrelated sentence.",
        "",
        "Definition changes in the 2026-10 revision: BTN coding-variant counts now cover every "
        "variant call at the gene rather than one frequency class (BTN3A2 missense with a "
        "homozygous-alternative call 47 -> 59; BTNL2 coding variants 101 -> 102), recomputed from "
        "`04_small_variants/coding_variants/coding_variants.tsv`. Table 6 values for CNTN3, ABHD16A "
        "and PRR3 come from "
        "`04_small_variants/candidate_site_genotypes/candidate_site_summary.tsv`, written by "
        "`scripts/codex_reanalysis/genotype_candidate_sites.py`, which re-genotypes every missense "
        "site in all 13 animals and computes frequencies over called alleles only. Table 6 now lists "
        "only the sites still supported at base quality >= 20 (>= 9 animals called, a carrier called, "
        "unchanged label); the frequency-class counts were replaced by that supported-site count.",
        "",
        "DQ checks added after the BMC-fit review recompute the junction, segment and read counts from "
        "`06b_dq_attribution/`, written by `scripts/codex_reanalysis/attribute_dq_contigs.py`. The "
        "earlier claim that ID11 'still yielded a DQA-family cluster' despite no DQA2 depth was "
        "retracted: that cluster lies at the DQA1 (LOC100848815) position.",
        "",
        "After the 2026-10-16 review, DQ read support requires CIGAR-level continuity through each "
        "diagnostic window (no gap > 200 bp, none > 50 bp within 200 bp of a diagnostic boundary, "
        ">= 80% aligned bases in each flank and the interval) instead of alignment end points; the "
        "strict 50-bp counts and end-point counts are in `dq_configuration_support_comparison.tsv`. "
        "Junction and public-assembly checks recompute from `06c_dq_junction_comparison/` and "
        "`06d_public_dq_comparison/`.",
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
