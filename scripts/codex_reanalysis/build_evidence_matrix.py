"""Claim-level evidence audit of the manuscript against ARS-UCD2.0 pipeline output.

One row per manuscript claim, linking its current wording to the pipeline artifact
that should support it and a classification (`accepted`/`limited`/`rejected`) per
`docs/superpowers/plans/2026-07-31-ucd2-paper-reanalysis.md` Task 12. The legacy
comparison baseline is `docs/INTEGRATED_MANUSCRIPT_EN.docx` (the plan's Global
Constraints name it as immutable), not `docs/CHATGPT_EVIDENCE_AUDIT.md` -- that file
documents a separate, independently built manuscript lineage with its own,
non-matching numbers, so diffing against it would report unrelated-lineage
differences as "drift" against nothing.

Anonymization boundary: manuscript claims are written with anonymized IDs
(ID01-ID13, `docs/ANONYMIZATION_KEY_INTERNAL.md`, gitignored); pipeline evidence is
keyed by real sample_id. `docs/codex_reanalysis/EVIDENCE_MATRIX.tsv` is tracked, so a
row may cite an anonymized ID or an evidence file/field, but never both the
anonymized ID and the real sample_id it decodes to -- doing so would republish the
gitignored key. `assert_no_identity_leak` enforces this before every write.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance

ACCEPTED = "accepted"
LIMITED = "limited"
REJECTED = "rejected"

#: Real sample_id prefixes and raw deliverable IDs used throughout results/data/.
#: Never legitimate in a tracked claim-to-evidence row: the manuscript side of every
#: row is anonymized (METHODS_EN.md:18 maps DELIVERY_14/15 to ID07/ID10, for example).
_REAL_SAMPLE_ID_RE = re.compile(r"\b(DELIVERY_|DELIVERY_|DELIVERY_|DELIVERY_\d+)\w*", re.IGNORECASE)


class EvidenceMatrixError(ValueError):
    """Raised when a claim row fails validation before it can be written."""


def classify_claim(
    source_exists: bool,
    reproducible: bool,
    independent_support: bool = False,
    limitation: bool = False,
) -> str:
    """Classify one claim as accepted, limited, or rejected.

    A claim is ``rejected`` if the new pipeline has no source for it or that
    source does not reproduce it -- ``independent_support`` cannot rescue a claim
    that fails either of those. Otherwise a declared ``limitation`` (a depth-gated
    marker call, an unconfirmed SV breakpoint, a claim explicitly hedged as
    hypothesis in the text) always yields ``limited``, even when independent
    support also exists: corroboration narrows a limitation, it does not erase it.
    A claim with a reproducing source and no declared limitation is ``accepted``.
    """
    if not source_exists or not reproducible:
        return REJECTED
    if limitation:
        return LIMITED
    return ACCEPTED


#: Sections of the manuscript this seed does NOT cover yet, with why. Named here so
#: an absent row reads as "not yet audited," never as "no finding" -- the same
#: principle pilot_gate.py applies by reporting an absent stage as `not_run` rather
#: than `passed`. Extend this dict (and CLAIMS) rather than treating the 10 seed
#: claims as a complete sweep of Abstract/Results/Discussion/Conclusions/tables/figures.
SCOPE_LIMITATIONS: dict[str, str] = {}

#: Curated claim registry. Seeded from docs/codex_reanalysis/EVIDENCE_REVIEW_2026-08-19.md's
#: hand-verified findings (main pass + addenda) rather than an automated sweep of the
#: manuscript files, so every row traces to analysis already checked against its cited
#: source. `legacy_excerpt` is verified against the docx baseline by
#: compare_legacy_claims.py; `None` means the claim has no legacy counterpart (new in
#: the ARS-UCD2.0 reanalysis) or its legacy source is not the docx (noted in `note`).
#: See SCOPE_LIMITATIONS above for what this seed does not cover yet.
CLAIMS: tuple[dict, ...] = (
    {
        "claim_id": "class_i_gene_coordinates",
        "section": "Results — structural map / Table S1",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/METHODS_EN.md"),
        "legacy_excerpt": "the evaluated classical class I region at approximately 28.49–28.72 Mb",
        "legacy_value": "class I region ~28.49-28.72 Mb (ARS-UCD1.2, corrected BED data/meta/bola_cnv_genes.v2.bed)",
        "new_evidence_paths": ("results/codex_reanalysis_v1/01_reference/coordinate_registry.tsv",),
        "new_value": "class I 28.49-28.72 Mb (ARS-UCD2.0 coordinate_registry.tsv, GFF-derived)",
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": False,
        "note": (
            "BTA23 is byte-identical between ARS-UCD1.2 and ARS-UCD2.0 (SHA-256, "
            "assembly_comparison/sequence_comparison.tsv), so the docx's already-corrected "
            "coordinates carry over unchanged; coordinate_registry.tsv independently "
            "reproduces them from the pinned GFF."
        ),
        "action": "keep as-is",
    },
    {
        "claim_id": "btnl2_depth_range",
        "section": "Results — butyrophilin depth",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": "BTNL2 ranged from 0.962 to 3.355 (median 1.361)",
        "legacy_value": "BTNL2 0.962-3.355, median 1.361 (3 fixed autosomal controls, fixed thresholds)",
        "new_evidence_paths": (
            "docs/codex_reanalysis/DEPTH_CNV_REPORT.md",
            "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv",
        ),
        "new_value": "BTNL2 0.9846-3.4667, median 1.3239 (6 distributed single-copy controls, MAD-based z)",
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": False,
        "note": (
            "Range/median shift ~2-3% given the different control-set methodology "
            "(disclosed in Methods' two-generations paragraph), not measurement drift; "
            "RESULTS_EN.md/PAPER_DRAFT.md Tables 3/4 were recomputed from "
            "normalized_depth_matrix.tsv to match exactly (EVIDENCE_REVIEW addendum, "
            "\"Numeric drift\")."
        ),
        "action": "already corrected",
    },
    {
        "claim_id": "class_i_window_depth_ranges",
        "section": "Results — class I depth",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": "NC1 ratios ranged from 0.67 to 1.69 (median 1.23), JSP.1 from 0.16 to 2.20 (median 0.89), and classI_BOLA from 0.57 to 2.39 (median 1.33)",
        "legacy_value": "NC1 0.67-1.69 (median 1.23); JSP.1 0.16-2.20 (median 0.89); classI_BOLA 0.57-2.39 (median 1.33)",
        "new_evidence_paths": (
            "docs/codex_reanalysis/DEPTH_CNV_REPORT.md",
            "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv",
        ),
        "new_value": "NC1 0.5882-1.8148 (median 1.2113); JSP1 0.1333-2.1852 (median 0.8193); classI_BoLA_A 0.5846-1.6154 (median 1.012), all MAPQ-pass breakdown in DEPTH_CNV_REPORT.md",
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": False,
        "note": (
            "Same loci, same qualitative heterogeneity pattern under the new control "
            "methodology; RESULTS_EN.md now cites per-animal MAPQ>=20 figures directly "
            "from DEPTH_CNV_REPORT.md rather than the docx's cohort range, so the two "
            "are not compared as the same statistic -- flagged here for a reader who "
            "wants the cohort-range form too."
        ),
        "action": "no manuscript action needed (cohort-range form is not written in the manuscript; documented here for reviewers who want it)",
    },
    {
        "claim_id": "dq_dr_no_carrier_pattern",
        "section": "Results — DQ/DR local assembly",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": "No consistent DQ copy-number pattern distinguished DRB3*016:01 carriers from non-carriers in this ascertained cohort",
        "legacy_value": "no consistent DQA/DQB recovery-count pattern by carrier status (Flye 2.9.3-b1797, docx Table 2B)",
        "new_evidence_paths": (
            "docs/codex_reanalysis/ASSEMBLY_DQ_DR_REPORT.md",
            "results/codex_reanalysis_v1/06_local_assembly/dq_dr/dq_dr_recovery.tsv",
        ),
        "new_value": "no consistent DQ cluster pattern distinguishes carriers from non-carriers (Flye 2.9.6-b1802, same NC_037350.1:25,400,000-26,000,000 interval)",
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": False,
        "note": "Same interval, same qualitative conclusion under the newer pinned Flye/pipeline.",
        "action": "keep as-is",
    },
    {
        "claim_id": "rs29026690_016_01_correlation",
        "section": "Results / Discussion — GWAS marker linkage",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": "The squared binary correlation was r² = 0.7347, reported as 0.74",
        "legacy_value": "r²=0.74 (6/6 carriers + 1/7 non-carrier T-positive); explicitly not population LD",
        "new_evidence_paths": (
            "docs/codex_reanalysis/MARKER_REPORT.md",
            "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv",
        ),
        "new_value": "same r²=0.74 genotype pattern reproduced; only 2/13 underlying rs29026690 calls meet the two-pileup/depth `accepted` bar, 11/13 `limited`",
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "The correlation itself is reproduced, but MARKER_REPORT.md shows most of "
            "the genotypes it is built from are depth-limited; Discussion now states "
            "the 2/13-accepted caveat next to both r²=0.74 statements "
            "(EVIDENCE_REVIEW, \"Discussion staleness\")."
        ),
        "action": "already corrected",
    },
    {
        "claim_id": "tnf_824_marker_genotyping",
        "section": "Results / Discussion — TNF -824A/G",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": "the TNF −824 genotype was unstable in the present R9 reads and cannot be linked to *016:01 here",
        "legacy_value": "TNF_824 genotype called unstable, not usable (no accepted/limited tiering reported)",
        "new_evidence_paths": ("docs/codex_reanalysis/MARKER_REPORT.md",),
        "new_value": "TNF_824 10/13 animals meet `accepted` evidence tier; no genotype pattern separates carriers from non-carriers",
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": False,
        "note": (
            "The two-pileup/depth-gated protocol resolves most TNF_824 genotypes to "
            "`accepted`, an improvement on the docx's blanket 'unstable'; Discussion "
            "still correctly declines to draw a carrier-association conclusion from it."
        ),
        "action": "keep as-is",
    },
    {
        "claim_id": "drb3_rs29026690_recombination_boundary",
        "section": "Results / Discussion — phase-block boundary",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": "A previously generated broad phase-block report used the superseded class I/class II coordinate compartments and was excluded; no recombination hotspot was inferred from it",
        "legacy_value": "prior phase-block/hotspot analysis explicitly excluded (used retracted coordinate compartments); no boundary claimed",
        "new_evidence_paths": (
            "docs/codex_reanalysis/MHC_PHASE_BLOCKS_REPORT.md",
            "results/codex_reanalysis_v1/04b_mhc_phase_blocks/all_boundary_hotspots_100kb.tsv",
            "results/codex_reanalysis_v1/04b_mhc_phase_blocks/landmark_co_membership.tsv",
            "docs/codex_reanalysis/ASSEMBLY_WIDE_REPORT.md",
            "results/codex_reanalysis_v1/06_local_assembly/wide/wide_linkage.tsv",
        ),
        "new_value": (
            "ARS-UCD2.0 WhatsHap re-run on target BAMs: strongest boundary window "
            "26.20-26.30 Mb with 21 breaks (legacy 22 at the same 100-kb bin); secondary "
            "peaks 25.60-25.70 (13) and 26.80-26.90 (11). DRB3 x rs29026690 same block 0/13; "
            "DRB3 x rs17872126 0/13; rs29026690 x rs17872126 7/13 (legacy 4/13). Corrected "
            "compartments label the peak as proximal BTN / DR-BTN interval (not the retracted "
            "between_class_I__class_II_DQ label). Wide-interval assembly independently confirms "
            "non-bridging in all 13 animals."
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "Ported this pass: `scripts/codex_reanalysis/mhc_phase_blocks.py` re-called SNPs "
            "and re-phased on UCD2.0 target BAMs (WhatsHap 2.2), with corrected compartment "
            "windows. The ~26.25 Mb location is no longer single-source legacy phasing. Still "
            "`limited` (not `accepted`): phase-block gaps are discontinuities in read-backed "
            "phasing, not a measured meiotic recombination rate. Fig. 7 HTML haploblock data "
            "and chatgpt SVG composites were regenerated from UCD2.0 outputs in item 6."
        ),
        "action": (
            "corrected this pass (ported WhatsHap to ARS-UCD2.0; METHODS/RESULTS/PAPER_DRAFT/"
            "two-generations paragraph updated; figure data regenerated item 6)"
        ),
    },
    {
        "claim_id": "extended_haplotype_mechanism_hypothesis",
        "section": "Discussion — ancestral extended-haplotype mechanism",
        "manuscript_files": ("docs/DISCUSSION_EN.md",),
        "legacy_excerpt": "The figure must not label a recombination hotspot, population LD block, or demonstrated extended haplotype",
        "legacy_value": "extended-haplotype framing explicitly forbidden in the docx (no supporting evidence at the time)",
        "new_evidence_paths": ("docs/DISCUSSION_EN.md",),
        "new_value": (
            "framed as a candidate mechanism analogous to human MHC ancestral haplotypes, "
            "explicitly hedged: \"this discussion remains at the level of hypothesis\", with the "
            "Site-A relatedness confound disclosed"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "Distinct from drb3_rs29026690_recombination_boundary: the boundary finding itself "
            "is well-supported (see that claim), but the *evolutionary mechanism* explaining why "
            "it persists (NAHR, closed herdbook, ancestral haplotype) is explicitly presented as "
            "hypothesis, not a result -- correct to keep as `limited`, not to promote to accepted."
        ),
        "action": "keep as-is (already correctly hedged as hypothesis)",
    },
    {
        "claim_id": "sv_breakpoint_evidence_scope",
        "section": "Results / Discussion — BTN sub-locus structural claims",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/METHODS_EN.md"),
        "legacy_excerpt": None,
        "legacy_value": "not present in legacy baseline (docx does not perform or discuss SV calling)",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/05_sv/sv_evidence_status.json",
            "docs/codex_reanalysis/SV_REPORT.md",
        ),
        "new_value": "0/1252 breakpoint groups meet the `supported` bar cohort-wide (Sniffles2 2.8.0, no second caller yet)",
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "New in this reanalysis generation. Sub-locus depth-shape claims (BTNL2/BTN2A1-like "
            "windowed profiles) already cite sv_evidence_status.json to flag themselves as "
            "depth-based, SV-caller-unconfirmed observations, which is the correct manuscript "
            "action for a `limited` claim."
        ),
        "action": "keep as-is (already correctly flagged as SV-unconfirmed)",
    },
    {
        "claim_id": "ncr3_variant_attribution",
        "section": "Not in manuscript — exploratory memo only",
        "manuscript_files": (),
        "legacy_excerpt": None,
        "legacy_value": "docs/CLASSIII_016_TAG_SNP_EXPLORATION.md §5 originally attributed chr23:27,713,795 to NCR3 (receptor-ligand narrative with BAG6)",
        "new_evidence_paths": ("results/codex_reanalysis_v1/01_reference/coordinate_registry.tsv",),
        "new_value": "pinned ARS-UCD2.0 GFF: NCR3 ends at 27,706,075; the variant is ~7.7 kb downstream in LTB, not NCR3",
        "source_exists": True,
        "reproducible": False,
        "independent_support": False,
        "limitation": False,
        "note": (
            "Old-coordinate gene assignment not reproduced by UCD2.0 targets, per Task 12's "
            "explicit reject criterion. Caught and corrected in the source memo before it ever "
            "entered manuscript prose (EVIDENCE_REVIEW addendum); listed here as the matrix's "
            "required rejected-claim example, not as an open manuscript risk."
        ),
        "action": "no manuscript action needed (never cited); memo already corrected in place",
    },
    {
        "claim_id": "btn_coding_variants",
        "section": "Results — BTN coding-sequence screen",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": "The coding screen detected many candidate differences at BTN3A2 and fewer at BTN2A2 and BTN3A3",
        "legacy_value": "qualitative only: most candidate coding differences at BTN3A2, fewer at BTN2A2/BTN3A3",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv",
            "docs/codex_reanalysis/CODING_VARIANTS_REPORT.md",
        ),
        "new_value": "BTN3A2 216 total (14 breed_divergent, 202 segregating: 119 missense); BTN2A2 18 segregating (9 missense); BTN3A3 18 segregating (12 missense)",
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "Same qualitative ranking as the docx (BTN3A2 >> BTN2A2/BTN3A3). Ported to the "
            "ARS-UCD2.0 codex_reanalysis pipeline this pass (`call_coding_variants.py`, "
            "13/13 samples, `02_mapping` target BAMs); the exact per-gene counts match the "
            "ARS-UCD1.2-generation screen (`results/btn/btn_coding_variants.tsv`) once the "
            "new pipeline's finer breed_divergent/segregating split is collapsed back to a "
            "single 'candidate' tier -- independently confirmed via "
            "`validate_cds_models.py`: the legacy and GFF-derived single-canonical-transcript "
            "CDS models for these 5 BTN genes are exon-for-exon identical (translate to a "
            "clean single ORF, both sides). One additional BTN paralog, LOC786760/BTN3A3's "
            "sister locus is fine, but its own transcript in the pinned GFF does not translate "
            "cleanly (2 stop codons, only one annotated isoform exists) -- reading frame "
            "unresolved for that one locus, noted in METHODS_EN.md."
        ),
        "action": "already corrected (Methods and two-generations paragraph moved this out of the ARS-UCD1.2-generation list this pass)",
    },
    {
        "claim_id": "cntn3_abhd16a_prr3_coding_variants",
        "section": "Results — coding variation flanking the GWAS SNPs",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": None,
        "legacy_value": "not present in legacy baseline (docx does not screen CNTN3/ABHD16A/PRR3 coding sequence)",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/04_small_variants/coding_variants/coding_variants.tsv",
            "docs/codex_reanalysis/CODING_VARIANTS_REPORT.md",
        ),
        "new_value": (
            "CNTN3 5 missense (two at 13/13, one at 3/13, one breed-fixed 13/13 hom "
            "27,374,043 Ile>Val, one rare 2/13 27,414,294 Ser>Phe); ABHD16A 4 missense "
            "(one at 10/13, 27,617,747 Thr>Ala); PRR3 4 missense"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "The previously reported ARS-UCD1.2-generation numbers for these three genes "
            "were wrong, not merely stale, discovered while porting this stage to the "
            "codex_reanalysis pipeline: the legacy `data/meta/cntn3_cds.bed` and "
            "`bta23_gwas_cds.bed` each mix CDS exons from more than one annotated transcript "
            "isoform of the same gene, which `validate_cds_models.py` shows does not "
            "translate as a clean single open reading frame (CNTN3: translation doesn't "
            "start with Met; ABHD16A: same; PRR3: length not a multiple of 3) -- confirmed "
            "independently of any prose claim, by translation alone. Rebuilt each gene's CDS "
            "from a single canonical GFF transcript (CNTN3 XM_024983152.2, ABHD16A "
            "NM_001075605.1, PRR3 NM_001034554.1; each verified to pass the clean-ORF check) "
            "and re-ran `bcftools mpileup`/`call -mv` against the corrected restriction BED "
            "for all 13 animals. Previously reported: CNTN3 'Leu>Phe' (27,391,361, now "
            "Ser>Asn), CNTN3 'Glu>Lys, charge-changing' (27,414,294, now Ser>Phe -- not a "
            "charge change), CNTN3 'Met>Thr' (27,374,043, now Ile>Val), and an ABHD16A row "
            "at 27,599,354 that falls outside the correct transcript's CDS entirely and does "
            "not exist as a real coding variant. ABHD16A's other four rows and all four PRR3 "
            "rows matched exactly once the exon model was corrected -- an independent "
            "confirmation that those specific calls were right even though the underlying "
            "legacy CDS file had the same isoform-mixing defect; it happened not to change "
            "their frame. The corrected screen also found CNTN3 has 5 real missense variants "
            "against the legacy count of 3 (a completeness gain, not just a correction), and "
            "the previously highlighted 'charge-changing functional candidate' framing no "
            "longer applies to any CNTN3 variant. Confirmed with the user before publishing "
            "corrected amino-acid calls, given this touches a named functional-candidate claim."
        ),
        "action": "corrected this pass (RESULTS_EN.md/DISCUSSION_EN.md/FULL_MANUSCRIPT_EN.md/PAPER_DRAFT.md prose and Table 6; user-confirmed)",
    },
    {
        "claim_id": "tlr_trim_depth",
        "section": "Results — TLR/TRIM loci (Table S2)",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": None,
        "legacy_value": "not present in legacy baseline (docx does not discuss TLR/TRIM)",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv",
            "docs/codex_reanalysis/TLR_TRIM_CNV_REPORT.md",
        ),
        "new_value": (
            "TLR2-10 (TLR1 has no standalone ARS-UCD2.0 locus, not screened) all consistent with "
            "cohort at MAPQ 20 -- zero autosomal TLR candidates, versus the legacy TLR3/TLR4 pair. "
            "12/18 TRIM-panel genes consistent in every animal; 6 flagged (TRIM10 x2 animals, "
            "TRIM11, TRIM14, TRIM25, TRIM27, TRIM40 x1 each) -- a different, non-overlapping "
            "candidate set from the legacy TRIM21-only finding"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "Ported to the ARS-UCD2.0 codex_reanalysis pipeline this pass: added TLR2-10 and the "
            "18-gene TRIM panel to `configs/codex_reanalysis.targets.tsv` (TLR1 excluded -- it "
            "survives only as a `gene_synonym` on TLR6/TLR10 in the pinned annotation, no invented "
            "coordinate), which required a full 13-animal cohort remapping since these loci are "
            "scattered across 14 chromosomes outside the previously mapped BoLA-extended region. "
            "`compute_depth.py`/`normalize_depth.py` are fully generic over `targets.tsv`, so no "
            "new depth code was needed; wrote `summarize_tlr_trim.py` for the one TLR/TRIM-specific "
            "piece -- sex inference from mean(TLR7,TLR8) ratio, since the generic cohort-median/MAD "
            "classifier folds male X-hemizygosity into the mixed-sex column's own spread rather "
            "than flagging it. The new candidate set is not comparable to the legacy one: this "
            "pipeline's cohort-relative MAD statistic is more sensitive than the legacy fixed "
            "0.15/0.60/1.60 absolute thresholds, and the two flag entirely non-overlapping loci "
            "(TLR3/TLR4/TRIM21 here are all `consistent_with_cohort`; the 6 flagged TRIM genes "
            "here were not flagged by legacy). Cross-checked the new sex-inference against "
            "RESULTS_EN.md Table 1's inferred_sex column (validated against primary blood-sampling "
            "records, 9/9 concordance): 12/13 agree; ID02 disagrees (Table 1: male, "
            "primary-record-confirmed; this pipeline: female) -- its TLR7 ratio is elevated under "
            "both normalizations (legacy 0.80, near their own 0.85 female cutoff; 0.96 here), so "
            "this is a pre-existing borderline signal common to both methods, not a new artifact. "
            "Table 1's primary-record-validated call was not changed or challenged; this claim's "
            "sex layer exists only to normalize the TLR/TRIM candidate list and explicitly defers "
            "to Table 1 (`TLR_TRIM_CNV_REPORT.md`)."
        ),
        "action": "corrected this pass (RESULTS_EN.md/FULL_MANUSCRIPT_EN.md/METHODS_EN.md/PAPER_DRAFT.md; ported to ARS-UCD2.0, two-generations paragraph updated)",
    },
    {
        "claim_id": "dq_read_based_in_cis_counting",
        "section": "Results — coordinate-independent copy number (read-based, in cis)",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": None,
        "legacy_value": "not present in legacy baseline (docx does not report a read-based in-cis DQ count)",
        "new_evidence_paths": (
            "docs/codex_reanalysis/DQ_READ_COPYCOUNT_REPORT.md",
            "results/codex_reanalysis_v1/03b_dq_read_copycount/dq_readbased_copynumber.tsv",
            "results/codex_reanalysis_v1/03b_dq_read_copycount/dq_read_copycount_status.json",
        ),
        "new_value": (
            "ARS-UCD2.0 target BAMs: DQA=1/DQB=1 for nearly all long reads (max in-cis DQA>=2 "
            "in 1/13 animals at 0.5% of that animal's long reads; DQB>=2 in 0/13); DQ paralogs "
            "span 49-141 kb vs. 22-30 kb median read length, so reads rarely span the array"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "New relative to the docx baseline and now re-run on ARS-UCD2.0 target BAMs "
            "(`scripts/codex_reanalysis/dq_read_copycount.py`). The manuscript itself states "
            "this method cannot resolve DQ copy number and immediately pivots to "
            "haplotype-keeping assembly for the real result -- honestly self-limiting, not an "
            "overclaim. Qualitative agreement with the legacy "
            "`results/cnv_corrected/dq_readbased_copynumber.tsv` (rare multicopy reads only)."
        ),
        "action": "corrected this pass (ported to ARS-UCD2.0; METHODS/RESULTS/PAPER_DRAFT/two-generations paragraph updated)",
    },
    {
        "claim_id": "wgs_drb3_typing",
        "section": "Results / Discussion — WGS-based DRB3 typing vs. PCR-SBT",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/DISCUSSION_EN.md", "docs/METHODS_EN.md"),
        "legacy_excerpt": "The homozygous animal ID08 provided an internal control, but the heterozygous results did not support high-resolution DRB3 typing from this dataset alone",
        "legacy_value": (
            "PCR-SBT could not be reliably reproduced from WGS alone; homozygous animal "
            "served as internal control confirming the failure is heterozygote-specific"
        ),
        "new_evidence_paths": (
            "docs/codex_reanalysis/DRB3_EXON2_TYPING_REPORT.md",
            "results/codex_reanalysis_v1/03d_drb3_exon2_typing/drb3_exon2_typing.tsv",
            "docs/codex_reanalysis/DRB3_WGS_TYPING_REPORT.md",
            "results/codex_reanalysis_v1/03c_drb3_wgs_typing/drb3_wgs_readbased_concordance.tsv",
            "results/codex_reanalysis_v1/06_local_assembly/dq_dr/dq_dr_recovery.tsv",
            "results/codex_reanalysis_v1/06_local_assembly/wide/wide_linkage.tsv",
        ),
        "new_value": (
            "Whole-read scoring (03c) did not reproduce PCR-SBT: top1 among PCR-SBT alleles in "
            "8/13, full diplotype from top alleles in 4/13, *016:01/*027:03 top1 despite absence "
            "from PCR-SBT in 5/13; these are the two full-length (801/780 bp) panel entries, so "
            "ranking by matching bases favoured length. Exon-2 typing (03d: 221-nt window shared "
            "by all 384 panel alleles, joint diplotype fit): best-cost diplotype equals PCR-SBT "
            "in 12/13 at every window length 145-221 nt, at the resolution of exon-2 sequence "
            "groups; the exception is the homozygous ID08, fitted with a second allele one "
            "nucleotide away (cost gap 7). DQ/DR assembly best calls were all *016:01/*027:03 "
            "and share the length bias. Developed on the same 13 animals, no held-out set."
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "Superseded interpretation: the earlier statement that WGS could not reproduce "
            "PCR-SBT came from a scoring artifact in `drb3_wgs_typing.py` (kept unchanged as 03c). "
            "`scripts/codex_reanalysis/drb3_exon2_typing.py` (03d) restricts to exon 2 and fits "
            "a diplotype; PCR-SBT is used only for comparison. Wide-interval `homozygous_control` "
            "remains a probe-placement label, and the 2.4-kb *016:01 probe offset is close to "
            "the exon 2-exon 3 distance (2,441 bp), not tested directly."
        ),
        "action": (
            "corrected (exon-2 typing replaces the length-biased whole-read scoring; "
            "METHODS/RESULTS/DISCUSSION/ABSTRACT updated)"
        ),
    },
    {
        "claim_id": "figure1_structural_map_legend",
        "section": "Figure 1 legend — structural map of the BoLA region",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md"),
        "legacy_excerpt": (
            "The figure should distinguish reference depth, local assembly, and marker "
            "genotyping and should mark rs29026690 within ABHD16A."
        ),
        "legacy_value": (
            "docx directive only (schematic not yet drawn): distinguish reference depth / "
            "local assembly / marker genotyping per interval; mark rs29026690 within ABHD16A"
        ),
        "new_evidence_paths": (
            "data/meta/bola_cnv_genes.v2.bed",
            "data/meta/btn_gene_coordinates.bed",
            "results/cnv_corrected",
        ),
        "new_value": (
            "current legend names class II DR/DQ, class I, TRIM, and both BTN clusters with "
            "their coordinates and the two confirmation methods (NCBI annotation + IPD local "
            "alignment for class I/II/TRIM, NCBI Gene E-utilities for BTN); rs29026690 itself "
            "is not depicted in Fig. 1 (it is a class III GWAS marker, not part of the "
            "structural-map interval set) but is repositioned within ABHD16A in Fig. 6/Table 5"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "First of the figure_legends SCOPE_LIMITATIONS rows (see EVIDENCE_REVIEW addendum "
            "\"Figure-legend compliance pass\"). The docx directive's reference-depth/local-"
            "assembly/marker-genotyping distinction is satisfied structurally, one method per "
            "interval, but rs29026690's specific ABHD16A placement lives in Fig. 6, not Fig. 1 "
            "-- the docx wrote both instructions into one schematic that the current figure set "
            "splits across two; not a gap, but worth recording since a literal reading of the "
            "directive expects it in Fig. 1."
        ),
        "action": "keep as-is (rs29026690/ABHD16A placement is correctly shown in Fig. 6 instead)",
    },
    {
        "claim_id": "figure2_dq_legend",
        "section": "Figure 2 legend — DQ region reference-paralog absence and assembly resolution",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md"),
        "legacy_excerpt": "ID11 should be marked as under-assembled.",
        "legacy_value": (
            "docx directive: panel (A)-(C) depth/assembly/conceptual schematic, with ID11 "
            "explicitly flagged as under-assembled rather than left as an ambiguous data point"
        ),
        "new_evidence_paths": (
            "results/cnv_corrected/dq_readbased_copynumber.tsv",
            "results/cnv_corrected/dq_hapasm_copynumber.tsv",
            "docs/codex_reanalysis/DEPTH_CNV_REPORT.md",
        ),
        "new_value": (
            "Fig. 2 legend panel (A) states \"ID11 DQA2 ratio 0.00\" as a depth value, not an "
            "under-assembly annotation; the under-assembly flag itself is carried in Table 2 "
            "(RESULTS_EN.md: \"ID11 | ... | under-assembly (excluded)\"), which Fig. 2's own "
            "section header cites alongside the figure (\"Fig. 2A, Table 2\")"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "The directive is satisfied by the paired Table 2, not by the figure legend text "
            "itself -- the legend gives ID11's raw depth ratio (0.00) while the under-assembly "
            "judgment call sits in the table it is presented alongside. Functionally compliant "
            "(a reader cannot see the 0.00 without also seeing the table's exclusion label), but "
            "not literally the same location the docx specifies; recorded rather than silently "
            "treated as fully satisfied."
        ),
        "action": "keep as-is (under-assembly flag present in the paired Table 2)",
    },
    {
        "claim_id": "figure3_btn_legend",
        "section": "Figure 3 legend — butyrophilin copy-number variation",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md", "docs/METHODS_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": "Values represent normalized depth, not integer copy number.",
        "legacy_value": "docx directive: heat map of 5 loci, default + MAPQ20, explicit not-integer-CN disclaimer required",
        "new_evidence_paths": (
            "results/btn/btn_cohort_cn.tsv",
            "results/btn/btn_cn_summary.tsv",
            "data/meta/btn_gene_coordinates.bed",
        ),
        "new_value": (
            "the required disclaimer is not in the Fig. 3 legend itself (whose own title and "
            "panel (D) both say \"copy-number\"/\"copy-number heatmap\") but is stated globally: "
            "METHODS_EN.md \"that report also states explicitly that a depth ratio is not a "
            "validated integer copy number\", and DISCUSSION_EN.md's Limitations paragraph, "
            "specific to BTN, \"the depth ratio represents a lower bound rather than an exact "
            "integer copy number\""
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "Same pattern as figure2_dq_legend: the docx's required disclaimer exists and is "
            "locus-specific to BTN (stronger than a generic caveat), but lives in Methods/"
            "Limitations rather than inside the Fig. 3 legend text itself. \"Copy-number "
            "variation\"/\"copy-number heatmap\" as the figure's own title and panel label is "
            "conventional shorthand in this literature (matching Fig. 4's parallel title), not "
            "a claim of validated integer CN on its own -- the disclaimer sits one click away "
            "in Methods/Discussion rather than restated per-figure."
        ),
        "action": "keep as-is (disclaimer present in Methods and BTN-specific Limitations text)",
    },
    {
        "claim_id": "figure4_classI_legend",
        "section": "Figure 4 legend — class I copy-number variation",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md"),
        "legacy_excerpt": (
            "The figure should use the corrected 28.49–28.72-Mb coordinate framework and avoid "
            "high-resolution allele claims."
        ),
        "legacy_value": "docx directive: per-animal NC1/JSP.1/classI_BOLA ratios, corrected 28.49-28.72 Mb frame, no allele-level claims",
        "new_evidence_paths": (
            "results/cnv_corrected/btn_classI_cohort_cn.tsv",
            "results/codex_reanalysis_v1/01_reference/coordinate_registry.tsv",
        ),
        "new_value": (
            "legend text: \"Copy-number variation in the true class I region (13 animals)... "
            "at JSP.1 and classical BoLA\"; the 28.49-28.72 Mb frame is stated in Fig. 1's "
            "legend and RESULTS_EN.md's structural-map section, which class_i_gene_coordinates "
            "already verifies against coordinate_registry.tsv; no allele-level (e.g. BoLA-*NN "
            "typing) claim appears anywhere in the Fig. 4 legend or its Results paragraph"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": False,
        "note": (
            "Both directive conditions hold: the coordinate frame is the same one already "
            "verified by class_i_gene_coordinates (independent support), and the legend/prose "
            "stay at the depth-ratio level throughout, never naming a specific classical allele."
        ),
        "action": "keep as-is",
    },
    {
        "claim_id": "figure5_016_haplotype_legend",
        "section": "Figure 5 legend — *016:01 haplotype structure",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md", "docs/METHODS_EN.md"),
        "legacy_excerpt": (
            "The descriptive conclusion is the absence of a consistent carrier-specific "
            "pattern; no inferential P value should be shown."
        ),
        "legacy_value": "docx directive: recovery counts for 6 carriers/7 non-carriers, descriptive-only conclusion, no P value",
        "new_evidence_paths": ("results/cnv_corrected/dq_hapasm_copynumber.tsv",),
        "new_value": (
            "Fig. 5 legend states \"Diploid-normal DQA2/DRB3 in *016:01 carriers... n = 6/13\"; "
            "no P value appears anywhere in the legend, and METHODS_EN.md states this globally "
            "for the whole manuscript: \"This is a descriptive study, and no formal statistical "
            "hypothesis testing or p-values are used\""
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": False,
        "note": "Both directive conditions hold; the no-P-value rule is manuscript-wide, not just this figure, so Fig. 5 satisfies it by construction.",
        "action": "keep as-is",
    },
    {
        "claim_id": "figure6_gwas_snp_legend",
        "section": "Figure 6 legend — GWAS SNP repositioning and genotyping",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md"),
        "legacy_excerpt": (
            "Within-cohort binary correlation between *016:01 carriage and rs29026690-T (r² "
            "= 0.74), explicitly labeled as descriptive and ascertained."
        ),
        "legacy_value": "docx directive: panel (A) positions, (B) genotype matrix, (C) r²=0.74 correlation labeled descriptive/ascertained",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv",
            "results/cntn/rs110616206_genotypes.tsv",
            "docs/tables/table5_gwas_snp_genotypes.tsv",
        ),
        "new_value": (
            "current Fig. 6 legend keeps only the docx's panels (A)/(B) (position + genotype "
            "matrix); the r²=0.74 correlation panel (docx panel C) was moved to Fig. 7 panel "
            "(C) instead, where it is labeled descriptive/ascertained (Discussion, Results \"A "
            "note on interpretation\" paragraph, and rs29026690_016_01_correlation's own matrix "
            "row) -- present and correctly hedged, just relocated to a different figure"
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "Panel migration, not data loss: a reader checking Fig. 6C specifically for the "
            "r²=0.74 correlation will not find it there, only in Fig. 7C. Recorded so the "
            "migration is documented rather than silently invisible to a docx-vs-current "
            "reviewer -- the hedging language itself (descriptive/ascertained) is intact at "
            "its new location and already covered by rs29026690_016_01_correlation."
        ),
        "action": "keep as-is (r² panel relocated to Fig. 7C, not dropped; documented here for traceability)",
    },
    {
        "claim_id": "figure7_extended_haplotype_legend",
        "section": "Figure 7 legend — extended MHC haplotype and recombination boundary",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/FULL_MANUSCRIPT_EN.md", "docs/DISCUSSION_EN.md"),
        "legacy_excerpt": (
            "The figure must not label a recombination hotspot, population LD block, or "
            "demonstrated extended haplotype."
        ),
        "legacy_value": (
            "docx directive: schematic separating genotype correlation from assembly "
            "non-bridging only; explicitly forbids labeling a recombination hotspot, "
            "population LD block, or demonstrated extended haplotype"
        ),
        "new_evidence_paths": (
            "docs/figures/haploblock_figdata.json",
            "docs/codex_reanalysis/MHC_PHASE_BLOCKS_REPORT.md",
            "results/codex_reanalysis_v1/04b_mhc_phase_blocks/all_boundary_hotspots_100kb.tsv",
            "results/codex_reanalysis_v1/04b_mhc_phase_blocks/landmark_co_membership.tsv",
            "docs/chatgpt_figures_submission/Figure7_correlation_phase_unresolved.svg",
        ),
        "new_value": (
            "before this pass, the Fig. 7 title (\"Extended MHC haplotype and recombination "
            "boundary...\"), its paired Results section header, and Discussion both stated the "
            "interpretation as an unhedged label, and Discussion used the literal word "
            "\"hotspot\" twice for this study's own ~26.25 Mb boundary -- a direct textual match "
            "to what the docx forbids. Fixed this pass: Fig. 7's title now reads \"Candidate "
            "extended MHC haplotype...\", its Results header reads \"...consistent with a "
            "candidate extended MHC haplotype\", the adjoining Results sentence is hedged as an "
            "interpretation rather than an assertion, and both Discussion/Results \"hotspot\" "
            "instances naming this study's boundary are now \"recombination-prone boundary\" "
            "with an explicit terminology note (mirrored in PAPER_DRAFT.md). Item 6 regenerated "
            "haploblock_figdata.json and the chatgpt Fig. 7 SVG from UCD2.0 phase-block outputs "
            "(peak 26.25 Mb, 21 breaks; rs29026690 x rs17872126 same-block 7/13)."
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": True,
        "limitation": True,
        "note": (
            "This is distinct from drb3_rs29026690_recombination_boundary and "
            "extended_haplotype_mechanism_hypothesis, which already judged the underlying "
            "*evidence* (WhatsHap phase-block discontinuity independently corroborated by Task "
            "11's reference-free wide-interval assembly, all 13 animals) as real and sufficient "
            "to supersede the docx's older exclusion -- that judgment is not revisited here and "
            "this row does not relitigate it. What this row catches is narrower: the figure "
            "*title* and its immediately paired Results header stated the interpretation as a "
            "settled label with no hedge, and Discussion twice used the specific word the docx "
            "names as forbidden (\"hotspot\") for this study's own finding, rather than only in "
            "its generic description of prior human-MHC literature (left as-is, since that use "
            "describes established literature, not this study's claim). Confirmed with the user "
            "before editing given this touches the manuscript's headline finding; user selected "
            "hedging to comply strictly with the docx's explicit red line. Abstract/Conclusions/"
            "Keywords already used softer language (\"interpreted as evidence of\", \"providing "
            "evidence (descriptive) for\") and were left untouched -- only the unhedged title/"
            "header/hotspot instances needed the fix."
        ),
        "action": "figure title, Results header, adjoining sentence, and two Discussion \"hotspot\" instances hedged this pass in RESULTS_EN.md/FULL_MANUSCRIPT_EN.md/DISCUSSION_EN.md/PAPER_DRAFT.md",
    },
    {
        "claim_id": "figureS1_class3_scan_not_retained",
        "section": "Supplementary Figure S1 — exploratory class III candidate-variant scan",
        "manuscript_files": ("docs/RESULTS_EN.md", "docs/PAPER_DRAFT.md"),
        "legacy_excerpt": (
            "If retained, show the distribution of candidate sites across NCR3–TNF–LTA and "
            "state prominently that the mismatch-only calls were not orthogonally validated, "
            "were selected in the same 13 animals, and cannot establish population association."
        ),
        "legacy_value": "docx directive is conditional (\"if retained\"); no unconditional requirement to include Fig. S1",
        "new_evidence_paths": (
            "results/codex_reanalysis_v1/08_manuscript/figures/FigureS1_classIII_exploratory.svg",
        ),
        "new_value": (
            "Final decision (2026-09-03): the manuscript's AI peer review "
            "(docs/AI_REVIEWER_COMMENTS.md, Major 5) judged the TNF/Fig. S1 exploratory "
            "candidate-variant content out of place in the main Results of a hypothesis-generating "
            "descriptive study. RESULTS_EN.md/DISCUSSION_EN.md were revised accordingly to drop the "
            "\"(Fig. S1)\" citation and the Figure S1 legend, and to end the TNF -824A/G paragraph "
            "with \"Complementary exploratory pileups outside the validated calling pipeline are "
            "not reported as results of this study.\" Presented with the resulting conflict against "
            "an earlier user-requested retention, the user explicitly chose to support the AI "
            "reviewer's judgment and remove Fig. S1 from the main manuscript; PAPER_DRAFT.md was "
            "brought back into sync with this English chain in the same pass."
        ),
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": True,
        "note": (
            "The underlying evidence is unchanged: the class III exploratory scan remains the same "
            "ad hoc, non-canonical, unvalidated pileup result. Only its visibility changed -- it is "
            "no longer cited as a figure or discussed in the main text -- not its evidentiary "
            "weight, which was always sub-`accepted`-threshold and hypothesis-generating only."
        ),
        "action": "reverted this pass to not-retained per the AI reviewer's judgment and the user's explicit decision; Fig. S1 citation removed from RESULTS_EN.md/PAPER_DRAFT.md, no legend rendered as a cited figure",
    },
)


def contains_real_sample_id(text: str) -> bool:
    """True if ``text`` contains a real sample_id token (DELIVERY_/DELIVERY_/DELIVERY_)."""
    return bool(_REAL_SAMPLE_ID_RE.search(text))


def assert_no_identity_leak(rows: Sequence[dict]) -> None:
    """Raise if any matrix row would republish the anonymization key.

    The matrix is tracked; the manuscript side of every row uses anonymized IDs.
    A real sample_id anywhere in a row would let a reader pair it with the row's
    anonymized-ID claim and reconstruct `docs/ANONYMIZATION_KEY_INTERNAL.md`.
    """
    for row in rows:
        for field, value in row.items():
            text = value if isinstance(value, str) else " ".join(value) if isinstance(value, (list, tuple)) else ""
            if contains_real_sample_id(text):
                raise EvidenceMatrixError(
                    f"claim {row.get('claim_id')!r} field {field!r} contains a real sample_id "
                    "token; evidence must be cited by file/field, never by real sample_id"
                )


def validate_evidence_paths(paths: Sequence[str], repo_root: Path) -> None:
    """Raise if any cited evidence path does not exist on disk."""
    missing = [path for path in paths if not (repo_root / path).exists()]
    if missing:
        raise EvidenceMatrixError(f"evidence path(s) do not exist: {missing}")


def build_matrix(claims: Sequence[dict], repo_root: Path) -> list[dict]:
    """Validate each claim and attach its computed classification."""
    rows = []
    seen_ids: set[str] = set()
    for claim in claims:
        claim_id = claim["claim_id"]
        if claim_id in seen_ids:
            raise EvidenceMatrixError(f"duplicate claim_id: {claim_id!r}")
        seen_ids.add(claim_id)

        if claim["source_exists"]:
            validate_evidence_paths(claim["new_evidence_paths"], repo_root)

        classification = classify_claim(
            source_exists=claim["source_exists"],
            reproducible=claim["reproducible"],
            independent_support=claim["independent_support"],
            limitation=claim["limitation"],
        )
        row = dict(claim)
        row["classification"] = classification
        rows.append(row)

    assert_no_identity_leak(rows)
    return rows


def write_matrix_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "claim_id",
        "section",
        "classification",
        "legacy_excerpt",
        "legacy_value",
        "new_value",
        "manuscript_files",
        "new_evidence_paths",
        "note",
        "action",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            flat["legacy_excerpt"] = row["legacy_excerpt"] or ""
            flat["manuscript_files"] = "; ".join(row["manuscript_files"])
            flat["new_evidence_paths"] = "; ".join(row["new_evidence_paths"])
            writer.writerow(flat)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="07_evidence_audit directory for status/provenance")
    parser.add_argument("--matrix", type=Path, required=True, help="docs/codex_reanalysis/EVIDENCE_MATRIX.tsv destination")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    rows = build_matrix(CLAIMS, args.repo_root)

    args.matrix.parent.mkdir(parents=True, exist_ok=True)
    write_matrix_tsv(args.matrix, rows)

    counts = {
        classification: sum(1 for row in rows if row["classification"] == classification)
        for classification in (ACCEPTED, LIMITED, REJECTED)
    }
    args.root.mkdir(parents=True, exist_ok=True)
    # Never "passed" alone when SCOPE_LIMITATIONS is non-empty: acceptance must not
    # come from absent evidence (an unaudited manuscript section is not the same as
    # an audited section with no finding), matching pilot_gate.py's not_run/incomplete
    # convention rather than letting a partial sweep read as a complete one.
    status_payload = {
        "claim_count": len(rows),
        **{f"{k}_count": v for k, v in counts.items()},
        "scope_limitations": dict(SCOPE_LIMITATIONS),
        "status": "passed_with_scope_limitation" if SCOPE_LIMITATIONS else "passed",
    }
    (args.root / "evidence_matrix_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(args.root, {"input_hashes": {}})
    print(f"Evidence matrix written: {args.matrix} ({len(rows)} claims: {counts})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
