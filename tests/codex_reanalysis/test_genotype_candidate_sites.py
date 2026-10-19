from pathlib import Path

import pytest

from scripts.codex_reanalysis.genotype_candidate_sites import (
    CandidateSiteError,
    frequency_class,
    judge_genotype,
    load_candidate_sites,
    summarize_site,
)


def test_reference_call_needs_depth_and_agreement_but_no_gq():
    assert judge_genotype(20, 1, "0/0", None, 10, 20) == ("called", "0/0", "")


def test_missing_record_is_never_reference_when_depth_is_low():
    assert judge_genotype(7, 0, "0/0", None, 10, 20) == ("unresolved", "./.", "depth_below_minimum")
    assert judge_genotype(0, 0, "./.", None, 10, 20) == ("unresolved", "./.", "no_data")


def test_allele_fraction_and_caller_must_agree():
    # 4/24 ALT reads is het by the 15% rule but the caller says 0/0.
    assert judge_genotype(20, 4, "0/0", None, 10, 20)[2] == "caller_disagrees_with_allele_fraction"


def test_non_reference_call_needs_gq():
    assert judge_genotype(3, 20, "1/1", 3, 10, 20)[2] == "genotype_quality_below_minimum"
    assert judge_genotype(10, 14, "0/1", 87, 10, 20) == ("called", "0/1", "")


def test_frequency_classes_are_not_called_fixed():
    assert frequency_class(1.0) == "high"
    assert frequency_class(0.85) == "high"
    assert frequency_class(0.5) == "intermediate"
    assert frequency_class(0.149) == "low"


def test_summary_uses_called_alleles_as_denominator():
    site = {"gene": "ABHD16A", "contig": "c", "position": 1, "ref": "A", "alt": "G", "detail": ""}
    calls = (
        [{"status": "called", "genotype": "0/1", "depth": 20}] * 7
        + [{"status": "called", "genotype": "1/1", "depth": 20}] * 2
        + [{"status": "called", "genotype": "0/0", "depth": 20}] * 3
        + [{"status": "unresolved", "genotype": "./.", "depth": 23}]
    )
    summary = summarize_site(site, calls)
    assert summary["called_alleles"] == 24
    assert summary["alt_alleles"] == 11
    assert summary["alt_af"] == "0.458"
    assert summary["carriers"] == 9
    assert summary["animals_unresolved"] == 1
    assert summary["frequency_class"] == "intermediate"


def test_load_candidate_sites_keeps_only_target_missense(tmp_path: Path):
    path = tmp_path / "coding_variants.tsv"
    path.write_text(
        "gene\tpos\tref\talt\teffect\tdetail\n"
        "PRR3\t28420019\tC\tA\tMISSENSE\tQ>H\n"
        "PRR3\t28420020\tC\tA\tsynonymous\t\n"
        "BTNL2\t25860000\tC\tA\tMISSENSE\t\n",
        encoding="utf-8",
    )
    sites = load_candidate_sites(path)
    assert [(s["gene"], s["contig"], s["position"]) for s in sites] == [("PRR3", "NC_037350.1", 28420019)]


def test_load_candidate_sites_rejects_empty(tmp_path: Path):
    path = tmp_path / "coding_variants.tsv"
    path.write_text("gene\tpos\tref\talt\teffect\tdetail\n", encoding="utf-8")
    with pytest.raises(CandidateSiteError):
        load_candidate_sites(path)
