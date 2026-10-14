from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.codex_reanalysis.audit_published_markers import (
    as_carrier,
    classify_marker,
    is_016_carrier,
    legacy_table5_concordance,
    mapping_bam_hashes,
    normalize_diplotype,
    render_report,
    require_complete_marker_matrix,
    row_from_genotype,
)


def test_marker_is_accepted_only_when_callers_agree():
    assert classify_marker(primary="C/T", pileup="C/T", depth=30) == "accepted"


def test_marker_is_limited_when_depth_is_insufficient():
    assert classify_marker(primary="C/T", pileup="C/T", depth=4) == "limited"


def test_marker_is_conflicting_when_pileups_disagree():
    assert classify_marker(primary="C/T", pileup="T/T", depth=30) == "conflicting"


def test_marker_is_no_data_when_either_pileup_has_no_data():
    assert classify_marker(primary="no_data", pileup="C/T", depth=30) == "no_data"
    assert classify_marker(primary="C/T", pileup="no_data", depth=30) == "no_data"


def test_marker_uses_configured_minimum_depth():
    assert classify_marker(primary="C/T", pileup="C/T", depth=9, minimum_depth=10) == "limited"
    assert classify_marker(primary="C/T", pileup="C/T", depth=10, minimum_depth=10) == "accepted"


def test_is_016_carrier_detects_either_haplotype():
    assert is_016_carrier(("*015:01", "*016:01"))
    assert is_016_carrier(("*016:01", "*016:01"))


def test_is_016_carrier_false_for_non_carrier():
    assert not is_016_carrier(("*010:01", "*010:01"))


def test_marker_is_limited_when_only_one_pileup_clears_depth():
    result = classify_marker(
        primary="C/T", pileup="C/T", depth=30, minimum_depth=10, pileup_depth=4
    )
    assert result == "limited"


def test_row_from_genotype_keeps_allele_depths_and_exploratory_caller():
    marker = {
        "marker_id": "rs29026690",
        "ucd2_0_contig": "NC_037350.1",
        "ucd2_0_position": "27603651",
        "resolution_status": "verified_identical_chromosome",
    }
    sample = {
        "sample_id": "ID02",
        "pcr_sbt_drb3_hap1": "*015:01",
        "pcr_sbt_drb3_hap2": "*016:01",
    }
    row = row_from_genotype(
        marker,
        sample,
        {
            "primary_genotype": "C/T",
            "primary_depth": 28,
            "primary_ref_count": 14,
            "primary_alt_count": 14,
            "pileup_genotype": "C/T",
            "pileup_depth": 27,
            "pileup_ref_count": 13,
            "pileup_alt_count": 14,
        },
        minimum_depth=10,
    )
    assert row["primary_ref_count"] == 14
    assert row["pileup_alt_count"] == 14
    assert row["caller_status"] == "exploratory"
    assert row["evidence_status"] == "accepted"
    assert row["carrier_016_01"] is True
    assert row["position"] == 27603651


def test_unresolved_marker_is_coordinate_unresolved_not_genotyped():
    marker = {
        "marker_id": "rs_test",
        "ucd2_0_contig": "NC_037350.1",
        "ucd2_0_position": "",
        "resolution_status": "decision_required",
    }
    sample = {"sample_id": "a", "pcr_sbt_drb3_hap1": "*010:01", "pcr_sbt_drb3_hap2": "*010:01"}
    row = row_from_genotype(marker, sample, None, minimum_depth=10)
    assert row["evidence_status"] == "coordinate_unresolved"
    assert row["primary_genotype"] == "no_data"
    assert row["primary_ref_count"] == 0


def test_require_complete_marker_matrix_rejects_a_missing_sample():
    rows = [
        {"marker_id": "rs1", "sample_id": "a"},
        {"marker_id": "rs1", "sample_id": "b"},
    ]
    with pytest.raises(ValueError, match="missing"):
        require_complete_marker_matrix(rows, ["a", "b", "c"], ["rs1"])


def test_require_complete_marker_matrix_rejects_duplicates():
    rows = [
        {"marker_id": "rs1", "sample_id": "a"},
        {"marker_id": "rs1", "sample_id": "a"},
    ]
    with pytest.raises(ValueError, match="duplicate"):
        require_complete_marker_matrix(rows, ["a"], ["rs1"])


def test_mapping_bam_hashes_read_status_digests(tmp_path: Path):
    digest = hashlib.sha256(b"bam").hexdigest()
    payload = {"artifact_hashes": {"target_bam": digest}, "status": "passed"}
    (tmp_path / "s1.mapping_status.json").write_text(json.dumps(payload), encoding="utf-8")
    assert mapping_bam_hashes(tmp_path, ["s1"]) == {"s1": digest}


def test_mapping_bam_hashes_reject_placeholder(tmp_path: Path):
    payload = {"artifact_hashes": {"target_bam": "s1.targets.bam"}, "status": "passed"}
    (tmp_path / "s1.mapping_status.json").write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256"):
        mapping_bam_hashes(tmp_path, ["s1"])


def test_report_states_exploratory_and_keeps_accepted_from_meaning_caller():
    markers = [
        {
            "marker_id": "rs29026690",
            "ucd2_0_contig": "NC_037350.1",
            "ucd2_0_position": "27603651",
            "resolution_status": "verified_identical_chromosome",
            "citation": "test",
            "notes": "",
        }
    ]
    rows = [
        {
            "marker_id": "rs29026690",
            "sample_id": "ID02",
            "primary_genotype": "C/T",
            "primary_depth": 28,
            "primary_ref_count": 14,
            "primary_alt_count": 14,
            "pileup_genotype": "C/T",
            "pileup_depth": 27,
            "evidence_status": "accepted",
            "carrier_016_01": True,
        }
    ]
    report = render_report(rows, markers)
    assert "Caller status: exploratory" in report
    assert "No Clair3 model verified" in report
    assert "not that a verified caller confirmed it" in report
    assert "14/14" in report


def test_normalize_diplotype_sorts_alleles():
    assert normalize_diplotype("T/C") == "C/T"
    assert normalize_diplotype("C/C") == "C/C"


def test_false_carrier_string_from_tsv_is_not_yes():
    assert as_carrier(False) is False
    assert as_carrier("False") is False
    assert as_carrier("True") is True


def test_report_states_reference_allele_so_hom_alt_ad_is_not_misread():
    markers = [
        {
            "marker_id": "rs29026690",
            "ucd2_0_contig": "NC_037350.1",
            "ucd2_0_position": "27603651",
            "ucd2_0_reference_allele": "T",
            "resolution_status": "verified_identical_chromosome",
            "citation": "test",
            "notes": "",
        }
    ]
    rows = [
        {
            "marker_id": "rs29026690",
            "sample_id": "ID01",
            "primary_genotype": "C/C",
            "primary_depth": 6,
            "primary_ref_count": 0,
            "primary_alt_count": 6,
            "pileup_genotype": "C/C",
            "pileup_depth": 6,
            "evidence_status": "limited",
            "carrier_016_01": "False",
        }
    ]
    report = render_report(rows, markers)
    assert "reference allele t" in report.lower()
    assert "ref/alt counts against the linear-reference base" in report
    assert "| ID01 | C/C (0/6) |" in report
    assert "| ID01 |" in report and "| no |" in report.split("ID01")[1].split("\n")[0]


def test_legacy_table5_concordance_matches_sorted_alleles():
    genotype_rows = [
        {
            "marker_id": "rs29026690",
            "sample_id": "ID02",
            "primary_genotype": "C/T",
            "evidence_status": "limited",
        },
        {
            "marker_id": "rs110616206",
            "sample_id": "ID02",
            "primary_genotype": "C/C",
            "evidence_status": "accepted",
        },
        {
            "marker_id": "rs17872126",
            "sample_id": "ID02",
            "primary_genotype": "T/T",
            "evidence_status": "accepted",
        },
    ]
    legacy_rows = [
        {
            "sample_id": "ID02",
            "rs110616206_CNTN3": "C/C",
            "rs29026690_ABHD16A": "T/C",
            "rs17872126_PRR3": "T/T",
        }
    ]
    summary = legacy_table5_concordance(genotype_rows, legacy_rows)
    assert summary["n_compared"] == 3
    assert summary["n_mismatch"] == 0
    assert summary["n_limited_matching"] == 1
    report = render_report(genotype_rows, [
        {
            "marker_id": marker_id,
            "ucd2_0_contig": "NC_037350.1",
            "ucd2_0_position": "1",
            "ucd2_0_reference_allele": "C",
            "resolution_status": "verified_identical_chromosome",
            "citation": "",
            "notes": "",
        }
        for marker_id in ("rs110616206", "rs29026690", "rs17872126")
    ], legacy_rows=legacy_rows)
    assert "39/39" not in report
    assert "3/3" in report
    assert "not a claim re-established from accepted UCD2.0 pileups" in report
    assert "| rs29026690 | 1 | 1 | 0 | 1 | 0 |" in report
