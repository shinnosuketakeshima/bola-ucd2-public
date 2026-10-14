from __future__ import annotations

from pathlib import Path

import pytest

from scripts.codex_reanalysis.check_manuscript_numbers import (
    ManuscriptNumberError,
    check_no_rejected_claims_in_manuscript,
    compare_numeric_claim,
    contains_rejected_claim,
    literal_found,
    phi_squared,
    render_report,
    run_checks,
)


def test_manuscript_table_count_matches_evidence_table():
    assert compare_numeric_claim("r_squared", manuscript=0.7347, evidence=0.7346938775510203)


def test_compare_numeric_claim_rejects_out_of_tolerance():
    assert not compare_numeric_claim("x", manuscript=0.74, evidence=0.60, tolerance=0.006)


def test_rejected_claim_cannot_remain_in_revised_text():
    assert contains_rejected_claim("...NCR3 at 27,713,795 in a receptor-ligand pair...", {"27,713,795"})
    assert not contains_rejected_claim("no such coordinate mentioned here", {"27,713,795"})


def test_phi_squared_matches_known_2x2_table():
    # 6 carrier+T-positive, 0 carrier+T-negative, 1 non-carrier+T-positive, 6 non-carrier+T-negative
    # -- the real rs29026690 x *016:01 table this project's manuscript reports as r^2 = 0.7347.
    assert phi_squared(6, 0, 1, 6) == pytest.approx(0.7346938775510203)


def test_phi_squared_raises_on_zero_margin():
    with pytest.raises(ManuscriptNumberError):
        phi_squared(0, 0, 5, 5)


def test_literal_found_true_for_substring():
    assert literal_found("cohort median ratio of 1.32 (range...)", "1.32")


def test_literal_found_false_when_absent():
    assert not literal_found("no numbers here", "1.32")


def test_check_no_rejected_claims_returns_empty_when_clean(tmp_path):
    clean = tmp_path / "clean.md"
    clean.write_text("Nothing about that coordinate here.", encoding="utf-8")
    assert check_no_rejected_claims_in_manuscript(tmp_path, ["clean.md"]) == []


def test_check_no_rejected_claims_flags_leaked_file(tmp_path):
    leaked = tmp_path / "leaked.md"
    leaked.write_text("The variant at 27,713,795 lies in NCR3.", encoding="utf-8")
    assert check_no_rejected_claims_in_manuscript(tmp_path, ["leaked.md"]) == ["leaked.md"]


def test_run_checks_raises_on_missing_manuscript_file(tmp_path):
    checks = [
        {
            "name": "missing_file_check",
            "manuscript_files": ("does/not/exist.md",),
            "manuscript_literal": "1.0",
            "compute": lambda repo_root: 1.0,
            "tolerance": 0,
        }
    ]
    with pytest.raises(ManuscriptNumberError, match="not found"):
        run_checks(checks, tmp_path)


def test_run_checks_matches_when_literal_found_and_value_within_tolerance(tmp_path):
    manuscript = tmp_path / "results.md"
    manuscript.write_text("The ratio was 1.32 in this cohort.", encoding="utf-8")
    checks = [
        {
            "name": "ratio_check",
            "manuscript_files": ("results.md",),
            "manuscript_literal": "1.32",
            "compute": lambda repo_root: 1.3239,
            "tolerance": 0.006,
        }
    ]
    results = run_checks(checks, tmp_path)
    assert results[0]["matches"] is True
    assert results[0]["found_in"] == ["results.md"]


def test_run_checks_fails_when_literal_absent_even_if_value_matches(tmp_path):
    manuscript = tmp_path / "results.md"
    manuscript.write_text("No relevant number is written here.", encoding="utf-8")
    checks = [
        {
            "name": "ratio_check",
            "manuscript_files": ("results.md",),
            "manuscript_literal": "1.32",
            "compute": lambda repo_root: 1.3239,
            "tolerance": 0.006,
        }
    ]
    results = run_checks(checks, tmp_path)
    assert results[0]["matches"] is False
    assert results[0]["found_in"] == []


def test_run_checks_fails_when_value_out_of_tolerance_even_if_literal_found(tmp_path):
    manuscript = tmp_path / "results.md"
    manuscript.write_text("The ratio was 1.32 in this cohort.", encoding="utf-8")
    checks = [
        {
            "name": "ratio_check",
            "manuscript_files": ("results.md",),
            "manuscript_literal": "1.32",
            "compute": lambda repo_root: 2.50,
            "tolerance": 0.006,
        }
    ]
    results = run_checks(checks, tmp_path)
    assert results[0]["matches"] is False


def test_render_report_includes_pass_fail_and_rejected_section():
    results = [
        {
            "name": "a",
            "manuscript_literal": "1.0",
            "found_in": ["results.md"],
            "evidence_value": 1.0,
            "matches": True,
        }
    ]
    rendered = render_report(results, rejected_claim_hits=[])
    assert "a" in rendered
    assert "PASS" in rendered
    assert "Rejected-claim leak check" in rendered


def test_render_report_flags_rejected_claim_hits():
    rendered = render_report([], rejected_claim_hits=["docs/RESULTS_EN.md"])
    assert "FAIL" in rendered
    assert "docs/RESULTS_EN.md" in rendered
