from __future__ import annotations

from scripts.codex_reanalysis.validate_cds_models import (
    choose_canonical_transcript,
    translate,
    validate_orf,
)


def test_translate_stops_at_incomplete_trailing_codon():
    assert translate("ATGAAATAG") == "MK*"
    assert translate("ATGAAATA") == "MK"  # trailing 2 bases dropped, not padded


def test_validate_orf_accepts_clean_single_orf():
    is_clean, reason = validate_orf("MK*", 9)
    assert is_clean
    assert "clean" in reason


def test_validate_orf_rejects_non_multiple_of_three():
    is_clean, _ = validate_orf("MK*", 10)
    assert not is_clean


def test_validate_orf_rejects_missing_start_methionine():
    is_clean, reason = validate_orf("KM*", 9)
    assert not is_clean
    assert "Met" in reason


def test_validate_orf_rejects_internal_stop():
    is_clean, reason = validate_orf("M*K*", 12)
    assert not is_clean
    assert "one stop" in reason or "premature" in reason


def test_validate_orf_rejects_missing_trailing_stop():
    is_clean, reason = validate_orf("MK", 6)
    assert not is_clean
    assert "one stop" in reason


def test_choose_canonical_transcript_prefers_curated_over_gnomon():
    by_transcript = {
        "XM_1.1": {"exons": [(1, 3000)], "strand": "+"},
        "NM_1.1": {"exons": [(1, 300)], "strand": "+"},
    }
    assert choose_canonical_transcript(by_transcript) == "NM_1.1"


def test_choose_canonical_transcript_ties_broken_by_longest_cds():
    by_transcript = {
        "XM_1.1": {"exons": [(1, 100)], "strand": "+"},
        "XM_2.1": {"exons": [(1, 300)], "strand": "+"},
    }
    assert choose_canonical_transcript(by_transcript) == "XM_2.1"
