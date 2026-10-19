from __future__ import annotations

import pytest

from scripts.codex_reanalysis.summarize_tlr_trim import (
    FEMALE,
    MALE,
    UNCERTAIN,
    TlrTrimSummaryError,
    autosomal_outliers,
    classify_sex,
    compute_sex_inference_ratios,
    find_sex_inference_gap,
)


def _row(locus, mapq, sample_id, ratio, interpretation="consistent_with_cohort"):
    return {
        "locus": locus,
        "mapq": mapq,
        "sample_id": sample_id,
        "normalized_ratio": str(ratio),
        "cohort_median_ratio": "1.0",
        "interpretation": interpretation,
    }


def test_compute_sex_inference_ratios_averages_tlr7_tlr8():
    rows = [
        _row("TLR7", "20", "s1", 0.4),
        _row("TLR8", "20", "s1", 0.6),
        _row("TLR7", "10", "s1", 99.0),  # wrong mapq, ignored
        _row("TLR9", "20", "s1", 5.0),  # not X-linked, ignored
    ]
    ratios = compute_sex_inference_ratios(rows)
    assert ratios == {"s1": 0.5}


def test_compute_sex_inference_ratios_raises_when_locus_missing():
    rows = [_row("TLR7", "20", "s1", 0.4)]
    with pytest.raises(TlrTrimSummaryError, match="missing TLR7 or TLR8"):
        compute_sex_inference_ratios(rows)


def test_find_sex_inference_gap_picks_widest_gap_midpoint():
    threshold, width = find_sex_inference_gap([0.3, 0.4, 0.5, 0.9, 1.0, 1.1])
    assert threshold == pytest.approx(0.7)
    assert width == pytest.approx(0.4)


def test_find_sex_inference_gap_requires_two_values():
    with pytest.raises(TlrTrimSummaryError):
        find_sex_inference_gap([0.5])


def test_classify_sex_below_threshold_minus_band_is_male():
    assert classify_sex(0.5, threshold=0.7, uncertainty_band=0.05) == MALE


def test_classify_sex_above_threshold_plus_band_is_female():
    assert classify_sex(0.9, threshold=0.7, uncertainty_band=0.05) == FEMALE


def test_classify_sex_within_band_is_uncertain():
    assert classify_sex(0.71, threshold=0.7, uncertainty_band=0.05) == UNCERTAIN


def test_autosomal_outliers_excludes_x_linked_and_wrong_mapq_and_consistent_rows():
    rows = [
        _row("TRIM10", "20", "s1", 1.5, interpretation="elevated_relative_to_cohort"),
        _row("TLR7", "20", "s1", 0.4, interpretation="elevated_relative_to_cohort"),  # X-linked, excluded
        _row("TRIM10", "10", "s1", 1.5, interpretation="elevated_relative_to_cohort"),  # wrong mapq
        _row("TRIM11", "20", "s1", 1.0, interpretation="consistent_with_cohort"),  # not an outlier
        _row("BOLA", "20", "s1", 5.0, interpretation="elevated_relative_to_cohort"),  # not TLR/TRIM
    ]
    outliers = autosomal_outliers(rows)
    assert len(outliers) == 1
    assert outliers[0]["locus"] == "TRIM10"
