from __future__ import annotations

import pytest

from scripts.codex_reanalysis.summarize_mapping_qc import (
    CohortError,
    flag_outliers,
    validate_cohort_ids,
)


def test_cohort_table_has_each_expected_sample_once():
    assert validate_cohort_ids(["ID01", "ID02"], ["ID01", "ID02"])


def test_cohort_table_rejects_duplicate_sample():
    with pytest.raises(CohortError):
        validate_cohort_ids(["ID01", "ID01"], ["ID01"])


def test_cohort_table_rejects_duplicate_observed_sample():
    with pytest.raises(CohortError):
        validate_cohort_ids(["ID01", "ID02"], ["ID01", "ID01"])


def test_cohort_table_rejects_missing_sample():
    with pytest.raises(CohortError, match="missing"):
        validate_cohort_ids(["ID01", "ID02"], ["ID01"])


def test_cohort_table_rejects_unexpected_sample():
    with pytest.raises(CohortError, match="unexpected"):
        validate_cohort_ids(["ID01"], ["ID01", "ID02"])


def test_flag_outliers_uses_median_absolute_deviation():
    values = {"a": 10.0, "b": 10.0, "c": 10.0, "d": 1.0}
    assert flag_outliers(values) == {"a": False, "b": False, "c": False, "d": True}


def test_flag_outliers_reports_no_flags_when_values_are_identical():
    values = {"a": 5.0, "b": 5.0, "c": 5.0}
    assert flag_outliers(values) == {"a": False, "b": False, "c": False}
