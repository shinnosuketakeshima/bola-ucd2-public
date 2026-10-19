import pytest

from scripts.codex_reanalysis.depth_flag_sensitivity import compare


def rows_for(locus, ratios, flags):
    return [
        {"locus": locus, "mapq": "20", "sample_id": s, "normalized_ratio": str(r), "control_median_depth": "30",
         "interpretation": f}
        for (s, r), f in zip(ratios.items(), flags)
    ]


def test_flag_that_depends_on_the_excluded_animal_is_reported():
    ratios = {"a": 1.0, "b": 1.0, "c": 1.02, "d": 1.3, "x": 2.0}
    # with x present the MAD is 0.02 and d is flagged; without x, d still exceeds the median strongly
    rows = rows_for("L", ratios, ["consistent_with_cohort"] * 3 + ["elevated_relative_to_cohort", "elevated_relative_to_cohort"])
    changes = compare(rows, "x", "20", 10.0)
    assert all(c["sample_id"] != "x" for c in changes)


def test_unchanged_flags_are_not_listed():
    ratios = {"a": 1.0, "b": 1.01, "c": 0.99, "d": 1.0, "e": 3.0, "x": 1.0}
    rows = rows_for("L", ratios, ["consistent_with_cohort"] * 4 + ["elevated_relative_to_cohort", "consistent_with_cohort"])
    assert compare(rows, "x", "20", 10.0) == []


def test_flag_lost_without_excluded_animal():
    ratios = {"a": 1.0, "b": 1.0, "c": 1.0, "d": 1.0, "e": 1.3, "x": 1.1}
    # all animals: median 1.0, MAD 0 -> e flagged; without x the MAD is still 0, e still flagged
    rows = rows_for("L", ratios, ["consistent_with_cohort"] * 4 + ["elevated_relative_to_cohort", "elevated_relative_to_cohort"])
    assert compare(rows, "x", "20", 10.0) == []
    ratios = {"a": 1.0, "b": 1.0, "c": 1.1, "d": 1.1, "e": 1.3, "x": 1.0}
    rows = rows_for("L", ratios, ["consistent_with_cohort"] * 4 + ["elevated_relative_to_cohort", "consistent_with_cohort"])
    (change,) = compare(rows, "x", "20", 10.0)
    assert change["sample_id"] == "e" and change["flag_without_x"] == "consistent_with_cohort"


def test_unknown_exclusion_is_rejected():
    with pytest.raises(ValueError):
        compare(rows_for("L", {"a": 1.0}, ["consistent_with_cohort"]), "x", "20", 10.0)
