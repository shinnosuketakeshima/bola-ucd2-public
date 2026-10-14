from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from scripts.codex_reanalysis.compute_depth import (
    build_depth_command,
    parse_depth_output,
    require_complete_depth,
    sha256_paths,
    summarize_depth,
)
from scripts.codex_reanalysis.normalize_depth import (
    UNREPRESENTED_EXCEPTION,
    build_normalized_matrix,
    control_median_depth,
    interpretation_flag,
    locus_range_rows,
    normalized_ratio,
    render_report,
)


def test_zero_depth_positions_are_retained():
    assert summarize_depth([(1, 0), (2, 4)]) == {"median": 2.0, "zero_positions": 1}


def test_summarize_depth_handles_no_positions():
    assert summarize_depth([]) == {"median": 0.0, "zero_positions": 0}


def test_normalized_ratio_uses_control_median():
    assert normalized_ratio(locus_median=8, control_median=4) == 2.0


def test_normalized_ratio_rejects_nonpositive_control():
    with pytest.raises(ValueError):
        normalized_ratio(locus_median=8, control_median=0)


def test_control_median_depth_uses_median_of_controls():
    assert control_median_depth({"a": 10.0, "b": 30.0, "c": 20.0}) == 20.0


def test_control_median_depth_requires_at_least_one_control():
    with pytest.raises(ValueError):
        control_median_depth({})


def test_interpretation_flags_low_control_depth_first():
    result = interpretation_flag(
        ratio=1.0, cohort_median=1.0, cohort_mad=0.1, control_median=2.0, minimum_control_depth=10.0
    )
    assert result == "insufficient_control_signal"


def test_interpretation_flags_consistent_ratio():
    result = interpretation_flag(
        ratio=1.0, cohort_median=1.0, cohort_mad=0.1, control_median=30.0, minimum_control_depth=10.0
    )
    assert result == "consistent_with_cohort"


def test_interpretation_flags_elevated_outlier():
    result = interpretation_flag(
        ratio=5.0, cohort_median=1.0, cohort_mad=0.1, control_median=30.0, minimum_control_depth=10.0
    )
    assert result == "elevated_relative_to_cohort"


def test_interpretation_flags_reduced_outlier():
    result = interpretation_flag(
        ratio=0.0, cohort_median=1.0, cohort_mad=0.1, control_median=30.0, minimum_control_depth=10.0
    )
    assert result == "reduced_relative_to_cohort"


def test_depth_command_uses_one_based_region_and_mapq():
    cmd = build_depth_command("s.bam", "NC_1", 10, 20, 20, samtools="samtools")
    assert cmd == ["samtools", "depth", "-a", "-Q", "20", "-r", "NC_1:11-20", "s.bam"]


def test_parse_depth_output_reads_position_and_depth():
    text = "chr1\t5\t0\nchr1\t6\t3\n"
    assert parse_depth_output(text) == [(5, 0), (6, 3)]


def test_incomplete_depth_positions_are_rejected():
    with pytest.raises(ValueError, match="expected 3 positions"):
        require_complete_depth([(1, 0), (2, 4)], expected_length=3, locus="DQA5", mapq=20)


def test_complete_depth_positions_are_accepted():
    require_complete_depth([(1, 0), (2, 4), (3, 4)], expected_length=3, locus="DQA5", mapq=20)


def test_sha256_paths_returns_64_hex(tmp_path: Path):
    path = tmp_path / "example.bed"
    path.write_bytes(b"NC_1\t0\t10\tDQA5\n")
    hashes = sha256_paths({"bed": path})
    assert hashes == {"bed": hashlib.sha256(b"NC_1\t0\t10\tDQA5\n").hexdigest()}
    assert len(hashes["bed"]) == 64


def test_zero_median_nonzero_sample_is_not_elevated():
    result = interpretation_flag(
        ratio=0.4706, cohort_median=0.0, cohort_mad=0.0, control_median=25.5, minimum_control_depth=10.0
    )
    assert result == UNREPRESENTED_EXCEPTION
    assert result != "elevated_relative_to_cohort"


def test_zero_median_zero_sample_is_consistent():
    result = interpretation_flag(
        ratio=0.0, cohort_median=0.0, cohort_mad=0.0, control_median=25.5, minimum_control_depth=10.0
    )
    assert result == "consistent_with_cohort"


def test_mad_zero_positive_median_keeps_direction():
    elevated = interpretation_flag(
        ratio=2.0, cohort_median=1.0, cohort_mad=0.0, control_median=30.0, minimum_control_depth=10.0
    )
    reduced = interpretation_flag(
        ratio=0.5, cohort_median=1.0, cohort_mad=0.0, control_median=30.0, minimum_control_depth=10.0
    )
    assert elevated == "elevated_relative_to_cohort"
    assert reduced == "reduced_relative_to_cohort"


def _synthetic_cohort() -> dict[str, list[dict]]:
    def rows(dqa5_median: float, control_median: float = 30.0) -> list[dict]:
        return [
            {"locus": "DQA5", "mapq": 20, "median": dqa5_median, "breadth": 1.0 if dqa5_median else 0.0},
            {"locus": "BTNL2", "mapq": 20, "median": 30.0, "breadth": 1.0},
            {"locus": "control_ALB", "mapq": 20, "median": control_median, "breadth": 1.0},
        ]

    return {"a": rows(0.0), "b": rows(0.0), "c": rows(12.0)}


def test_build_normalized_matrix_marks_unrepresented_exception():
    rows = build_normalized_matrix(_synthetic_cohort())
    dqa5 = {row["sample_id"]: row for row in rows if row["locus"] == "DQA5"}
    btnl2 = {row["sample_id"]: row for row in rows if row["locus"] == "BTNL2"}
    assert dqa5["a"]["interpretation"] == "consistent_with_cohort"
    assert dqa5["c"]["normalized_ratio"] == 0.4
    assert dqa5["c"]["interpretation"] == UNREPRESENTED_EXCEPTION
    assert all(row["interpretation"] == "consistent_with_cohort" for row in btnl2.values())


def test_locus_ranges_and_report_separate_unrepresented_exceptions():
    rows = build_normalized_matrix(_synthetic_cohort())
    ranges = {row["locus"]: row for row in locus_range_rows(rows, mapq=20)}
    assert ranges["DQA5"]["minimum"] == 0.0
    assert ranges["DQA5"]["maximum"] == 0.4
    assert ranges["DQA5"]["n_zero"] == 2
    report = render_report(rows, ["a", "b", "c"])
    assert "present_at_otherwise_unrepresented_locus" in report
    assert "| DQA5 | c |" in report
    assert "elevated_relative_to_cohort" not in report
    assert "| DQA5 | 0.0 | 0.0 | 0.4 |" in report
    assert "modified z-score > 3.5" in report
    assert "not a validated integer copy number" in report.lower() or "not a validated integer" in report
