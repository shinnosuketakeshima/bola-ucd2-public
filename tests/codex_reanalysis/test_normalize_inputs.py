import gzip
import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.codex_reanalysis.normalize_inputs import (
    NormalizationError,
    filter_fastq_by_length,
    plan_normalization,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def qc_rows():
    return [
        {"sample_id": "keep_a", "min_read_length": "500", "read_count": "100"},
        {"sample_id": "keep_b", "min_read_length": "500", "read_count": "100"},
        {"sample_id": "short_a", "min_read_length": "12", "read_count": "180"},
    ]


def write_fastq(path, records):
    with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_plan_selects_only_samples_below_the_cohort_floor():
    plan = plan_normalization(qc_rows(), minimum_read_length=500)

    required = [row["sample_id"] for row in plan if row["action"] == "length_filter"]
    assert required == ["short_a"]


def test_plan_records_unaffected_samples_explicitly():
    """An unfiltered sample must be stated, never left out of the record."""
    plan = plan_normalization(qc_rows(), minimum_read_length=500)

    assert {row["sample_id"] for row in plan} == {"keep_a", "keep_b", "short_a"}
    unaffected = [row for row in plan if row["sample_id"] == "keep_a"][0]
    assert unaffected["action"] == "not_required"
    assert unaffected["observed_min_read_length"] == 500


def test_plan_rejects_a_sample_whose_minimum_length_is_unreadable():
    rows = [{"sample_id": "bad", "min_read_length": "", "read_count": "1"}]

    with pytest.raises(NormalizationError, match="minimum read length"):
        plan_normalization(rows, minimum_read_length=500)


def test_filter_drops_short_reads_and_keeps_the_rest_byte_for_byte(tmp_path):
    source = tmp_path / "in.fq.gz"
    destination = tmp_path / "out.fq.gz"
    write_fastq(source, [("keep1", "A" * 600), ("drop1", "C" * 499), ("keep2", "G" * 500)])

    summary = filter_fastq_by_length(source, destination, minimum_read_length=500)

    with gzip.open(destination, "rt", encoding="ascii", newline="") as handle:
        kept = handle.read()
    assert "keep1" in kept and "keep2" in kept
    assert "drop1" not in kept
    assert kept.count("@") == 2
    assert summary["reads_in"] == 3
    assert summary["reads_out"] == 2
    assert summary["reads_removed"] == 1
    assert summary["bases_in"] == 600 + 499 + 500
    assert summary["bases_out"] == 600 + 500
    assert summary["bases_removed"] == 499


def test_filter_keeps_a_read_exactly_at_the_threshold(tmp_path):
    source = tmp_path / "in.fq.gz"
    destination = tmp_path / "out.fq.gz"
    write_fastq(source, [("edge", "A" * 500)])

    summary = filter_fastq_by_length(source, destination, minimum_read_length=500)

    assert summary["reads_out"] == 1


def test_filter_rejects_a_truncated_record(tmp_path):
    source = tmp_path / "in.fq.gz"
    destination = tmp_path / "out.fq.gz"
    with gzip.open(source, "wt", encoding="ascii", newline="") as handle:
        handle.write("@r1\nACGT\n+\n")

    with pytest.raises(NormalizationError, match="truncated"):
        filter_fastq_by_length(source, destination, minimum_read_length=1)


def test_filter_rejects_a_quality_length_mismatch(tmp_path):
    source = tmp_path / "in.fq.gz"
    destination = tmp_path / "out.fq.gz"
    with gzip.open(source, "wt", encoding="ascii", newline="") as handle:
        handle.write("@r1\nACGT\n+\nIII\n")

    with pytest.raises(NormalizationError, match="quality length"):
        filter_fastq_by_length(source, destination, minimum_read_length=1)


def test_filter_does_not_leave_a_partial_output_when_the_source_is_broken(tmp_path):
    source = tmp_path / "in.fq.gz"
    destination = tmp_path / "out.fq.gz"
    with gzip.open(source, "wt", encoding="ascii", newline="") as handle:
        handle.write("@r1\n" + "A" * 600 + "\n+\n" + "I" * 600 + "\n@r2\nACGT\n+\n")

    with pytest.raises(NormalizationError):
        filter_fastq_by_length(source, destination, minimum_read_length=1)

    assert not destination.exists()


def test_normalize_script_runs_as_a_direct_command():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "normalize_inputs.py"),
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--audit-dir" in result.stdout
    assert "--minimum-read-length" in result.stdout
