import gzip
import random
import subprocess
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.codex_reanalysis.audit_inputs import FastqValidationError, validate_fastq


ASCII_WHITESPACE = tuple(code for code in range(128) if chr(code).isspace())


def naive_summary(records):
    """Independent per-character reference for the streaming summariser.

    This mirrors the original character-by-character implementation so that a
    faster implementation cannot silently change what the audit reports.
    """
    read_count = base_count = quality_total = 0
    min_length = None
    max_length = 0
    for _, sequence, quality in records:
        length = len(sequence)
        read_count += 1
        base_count += length
        min_length = length if min_length is None else min(min_length, length)
        max_length = max(max_length, length)
        quality_total += sum(ord(score) - 33 for score in quality)
    return {
        "read_count": read_count,
        "base_count": base_count,
        "min_read_length": min_length or 0,
        "max_read_length": max_length,
        "mean_read_length": (base_count / read_count) if read_count else 0.0,
        "mean_base_quality": (quality_total / base_count) if base_count else 0.0,
    }


def write_fastq(path, records):
    with path.open("w", encoding="ascii", newline="") as handle:
        for read_id, sequence, quality in records:
            handle.write(f"@{read_id}\n{sequence}\n+\n{quality}\n")


def test_audit_script_runs_as_a_direct_command():
    """The protocol and the Snakemake rule both invoke this file as a script.

    Without a path bootstrap the absolute package imports fail before argument
    parsing, so the audit cannot start at all.
    """
    result = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "audit_inputs.py"), "--help"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--deliver-root" in result.stdout


def test_validate_fastq_rejects_quality_length_mismatch(tmp_path):
    path = tmp_path / "bad.fastq"
    path.write_text("@r1\nACGT\n+\nIII\n", encoding="ascii")

    with pytest.raises(FastqValidationError, match="quality length"):
        validate_fastq(path)


def test_validate_fastq_streams_gzip_and_summarizes_quality(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
        handle.write("@read-1 extra\nACGT\n+\nIIII\n@read-2\nNN\n+\n!!\n")

    summary = validate_fastq(path)

    assert summary.read_count == 2
    assert summary.base_count == 6
    assert summary.min_read_length == 2
    assert summary.max_read_length == 4
    assert summary.mean_read_length == pytest.approx(3.0)
    assert summary.mean_base_quality == pytest.approx(26.6666666667)
    assert summary.duplicate_read_ids == 0


def test_validate_fastq_rejects_truncated_record(tmp_path):
    path = tmp_path / "truncated.fastq"
    path.write_text("@r1\nACGT\n+\n", encoding="ascii")

    with pytest.raises(FastqValidationError, match="truncated"):
        validate_fastq(path)


LINE_TERMINATORS = (10, 13)  # universal newlines end a line on LF and bare CR alike


@pytest.mark.parametrize(
    "code", [code for code in ASCII_WHITESPACE if code not in LINE_TERMINATORS]
)
def test_validate_fastq_rejects_any_ascii_whitespace_inside_a_sequence(tmp_path, code):
    path = tmp_path / "whitespace.fastq"
    path.write_text(f"@r1\nAC{chr(code)}GT\n+\nIIIII\n", encoding="ascii", newline="")

    with pytest.raises(FastqValidationError, match="invalid sequence characters"):
        validate_fastq(path)


def test_validate_fastq_rejects_a_bare_carriage_return_as_a_broken_record(tmp_path):
    """A bare CR ends a line for the reader, so it breaks record structure.

    It is still rejected, but as a malformed record rather than as an invalid
    base, and identically for plain and gzip inputs.
    """
    path = tmp_path / "carriage.fastq"
    path.write_text("@r1\nAC\rGT\n+\nIIIII\n", encoding="ascii", newline="")

    with pytest.raises(FastqValidationError):
        validate_fastq(path)

    gz_path = tmp_path / "carriage.fastq.gz"
    with gzip.open(gz_path, "wt", encoding="ascii", newline="") as handle:
        handle.write("@r1\nAC\rGT\n+\nIIIII\n")

    with pytest.raises(FastqValidationError):
        validate_fastq(gz_path)


@pytest.mark.parametrize("code", [0, 1, 7, 27, 127])
def test_validate_fastq_accepts_non_whitespace_control_characters_in_a_sequence(tmp_path, code):
    """Documents the current contract: only whitespace and non-ASCII are rejected.

    Control characters that are not whitespace are counted as bases today. This
    test exists so that a faster character scan cannot quietly tighten or
    loosen the rule without a deliberate decision.
    """
    path = tmp_path / "control.fastq"
    path.write_text(f"@r1\nAC{chr(code)}GT\n+\nIIIII\n", encoding="ascii", newline="")

    summary = validate_fastq(path)

    assert summary.read_count == 1
    assert summary.base_count == 5


def test_validate_fastq_reports_a_non_ascii_byte_as_a_malformed_stream(tmp_path):
    path = tmp_path / "latin1.fastq"
    path.write_bytes(b"@r1\nAC\xc3GT\n+\nIIIII\n")

    with pytest.raises(FastqValidationError, match="malformed or truncated"):
        validate_fastq(path)


@pytest.mark.parametrize("code", [32, 127, 0])
def test_validate_fastq_rejects_quality_outside_the_printable_range(tmp_path, code):
    path = tmp_path / "quality.fastq"
    path.write_text(f"@r1\nACGT\n+\nII{chr(code)}I\n", encoding="ascii", newline="")

    with pytest.raises(FastqValidationError, match="invalid quality characters"):
        validate_fastq(path)


def test_validate_fastq_accepts_the_whole_printable_quality_range(tmp_path):
    quality = "".join(chr(code) for code in range(33, 127))
    path = tmp_path / "range.fastq"
    path.write_text(f"@r1\n{'A' * len(quality)}\n+\n{quality}\n", encoding="ascii", newline="")

    summary = validate_fastq(path)

    expected = sum(code - 33 for code in range(33, 127)) / len(quality)
    assert summary.mean_base_quality == pytest.approx(expected)


def test_validate_fastq_summary_matches_a_naive_per_character_implementation(tmp_path):
    rng = random.Random(20260801)
    records = [
        (
            f"read-{index}",
            "".join(rng.choice("ACGTN") for _ in range(rng.randint(1, 400))),
            None,
        )
        for index in range(200)
    ]
    records = [
        (read_id, sequence, "".join(chr(rng.randint(33, 126)) for _ in sequence))
        for read_id, sequence, _ in records
    ]
    path = tmp_path / "random.fastq"
    write_fastq(path, records)

    summary = validate_fastq(path)
    expected = naive_summary(records)

    assert summary.read_count == expected["read_count"]
    assert summary.base_count == expected["base_count"]
    assert summary.min_read_length == expected["min_read_length"]
    assert summary.max_read_length == expected["max_read_length"]
    assert summary.mean_read_length == pytest.approx(expected["mean_read_length"])
    assert summary.mean_base_quality == pytest.approx(expected["mean_base_quality"])
