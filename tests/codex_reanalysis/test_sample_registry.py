import pytest

from scripts.codex_reanalysis.sample_registry import (
    SampleRegistryError,
    build_rows_from_metadata,
    build_sample_registry,
)


def write_metadata(tmp_path, deliver_root, fastq_name):
    """Write the three tracked metadata tables the audit joins."""
    samplesheet = tmp_path / "samplesheet.tsv"
    samplesheet.write_text("sample_id\tbreed\nID01\tJapanese_Black\n", encoding="utf-8")
    mapping = tmp_path / "deliverable_mapping.tsv"
    mapping.write_text(
        "deliverable_id\tsample_id\tfastq_path\n"
        f"DELIVERY_01\tID01\t{fastq_name}\n",
        encoding="utf-8",
    )
    drb3 = tmp_path / "drb3_alleles.tsv"
    drb3.write_text("sample_id\tdrb3_hap1\tdrb3_hap2\nID01\t*015:01\t*016:01\n", encoding="utf-8")
    delivered = deliver_root / "Filtered"
    delivered.mkdir(parents=True, exist_ok=True)
    (delivered / "Filtered_DELIVERY_01_pass.fq.gz").write_bytes(b"gz")
    return samplesheet, mapping, drb3


def test_metadata_join_reads_a_windows_delivery_path_on_any_host(tmp_path):
    """The tracked metadata records Windows drive paths for the delivered FASTQ.

    Mapping runs under WSL, so the same tables must resolve identically on a
    POSIX host; a POSIX parser treats the whole backslash string as one name
    and would silently drop every sample.
    """
    deliver_root = tmp_path / "delivery"
    samplesheet, mapping, drb3 = write_metadata(
        tmp_path, deliver_root, r"H:\Deliver_ADD\Filtered\Filtered_DELIVERY_01_pass.fq.gz"
    )

    rows, expected_ids = build_rows_from_metadata(samplesheet, mapping, drb3, deliver_root)

    assert expected_ids == ["ID01"]
    assert len(rows) == 1
    assert rows[0]["fastq"] == "Filtered/Filtered_DELIVERY_01_pass.fq.gz"


def test_metadata_join_also_accepts_a_posix_delivery_path(tmp_path):
    deliver_root = tmp_path / "delivery"
    samplesheet, mapping, drb3 = write_metadata(
        tmp_path, deliver_root, "/mnt/h/Deliver_ADD/Filtered/Filtered_DELIVERY_01_pass.fq.gz"
    )

    rows, _ = build_rows_from_metadata(samplesheet, mapping, drb3, deliver_root)

    assert rows[0]["fastq"] == "Filtered/Filtered_DELIVERY_01_pass.fq.gz"


def test_sample_registry_requires_exactly_one_input_per_expected_sample():
    rows = [{"sample_id": "ID01", "fastq": "a.fastq.gz"}]

    with pytest.raises(SampleRegistryError, match="missing expected samples"):
        build_sample_registry(rows, expected_ids=["ID01", "ID02"])


def test_sample_registry_rejects_multiple_inputs_for_one_expected_sample():
    rows = [
        {"sample_id": "ID01", "fastq": "a.fastq.gz"},
        {"sample_id": "ID01", "fastq": "b.fastq.gz"},
    ]

    with pytest.raises(SampleRegistryError, match="multiple inputs"):
        build_sample_registry(rows, expected_ids=["ID01"])


def test_sample_registry_adds_explicit_absent_deliverable_rows():
    rows = [
        {
            "deliverable_id": "DELIVERY_01",
            "sample_id": "ID01",
            "fastq": "Filtered/ID01.fq.gz",
            "manuscript_id": "ID01",
        }
    ]

    registry = build_sample_registry(rows, expected_ids=["ID01"])

    assert registry[0]["analysis_status"] == "analytical"
    absent = {row["deliverable_id"] for row in registry if row["analysis_status"] == "absent_cohort"}
    assert absent == {"DELIVERY_07", "DELIVERY_10"}
