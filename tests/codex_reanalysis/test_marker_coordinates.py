from __future__ import annotations

import pytest

from scripts.codex_reanalysis.resolve_marker_coordinates import (
    MarkerCoordinateError,
    annotate_reference_alleles,
    build_faidx_command,
    decide_resolution,
    hash_sequence,
    identity_from_comparison_rows,
    parse_fasta_sequence,
    parse_faidx_base,
    resolve_markers_from_identity,
)


def test_build_faidx_command_requests_one_contig():
    assert build_faidx_command("ref.fa", "chr23", samtools="samtools") == [
        "samtools",
        "faidx",
        "ref.fa",
        "chr23",
    ]


def test_parse_fasta_sequence_joins_wrapped_lines_uppercase():
    text = ">chr23\nacgt\nACGT\n"
    assert parse_fasta_sequence(text) == "ACGTACGT"


def test_parse_fasta_sequence_rejects_non_fasta_output():
    with pytest.raises(MarkerCoordinateError):
        parse_fasta_sequence("")


def test_hash_sequence_is_deterministic():
    assert hash_sequence("ACGT") == hash_sequence("ACGT")
    assert hash_sequence("ACGT") != hash_sequence("ACGG")


def test_decide_resolution_accepts_matching_hashes():
    assert decide_resolution("abc", "abc") == "verified_identical_chromosome"


def test_decide_resolution_flags_mismatched_hashes_as_decision_required():
    assert decide_resolution("abc", "def") == "decision_required"


def test_resolve_markers_reuses_position_for_identical_chromosome():
    markers = [
        {
            "marker_id": "rs_test",
            "ucd1_2_contig": "chr23",
            "ucd1_2_position": "27603651",
            "ucd2_0_contig": "NC_037350.1",
            "citation": "test",
            "notes": "",
        }
    ]
    checked = {
        "chr23->NC_037350.1": {"ucd1_2_sha256": "same", "ucd2_0_sha256": "same"},
    }
    resolved = resolve_markers_from_identity(markers, checked)
    assert resolved == [
        {
            "marker_id": "rs_test",
            "ucd1_2_contig": "chr23",
            "ucd1_2_position": 27603651,
            "ucd2_0_contig": "NC_037350.1",
            "ucd2_0_position": 27603651,
            "resolution_status": "verified_identical_chromosome",
            "citation": "test",
            "notes": "",
        }
    ]


def test_resolve_markers_leaves_position_unset_on_mismatch():
    markers = [
        {
            "marker_id": "rs_test",
            "ucd1_2_contig": "chr23",
            "ucd1_2_position": "1000",
            "ucd2_0_contig": "NC_037350.1",
            "citation": "",
            "notes": "",
        }
    ]
    checked = {
        "chr23->NC_037350.1": {"ucd1_2_sha256": "aaa", "ucd2_0_sha256": "bbb"},
    }
    resolved = resolve_markers_from_identity(markers, checked)
    assert resolved[0]["resolution_status"] == "decision_required"
    assert resolved[0]["ucd2_0_position"] is None


def test_comparison_identical_hashes_reuse_the_legacy_position():
    comparison = [
        {
            "old_accession": "NC_037349.1",
            "new_accession": "NC_037349.1",
            "old_fasta_name": "chr22",
            "old_sha256": "abc",
            "new_sha256": "abc",
            "status": "identical",
        }
    ]
    checked = identity_from_comparison_rows(comparison)
    markers = [
        {
            "marker_id": "rs110616206",
            "ucd1_2_contig": "chr22",
            "ucd1_2_position": "27162840",
            "ucd2_0_contig": "NC_037349.1",
            "citation": "",
            "notes": "",
        }
    ]
    resolved = resolve_markers_from_identity(markers, checked)
    assert resolved[0]["ucd2_0_position"] == 27162840
    assert resolved[0]["resolution_status"] == "verified_identical_chromosome"


def test_comparison_status_cannot_override_hash_mismatch():
    comparison = [
        {
            "old_accession": "NC_037350.1",
            "new_accession": "NC_037350.1",
            "old_fasta_name": "chr23",
            "old_sha256": "aaa",
            "new_sha256": "bbb",
            "status": "identical",
        }
    ]
    checked = identity_from_comparison_rows(comparison)
    markers = [
        {
            "marker_id": "rs29026690",
            "ucd1_2_contig": "chr23",
            "ucd1_2_position": "27603651",
            "ucd2_0_contig": "NC_037350.1",
            "citation": "",
            "notes": "",
        }
    ]
    resolved = resolve_markers_from_identity(markers, checked)
    assert resolved[0]["resolution_status"] == "decision_required"
    assert resolved[0]["ucd2_0_position"] is None


def test_resolve_markers_requires_a_chromosome_check():
    markers = [
        {
            "marker_id": "rs_test",
            "ucd1_2_contig": "chr23",
            "ucd1_2_position": "1",
            "ucd2_0_contig": "NC_037350.1",
            "citation": "",
            "notes": "",
        }
    ]
    with pytest.raises(MarkerCoordinateError, match="no chromosome comparison"):
        resolve_markers_from_identity(markers, {})


def test_parse_faidx_base_requires_exactly_one_base():
    assert parse_faidx_base(">NC_037350.1:27603651-27603651\nT\n") == "T"
    with pytest.raises(MarkerCoordinateError, match="one base"):
        parse_faidx_base(">NC_037350.1:1-2\nAT\n")


def test_annotate_reference_alleles_keeps_unresolved_empty():
    resolved = [
        {
            "marker_id": "rs_test",
            "ucd2_0_contig": "NC_037350.1",
            "ucd2_0_position": None,
            "resolution_status": "decision_required",
        }
    ]
    annotated = annotate_reference_alleles(resolved, lambda contig, position: "T")
    assert annotated[0]["ucd2_0_reference_allele"] == ""
    annotated_ok = annotate_reference_alleles(
        [
            {
                "marker_id": "rs29026690",
                "ucd2_0_contig": "NC_037350.1",
                "ucd2_0_position": 27603651,
                "resolution_status": "verified_identical_chromosome",
            }
        ],
        lambda contig, position: "T",
    )
    assert annotated_ok[0]["ucd2_0_reference_allele"] == "T"
