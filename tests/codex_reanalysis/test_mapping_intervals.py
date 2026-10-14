import pytest

from scripts.codex_reanalysis.coordinate_registry import (
    CoordinateError,
    build_mapping_intervals,
)


def record(name, contig, start, end, category, target_type="reference_feature"):
    return {
        "target_name": name,
        "contig_accession": contig,
        "bed_start": start,
        "bed_end": end,
        "category": category,
        "target_type": target_type,
    }


LENGTHS = {"NC_1": 40_000_000, "NC_2": 1_000_000}


def test_each_gene_gains_a_flank_on_both_sides():
    intervals = build_mapping_intervals(
        [record("g1", "NC_1", 1_000_000, 1_010_000, "control_candidate")],
        flank=50_000,
        contiguous_categories=set(),
        contig_lengths=LENGTHS,
    )

    assert intervals == [("NC_1", 950_000, 1_060_000, "g1")]


def test_a_flank_is_clamped_to_the_contig_bounds():
    intervals = build_mapping_intervals(
        [
            record("start", "NC_2", 1_000, 2_000, "control_candidate"),
            record("end", "NC_2", 990_000, 999_000, "control_candidate"),
        ],
        flank=50_000,
        contiguous_categories=set(),
        contig_lengths=LENGTHS,
    )

    assert intervals[0][1] == 0
    assert intervals[-1][2] == 1_000_000


def test_overlapping_flanked_genes_are_merged():
    intervals = build_mapping_intervals(
        [
            record("a", "NC_1", 1_000_000, 1_010_000, "control_candidate"),
            record("b", "NC_1", 1_050_000, 1_060_000, "control_candidate"),
        ],
        flank=50_000,
        contiguous_categories=set(),
        contig_lengths=LENGTHS,
    )

    assert len(intervals) == 1
    assert intervals[0][1] == 950_000
    assert intervals[0][2] == 1_110_000


def test_a_contiguous_category_fills_the_gap_between_its_genes():
    """Gene bodies alone leave the intergenic space unmapped, which is where
    local assembly and same-contig determination actually need reads."""
    genes = [
        record("far_left", "NC_1", 25_600_000, 25_610_000, "class_II"),
        record("far_right", "NC_1", 25_830_000, 25_840_000, "class_II"),
    ]

    without = build_mapping_intervals(genes, 10_000, set(), LENGTHS)
    with_span = build_mapping_intervals(genes, 10_000, {"class_II"}, LENGTHS)

    assert len(without) == 2
    assert len(with_span) == 1
    assert with_span[0][1] == 25_590_000
    assert with_span[0][2] == 25_850_000


def test_categories_outside_the_contiguous_set_are_not_spanned():
    """Single-copy controls sit on different chromosomes and must not be joined."""
    genes = [
        record("c1", "NC_1", 1_000_000, 1_010_000, "control_candidate"),
        record("c2", "NC_1", 30_000_000, 30_010_000, "control_candidate"),
    ]

    intervals = build_mapping_intervals(genes, 50_000, {"class_I"}, LENGTHS)

    assert len(intervals) == 2


def test_coordinate_free_probes_are_not_given_intervals():
    genes = [
        record("real", "NC_1", 1_000_000, 1_010_000, "class_I"),
        record("probe", "", None, None, "class_I", target_type="haplotype_probe"),
    ]

    intervals = build_mapping_intervals(genes, 50_000, {"class_I"}, LENGTHS)

    assert len(intervals) == 1
    assert "probe" not in intervals[0][3]
    assert "real" in intervals[0][3]


def test_an_unknown_contig_length_is_rejected():
    with pytest.raises(CoordinateError, match="contig length"):
        build_mapping_intervals(
            [record("g", "NC_missing", 1_000, 2_000, "class_I")],
            flank=50_000,
            contiguous_categories=set(),
            contig_lengths=LENGTHS,
        )


def test_intervals_are_sorted_and_non_overlapping():
    genes = [
        record("b", "NC_1", 3_000_000, 3_010_000, "class_I"),
        record("a", "NC_1", 1_000_000, 1_010_000, "class_I"),
        record("c", "NC_2", 100_000, 110_000, "control_candidate"),
    ]

    intervals = build_mapping_intervals(genes, 10_000, set(), LENGTHS)

    assert [i[0] for i in intervals] == ["NC_1", "NC_1", "NC_2"]
    assert intervals[0][2] <= intervals[1][1]
