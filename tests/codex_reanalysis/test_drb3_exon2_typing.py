from __future__ import annotations

import pytest

from scripts.codex_reanalysis.drb3_exon2_typing import (
    Drb3Exon2Error,
    diplotype_costs,
    exon_from_gff,
    expected_cost,
    group_identical,
    group_indices,
    normalize_allele,
    query_segment,
    rank_diplotypes,
)


def test_normalize_allele_reduces_to_two_fields():
    assert normalize_allele("BoLA-DRB3*015:01:02") == "*015:01"
    assert normalize_allele("*016:01") == "*016:01"
    with pytest.raises(Drb3Exon2Error):
        normalize_allele("BoLA-DQA*001:01")


def test_exon_from_gff_takes_second_exon_in_transcript_order(tmp_path):
    gff = tmp_path / "g.gff"
    gff.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tx\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=BOLA-DRB3",
                "chr1\tx\tmRNA\t100\t900\t.\t+\t.\tID=rna-1;Parent=gene-A",
                "chr1\tx\texon\t100\t200\t.\t+\t.\tID=e1;Parent=rna-1",
                "chr1\tx\texon\t400\t670\t.\t+\t.\tID=e2;Parent=rna-1",
                "chr1\tx\texon\t800\t900\t.\t+\t.\tID=e3;Parent=rna-1",
            ]
        ),
        encoding="utf-8",
    )
    assert exon_from_gff(gff) == ("chr1", 400, 670, "+")
    with pytest.raises(Drb3Exon2Error):
        exon_from_gff(gff, gene_name="BOLA-DQA1")


def test_query_segment_spans_interval_through_indels_and_clips():
    # 5 soft-clipped bases, 10M, 2D, 3M, 2I, 5M ; POS=101
    seq = "SSSSS" + "AAAAAAAAAA" + "CCC" + "II" + "GGGGG"
    cigar = "5S10M2D3M2I5M"
    # reference 100..110 is the 10M block (0-based), interval [102, 108)
    assert query_segment(101, cigar, seq, 102, 108) == "AAAAAA"
    # an interval spanning the deletion and the 3M block keeps the aligned bases
    assert query_segment(101, cigar, seq, 108, 115) == seq[13:18]


def test_query_segment_rejects_reads_that_do_not_span():
    assert query_segment(101, "10M", "A" * 10, 95, 105) is None
    assert query_segment(101, "10M", "A" * 10, 105, 115) is None


def test_group_identical_collapses_equal_sequences():
    seqs, names = group_identical(
        {"BoLA-DRB3*001:01": "AAA", "BoLA-DRB3*002:01": "AAA", "BoLA-DRB3*003:01:01": "CCC"}
    )
    assert len(seqs) == 2
    assert ["*001:01", "*002:01"] in names
    assert group_indices("*003:01", names) == [names.index(["*003:01"])]


def test_diplotype_cost_prefers_the_pair_that_explains_both_haplotypes():
    # groups 0,1,2; reads 0-2 look like group 0, reads 3-5 like group 1
    distances = [[0, 9, 9]] * 3 + [[9, 0, 9]] * 3
    costs = diplotype_costs(distances, 3)
    best, tied, second = rank_diplotypes(costs)
    assert best == 0 and tied == [(0, 1)]
    assert second == 9 * 3 and costs[(0, 0)] == 27
    assert expected_cost(costs, [0], [1]) == 0
    assert expected_cost(costs, [0], [0]) == 27


def test_tied_pairs_are_all_reported():
    distances = [[1, 1], [1, 1]]
    best, tied, second = rank_diplotypes(diplotype_costs(distances, 2))
    assert best == 2 and len(tied) == 3 and second is None


def test_panel_window_is_the_intersection_of_allele_placements():
    pytest.importorskip("edlib")
    from scripts.codex_reanalysis.drb3_exon2_typing import cut_to_exon, panel_window

    ref = "ACGTTGCAAGGCTTAACCGGATCCGATTACAGGCTAGCTAGGATCCA" * 2
    panel = {"BoLA-DRB3*001:01": ref, "BoLA-DRB3*002:01": ref[10:-8], "BoLA-DRB3*003:01": ref[4:-15]}
    start, end, rows = panel_window(cut_to_exon(panel, ref), ref)
    assert (start, end) == (10, len(ref) - 15)
    assert len(rows) == 3


def test_distance_matrix_and_trim_use_edit_distance():
    pytest.importorskip("edlib")
    from scripts.codex_reanalysis.drb3_exon2_typing import distance_matrix, trim_to

    window = "GATTACAGATTACA"
    trimmed, edits = trim_to(window, "TTT" + window + "CCC")
    assert trimmed == window and edits == 0
    assert distance_matrix([window], [window, "GATTACAGATTAGA"]) == [[0, 1]]
