from pathlib import Path

import pytest

from scripts.codex_reanalysis.attribute_dq_contigs import (
    Alignment,
    Anchor,
    AttributionError,
    attribute_feature,
    chain_anchors,
    count_spanning,
    dqa2_read_origin,
    find_events,
    load_features,
    parse_blast_hsps,
    parse_paf_blocks,
    parse_sam_alignments,
    parse_sam_spans,
    parse_source,
    project,
    supporting_reads,
)

FEATURES = {
    "BOLA-DQA2": (25583031, 25589201),
    "BOLA-DQB": (25607507, 25620857),
    "BOLA-DQA5": (25631958, 25643878),
    "BOLA-DQB1": (25672687, 25680507),
    "LOC100848815": (25691259, 25695296),
    "BOLA-DRB3": (25723717, 25734819),
}


def anchor(q1, q2, r1, r2, strand="+"):
    return Anchor("c", q1, q2, r1, r2, strand, 98.0)


def test_blast_hsps_are_shifted_to_chromosome_and_filtered():
    text = "1\t5000\t1\t5000\t98.0\t5000\n10\t900\t1\t890\t99.0\t890\n1\t3000\t9000\t6001\t90.0\t3000\n"
    anchors = parse_blast_hsps(text, "c", 25_400_000, 95.0, 2000)
    assert anchors == [Anchor("c", 1, 5000, 25_400_001, 25_405_000, "+", 98.0)]


def test_paf_blocks_become_one_based_closed():
    line = "c\t200000\t0\t9000\t+\tw\t600000\t100\t9100\t8000\t9000\t60\n"
    blocks = parse_paf_blocks(line, 25_400_000, 5000)
    assert blocks["c"][0].q_start == 1 and blocks["c"][0].q_end == 9000
    assert blocks["c"][0].r_start == 25_400_101 and blocks["c"][0].r_end == 25_409_100


def test_gap_dominated_paf_blocks_are_dropped():
    line = "c\t200000\t73930\t99472\t+\tw\t600000\t193252\t279536\t22000\t85000\t60\n"
    assert parse_paf_blocks(line, 25_400_000, 5000, 70.0) == {}


def test_dqa2_read_origin_counts_reads_on_placed_hit():
    hits = [{"contig": "c3", "contig_start": "25049", "contig_end": "30593"}]
    alignments = [("c3", 20000, 40000, 60, 0), ("c3", 90000, 95000, 60, 0), ("c9", 1, 5000, 60, 0),
                  ("c3", 20000, 40000, 5, 0)]
    assert dqa2_read_origin(alignments, hits, 2000) == {"reads_on_dqa2_placed_hit": 1, "reads_elsewhere": 2}


def test_chain_skips_off_diagonal_repeat_hits():
    chain = chain_anchors(
        [
            anchor(1, 90000, 25490000, 25580000),
            anchor(40000, 42000, 25900000, 25902000),  # interspersed repeat elsewhere
            anchor(95000, 180000, 25658000, 25743000),
        ]
    )
    assert [(a.q_start, a.r_start) for a in chain] == [(1, 25490000), (95000, 25658000)]


def test_chain_follows_minus_strand_contigs():
    chain = chain_anchors([anchor(1, 5000, 25700000, 25705000, "-"), anchor(6000, 9000, 25690000, 25693000, "-")])
    assert len(chain) == 2


def test_hit_inside_anchor_is_projected_and_named():
    chain = [anchor(114000, 180000, 25677700, 25743700)]
    placed = project(chain, 128753, 132279)
    assert placed["placement"] == "inside_anchor"
    assert attribute_feature(placed["reference_position"], None, None, FEATURES, 2000) == "DQA1 (LOC100848815)"


def test_minus_strand_projection_runs_backwards():
    placed = project([anchor(1, 10000, 25580000, 25590000, "-")], 4000, 4000)
    assert placed["reference_position"] == 25590000 - 3999


def test_hit_between_anchors_is_only_bracketed():
    chain = [anchor(1, 90000, 25480000, 25571270), anchor(190000, 200000, 25699670, 25709670)]
    placed = project(chain, 128000, 133000)
    assert placed["placement"] == "between_anchors"
    label = attribute_feature(None, 25571270, 25699670, FEATURES, 2000)
    assert label.startswith("bracketed:") and "DQA2" in label and "DQA1 (LOC100848815)" in label


def test_deletion_lists_features_inside_the_skipped_reference():
    chain = [anchor(1, 91774, 25490102, 25581631), anchor(96310, 179542, 25658061, 25742714)]
    (event,) = find_events(chain, FEATURES, 10000, "deletion")
    assert event["reference_gap_bp"] == 76429
    assert event["contig_gap_bp"] == 4535
    assert event["features_between"] == "DQA2|DQB (reference BOLA-DQB)|DQA5"
    assert find_events(chain, FEATURES, 10000, "insertion") == []


def test_insertion_is_contig_sequence_absent_from_reference():
    chain = [anchor(1, 72164, 25504869, 25576670), anchor(97927, 116539, 25576909, 25594622)]
    (event,) = find_events(chain, FEATURES, 10000, "insertion")
    assert event["net_bp"] == 97927 - 72164 - 1 - (25576909 - 25576670 - 1)


def aln(start, cigar, target="c", mapq=60, flag=0, read_id="r"):
    end = start + sum(int(n) for n, op in __import__("re").findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "MDN=X") - 1
    return Alignment(read_id, target, start, end, mapq, flag, cigar)


def test_spanning_reads_need_flanks_mapq_and_primary():
    alignments = [
        aln(80000, "30001M"),  # spans
        aln(90000, "20001M"),  # left flank too short
        aln(80000, "30001M", mapq=5),  # low MAPQ
        aln(80000, "30001M", flag=0x800),  # supplementary
        aln(80000, "30001M", target="d"),  # other contig
    ]
    assert count_spanning(alignments, "c", 91774, 96310, 2000) == 1


def test_large_deletion_across_the_diagnostic_is_not_support():
    # end points enclose 5,000-9,000 with 2 kb flanks, but no query base covers it
    skipped = aln(1, "2000M10000D2000M")
    assert skipped.end == 14000
    assert count_spanning([skipped], "c", 5000, 9000, 2000) == 0
    (row,) = supporting_reads([skipped], "c", 5000, 9000, 2000, (50, 200), 0.8)
    assert row["longest_deletion"] == 10000 and row["core_aligned"] == 0
    assert not row["continuous_gap50"] and not row["continuous_gap200"]


def test_small_ont_indels_still_count():
    cigar = "".join(["300M", "3D", "200M", "2I", "400M", "12D"] * 30) + "2000M"
    read = aln(1, cigar)
    assert count_spanning([read], "c", 10000, 14000, 2000) == 1


def test_large_insertion_at_a_point_diagnostic_is_not_support():
    # continuous with the reference on both sides of 5,000 but carrying 3 kb of extra sequence there
    read = aln(1, "4999M3000I5001M")
    assert count_spanning([read], "c", 5000, 5000, 2000) == 0
    assert count_spanning([read], "c", 5000, 5000, 2000, max_gap=5000) == 0  # still at the boundary
    assert count_spanning([read], "c", 5000, 5000, 2000, max_gap=5000, boundary_max_gap=5000) == 1


def test_tolerated_gap_straddling_a_point_diagnostic_is_not_support():
    read = aln(1, "4920M150D5000M")  # 150-bp deletion over 4,921-5,070 carries the read across 5,000
    assert count_spanning([read], "c", 5000, 5000, 2000) == 0


def test_indel_away_from_the_boundaries_counts_only_at_the_tolerant_gap():
    # a 91-bp deletion about 1 kb inside a 4.5-kb junction interval (a sequence difference, not a skip)
    read = aln(1, "12575M91D6000M")
    (row,) = supporting_reads([read], "c", 10001, 14537, 2000, (200, 50), 0.8)
    assert row["longest_gap_at_boundary"] == 0
    assert row["continuous_gap200"] and not row["continuous_gap50"]


def test_gap_in_a_flank_beyond_tolerance_is_counted_only_at_the_wider_gap():
    read = aln(1, "3500M120D6500M")  # 120-bp deletion inside the left flank of 5,000-6,000
    (row,) = supporting_reads([read], "c", 5000, 6000, 2000, (50, 200), 0.8)
    assert not row["continuous_gap50"] and row["continuous_gap200"]


def test_supporting_rows_carry_read_identity_and_clipped_cigar():
    read = aln(101, "500S2000M5D2000M", read_id="read-1", flag=0x10)
    (row,) = supporting_reads([read], "c", 2400, 2405, 1000, (50,), 0.8)
    assert (row["read_id"], row["strand"], row["mapq"]) == ("read-1", "-", 60)
    assert (row["window_left"], row["window_right"]) == (1400, 3405)
    assert row["window_cigar"] == "701M5D1300M"
    assert row["longest_deletion"] == 5 and row["continuous_gap50"]


def test_sam_spans_use_reference_consuming_operations():
    sam = "@HD\tVN:1.6\nr1\t0\tc\t100\t60\t10S50M5I20M3D10M\t*\t0\t0\tA\t*\n" "r2\t4\t*\t0\t0\t*\t*\t0\t0\tA\t*\n"
    assert parse_sam_spans(sam) == [("c", 100, 100 + 50 + 20 + 3 + 10 - 1, 60, 0)]
    (record,) = parse_sam_alignments(sam)
    assert (record.read_id, record.cigar, record.strand) == ("r1", "10S50M5I20M3D10M", "+")


def test_source_needs_sample_and_contig():
    assert parse_source("ID11:contig_22") == ("ID11", "contig_22")
    with pytest.raises(AttributionError):
        parse_source("ID11")


def test_features_must_all_be_annotated(tmp_path: Path):
    gff = tmp_path / "a.gff"
    gff.write_text("NC_037350.1\tRefSeq\tgene\t25583031\t25589201\t.\t-\t.\tID=g;Name=BOLA-DQA2\n", encoding="utf-8")
    with pytest.raises(AttributionError):
        load_features(gff)


def test_window_identity_needs_eqx_cigar():
    plain = aln(1, "10000M")
    (row,) = supporting_reads([plain], "c", 5000, 5000, 2000, (200,), 0.8)
    assert row["window_identity"] == ""
    eqx = aln(1, "4000=10X5986=")  # 10 mismatches inside the 4,001-bp window 3,000-7,000
    (row,) = supporting_reads([eqx], "c", 5000, 5000, 2000, (200,), 0.8)
    assert row["window_identity"] == round(3991 / 4001, 4)
