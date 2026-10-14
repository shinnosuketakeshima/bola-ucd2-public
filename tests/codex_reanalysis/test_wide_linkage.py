from scripts.codex_reanalysis.wide_linkage import (
    MapHit,
    classify_drb3_resolution,
    extract_drb3_panel_sequences,
    summarize_linkage,
)


PANEL = """>BoLA-DRB3*016:01
ACGT
>BoLA-DRB3*015:01
TGCA
>BoLA-DQA*026:03
NNNN
"""


def test_extract_drb3_panel_sequences_keeps_both_haplotypes():
    fasta = extract_drb3_panel_sequences(PANEL, "*015:01", "*016:01")
    assert "BoLA-DRB3*016:01" in fasta
    assert "BoLA-DRB3*015:01" in fasta
    assert "BoLA-DQA" not in fasta


def test_homozygous_control():
    resolution, note = classify_drb3_resolution(
        "*010:01",
        "*010:01",
        [MapHit("BoLA-DRB3*010:01", "contig_2", 100, 60)],
    )
    assert resolution == "homozygous_control"
    assert "mapq60" in note


def test_collapsed_identical_position():
    hits = [
        MapHit("BoLA-DRB3*002:01", "contig_13", 500, 60),
        MapHit("BoLA-DRB3*015:01", "contig_13", 500, 60),
    ]
    resolution, note = classify_drb3_resolution("*002:01", "*015:01", hits)
    assert resolution == "collapsed"
    assert "identical position" in note


def test_offset_artifact():
    hits = [
        MapHit("BoLA-DRB3*016:01", "contig_12", 1000, 57),
        MapHit("BoLA-DRB3*015:01", "contig_12", 3463, 60),
    ]
    resolution, note = classify_drb3_resolution("*015:01", "*016:01", hits)
    assert resolution == "collapsed"
    assert "offset artifact" in note


def test_summarize_linkage_marks_different_contigs():
    row = summarize_linkage(
        "*015:01",
        "*016:01",
        [MapHit("BoLA-DRB3*016:01", "contig_12", 1000, 57)],
        [MapHit("rs29026690_flank", "contig_25", 200, 60)],
    )
    assert row["same_contig"] == "no"
    assert row["drb3_contig"] == "contig_12"
    assert row["rs29026690_contig"] == "contig_25"
