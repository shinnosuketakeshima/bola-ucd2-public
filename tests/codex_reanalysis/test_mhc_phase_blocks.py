from __future__ import annotations

from scripts.codex_reanalysis.mhc_phase_blocks import (
    BlockBoundary,
    PhaseBlock,
    PhasedHet,
    build_blocks,
    build_boundaries,
    compartment_at,
    landmark_same_block,
    parse_phased_hets,
    window_hotspots,
)


def test_compartment_at_uses_corrected_windows():
    assert compartment_at(25_720_000) == "class_II_DQ_DR"
    assert compartment_at(26_250_000) == "btn"
    assert compartment_at(27_600_000) == "class_III"
    assert compartment_at(28_600_000) == "class_I"


def test_build_blocks_and_boundaries():
    hets = [
        PhasedHet(25_700_000, "0|1", "A"),
        PhasedHet(25_710_000, "1|0", "A"),
        PhasedHet(26_300_000, "0|1", "B"),
        PhasedHet(26_310_000, "1|0", "B"),
    ]
    blocks = build_blocks("s1", hets)
    assert len(blocks) == 2
    assert blocks[0].start == 25_700_000 and blocks[0].end == 25_710_000
    boundaries = build_boundaries("s1", blocks)
    assert len(boundaries) == 1
    assert boundaries[0].midpoint == 25_710_000 + (26_300_000 - 25_710_000) // 2
    assert "between_" in boundaries[0].boundary_type


def test_landmark_same_block_requires_span():
    blocks = [
        PhaseBlock("s", "A", 25_700_000, 25_750_000, 10),
        PhaseBlock("s", "B", 27_500_000, 27_700_000, 10),
    ]
    assert landmark_same_block(blocks, 25_728_000, 27_603_651) is False
    spanning = [PhaseBlock("s", "C", 25_700_000, 27_700_000, 50)]
    assert landmark_same_block(spanning, 25_728_000, 27_603_651) is True


def test_window_hotspots_bins_midpoints():
    boundaries = [
        BlockBoundary("s", "a", "b", 0, 0, 0, 26_250_100, "x", "y", "between_x__y"),
        BlockBoundary("s", "a", "b", 0, 0, 0, 26_260_000, "x", "y", "between_x__y"),
        BlockBoundary("s", "a", "b", 0, 0, 0, 25_650_000, "x", "x", "within_x"),
    ]
    rows = window_hotspots(boundaries, win=100_000)
    assert rows[0]["window_start"] == 26_200_000
    assert rows[0]["n_breaks"] == 2


def test_parse_phased_hets_filters(tmp_path):
    vcf = tmp_path / "x.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "NC_037350.1\t25700000\t.\tA\tG\t40\tPASS\t.\tGT:PS\t0|1:100\n"
        "NC_037350.1\t25700001\t.\tA\tG\t40\tPASS\t.\tGT:PS\t0/1:100\n"  # unphased
        "NC_037350.1\t25700002\t.\tA\tG\t40\tPASS\t.\tGT:PS\t0|0:100\n"  # hom
        ,
        encoding="utf-8",
    )
    hets = parse_phased_hets(vcf, start=25_000_000, end=29_000_000)
    assert len(hets) == 1
    assert hets[0].ps == "100"
