import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.codex_reanalysis.pilot_smoke import (
    SmokeError,
    build_smoke_status,
    coverage_from_paf,
    distance_to_intervals,
    evaluate_assembly,
    evaluate_secondary_saturation,
    evaluate_split_reads,
    load_intervals,
    parse_assembly_info,
    parse_sa_segments,
    secondary_histogram,
    segment_inside,
)


# --- secondary alignment saturation -----------------------------------------


def test_a_decaying_tail_is_not_saturated():
    """Fewer reads need many hits than few: the ceiling was not binding."""
    histogram = {1: 5000, 2: 800, 18: 40, 19: 20, 20: 8}

    observation = evaluate_secondary_saturation(histogram, secondary_max=20)

    assert observation.status == "passed"
    assert observation.metrics["at_ceiling"] == 8


def test_a_spike_at_exactly_the_ceiling_is_saturated():
    """minimap2 clips at -N, so reads that wanted more all land on N."""
    histogram = {1: 5000, 2: 800, 18: 40, 19: 20, 20: 900}

    observation = evaluate_secondary_saturation(histogram, secondary_max=20)

    assert observation.status == "failed"
    assert "ceiling" in observation.detail


def test_no_read_reaching_the_ceiling_passes_without_a_ratio():
    histogram = {1: 5000, 2: 800}

    observation = evaluate_secondary_saturation(histogram, secondary_max=20)

    assert observation.status == "passed"
    assert observation.metrics["at_ceiling"] == 0


def test_reads_at_the_ceiling_with_none_below_it_are_saturated():
    """An isolated pile at N cannot be the tail of a decaying distribution."""
    histogram = {1: 5000, 20: 30}

    observation = evaluate_secondary_saturation(histogram, secondary_max=20)

    assert observation.status == "failed"


def test_more_secondary_alignments_than_the_ceiling_allows_is_an_error():
    with pytest.raises(SmokeError, match="exceeds"):
        evaluate_secondary_saturation({21: 3}, secondary_max=20)


def test_the_histogram_counts_reads_not_alignments():
    qnames = ["r1", "r1", "r1", "r2", "r3", "r3"]

    assert secondary_histogram(qnames) == {3: 1, 1: 1, 2: 1}


# --- split reads across the retention boundary ------------------------------


def test_sa_segments_are_parsed_from_the_tag():
    tag = "NC_037350.1,25600000,+,100M,60,2;NC_037329.1,6300000,-,80M,60,1;"

    assert parse_sa_segments(tag) == [
        ("NC_037350.1", 25600000),
        ("NC_037329.1", 6300000),
    ]


def test_a_malformed_sa_tag_is_an_error_not_a_silent_skip():
    with pytest.raises(SmokeError, match="SA"):
        parse_sa_segments("NC_037350.1,notanumber,+,100M,60,2;")


def test_intervals_are_loaded_from_the_mapping_bed(tmp_path):
    bed = tmp_path / "mapping_targets.bed"
    bed.write_text(
        "NC_037350.1\t25533030\t31705624\tspan_MHC\n"
        "NC_037329.1\t6228863\t6335491\tcontrol_MSTN\n",
        encoding="utf-8",
    )

    intervals = load_intervals(bed)

    assert intervals["NC_037350.1"] == [(25533030, 31705624)]
    assert segment_inside(intervals, "NC_037350.1", 25600000) is True
    assert segment_inside(intervals, "NC_037350.1", 25000000) is False
    assert segment_inside(intervals, "NC_999999.1", 1) is False


def test_split_reads_contained_by_the_retained_intervals_pass(tmp_path):
    intervals = {"NC_037350.1": [(25000000, 29000000)]}
    records = [
        ("r1", "NC_037350.1,25600000,+,100M,60,2;"),
        ("r2", "NC_037350.1,28000000,+,100M,60,2;"),
    ]

    observation = evaluate_split_reads(records, intervals, max_external_fraction=0.10)

    assert observation.status == "passed"
    assert observation.metrics["boundary_fraction"] == 0.0


def test_split_reads_cut_just_outside_the_retained_intervals_fail():
    """A segment 100 kb past the edge is local structure the window should hold."""
    intervals = {"NC_037350.1": [(25000000, 29000000)]}
    records = [("r%d" % i, "NC_037350.1,24900000,+,100M,60,2;") for i in range(3)] + [
        ("kept", "NC_037350.1,25600000,+,100M,60,2;")
    ]

    observation = evaluate_split_reads(records, intervals, max_external_fraction=0.10)

    assert observation.status == "failed"
    assert observation.metrics["boundary_fraction"] == pytest.approx(0.75)


def test_a_distant_paralog_is_reported_but_does_not_fail_the_check():
    """Reads whose other segment sits on another chromosome are ordinary
    paralogy; no target-only mapping can retain them, and counting them as a
    defect would fail the check on healthy data."""
    intervals = {"NC_037350.1": [(25000000, 29000000)]}
    records = [("r%d" % i, "NC_037329.1,6300000,+,100M,60,2;") for i in range(9)] + [
        ("kept", "NC_037350.1,25600000,+,100M,60,2;")
    ]

    observation = evaluate_split_reads(records, intervals, max_external_fraction=0.10)

    assert observation.status == "passed"
    assert observation.metrics["boundary_fraction"] == 0.0
    assert observation.metrics["distant_fraction"] == pytest.approx(0.9)


def test_a_far_segment_on_the_same_contig_is_distant_not_a_boundary_cut():
    intervals = {"NC_037350.1": [(25000000, 29000000)]}

    observation = evaluate_split_reads(
        [("r1", "NC_037350.1,1000,+,100M,60,2;")],
        intervals,
        max_external_fraction=0.10,
        adjacency_window=1_000_000,
    )

    assert observation.status == "passed"
    assert observation.metrics["reads_with_distant_segments"] == 1


def test_a_read_is_counted_once_however_many_segments_are_cut():
    intervals = {"NC_037350.1": [(25000000, 29000000)]}
    records = [
        ("r1", "NC_037350.1,24900000,+,100M,60,2;NC_037350.1,29100000,+,100M,60,2;"),
        ("r2", "NC_037350.1,25600000,+,100M,60,2;"),
    ]

    observation = evaluate_split_reads(records, intervals, max_external_fraction=0.60)

    assert observation.metrics["reads_with_supplementary"] == 2
    assert observation.metrics["boundary_fraction"] == pytest.approx(0.5)


def test_the_adjacency_window_is_measured_from_the_interval_edge():
    intervals = {"NC_037350.1": [(25000000, 29000000)]}

    assert distance_to_intervals(intervals, "NC_037350.1", 25600000) == 0
    assert distance_to_intervals(intervals, "NC_037350.1", 24900001) == 100000
    assert distance_to_intervals(intervals, "NC_037350.1", 29100000) == 100000
    assert distance_to_intervals(intervals, "NC_999999.1", 1) is None


def test_no_split_reads_at_all_is_not_evidence_of_success():
    """A target BAM with no supplementary alignments cannot demonstrate that
    structural evidence survives the retention filter."""
    observation = evaluate_split_reads([], {"x": [(0, 10)]}, max_external_fraction=0.10)

    assert observation.status == "not_run"


# --- assembly ---------------------------------------------------------------

ASSEMBLY_INFO = """#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path
contig_1\t2100000\t35\tN\tN\t1\t*\t1
contig_2\t900000\t31\tN\tN\t1\t*\t2
"""


def test_assembly_info_is_parsed(tmp_path):
    contigs = parse_assembly_info(ASSEMBLY_INFO)

    assert [c["length"] for c in contigs] == [2100000, 900000]


def test_paf_coverage_measures_breadth_and_identity():
    # query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, block, mapq
    paf = [
        "contig_1\t2100000\t0\t2100000\t+\tNC_037350.1\t3400000\t0\t2100000\t2079000\t2100000\t60",
        "contig_2\t900000\t0\t900000\t+\tNC_037350.1\t3400000\t2100000\t3000000\t891000\t900000\t60",
    ]

    breadth, identity, covered = coverage_from_paf(paf, region_span=3400000)

    assert covered == 3000000
    assert breadth == pytest.approx(3000000 / 3400000)
    assert identity == pytest.approx((2079000 + 891000) / 3000000)


def test_overlapping_paf_intervals_are_not_double_counted():
    paf = [
        "c1\t100\t0\t100\t+\tR\t1000\t0\t600\t600\t600\t60",
        "c2\t100\t0\t100\t+\tR\t1000\t400\t900\t500\t500\t60",
    ]

    _, _, covered = coverage_from_paf(paf, region_span=1000)

    assert covered == 900


def test_a_faithful_assembly_passes():
    observation = evaluate_assembly(
        contigs=[{"length": 2100000}, {"length": 900000}],
        region_span=3400000,
        breadth=0.88,
        identity=0.995,
        genes_expected=12,
        genes_recovered=12,
        thresholds={"min_breadth": 0.80, "min_identity": 0.98, "min_gene_recovery": 0.90},
    )

    assert observation.status == "passed"


def test_a_collapsed_assembly_fails_on_span():
    """Flye folding the two haplotypes together shortens the assembly."""
    observation = evaluate_assembly(
        contigs=[{"length": 900000}],
        region_span=3400000,
        breadth=0.26,
        identity=0.995,
        genes_expected=12,
        genes_recovered=4,
        thresholds={"min_breadth": 0.80, "min_identity": 0.98, "min_gene_recovery": 0.90},
    )

    assert observation.status == "failed"
    assert "breadth" in observation.detail


def test_a_chimeric_assembly_fails_on_identity():
    observation = evaluate_assembly(
        contigs=[{"length": 3400000}],
        region_span=3400000,
        breadth=0.95,
        identity=0.91,
        genes_expected=12,
        genes_recovered=12,
        thresholds={"min_breadth": 0.80, "min_identity": 0.98, "min_gene_recovery": 0.90},
    )

    assert observation.status == "failed"
    assert "identity" in observation.detail


def test_missing_expected_genes_fail_even_when_the_span_looks_right():
    observation = evaluate_assembly(
        contigs=[{"length": 3400000}],
        region_span=3400000,
        breadth=0.95,
        identity=0.995,
        genes_expected=12,
        genes_recovered=7,
        thresholds={"min_breadth": 0.80, "min_identity": 0.98, "min_gene_recovery": 0.90},
    )

    assert observation.status == "failed"
    assert "gene" in observation.detail


def test_excluded_footprint_is_dropped_from_identity_not_breadth():
    """A paralog-cluster block scoring poorly must not sink identity, but a
    dropout inside that same footprint is still a real dropout for breadth."""
    paf = [
        # clean block outside the excluded footprint
        "c1\t1000000\t0\t1000000\t+\tR\t3400000\t0\t1000000\t995000\t1000000\t60",
        # poor block entirely inside the excluded footprint [1000000, 1500000)
        "c2\t500000\t0\t500000\t+\tR\t3400000\t1000000\t1500000\t300000\t500000\t60",
    ]

    breadth, identity, covered = coverage_from_paf(
        paf, region_span=3400000, exclude=[(1000000, 1500000)]
    )

    assert covered == 1500000
    assert breadth == pytest.approx(1500000 / 3400000)
    assert identity == pytest.approx(995000 / 1000000)


def test_a_block_straddling_the_excluded_boundary_is_partially_kept():
    paf = ["c1\t100000\t0\t100000\t+\tR\t3400000\t900000\t1000000\t80000\t100000\t60"]

    _, identity, _ = coverage_from_paf(paf, region_span=3400000, exclude=[(950000, 1000000)])

    # half the block (50000 of 100000 bp) falls inside the excluded region
    assert identity == pytest.approx((80000 * 0.5) / (100000 * 0.5))


def test_evaluate_assembly_names_the_excluded_genes_for_audit():
    observation = evaluate_assembly(
        contigs=[{"length": 3400000}],
        region_span=3400000,
        breadth=0.95,
        identity=0.995,
        genes_expected=10,
        genes_recovered=10,
        thresholds={"min_breadth": 0.80, "min_identity": 0.90, "min_gene_recovery": 0.90},
        excluded_genes=["BTN_LOC104969893", "BTN_LOC112443860"],
    )

    assert observation.status == "passed"
    assert observation.metrics["excluded_paralog_genes"] == [
        "BTN_LOC104969893",
        "BTN_LOC112443860",
    ]


def test_an_assembly_with_no_contigs_is_a_failure_not_a_perfect_score():
    observation = evaluate_assembly(
        contigs=[],
        region_span=3400000,
        breadth=0.0,
        identity=0.0,
        genes_expected=12,
        genes_recovered=0,
        thresholds={"min_breadth": 0.80, "min_identity": 0.98, "min_gene_recovery": 0.90},
    )

    assert observation.status == "failed"


# --- status document --------------------------------------------------------


def test_the_status_names_the_checks_the_gate_consumes(tmp_path):
    from scripts.codex_reanalysis.pilot_smoke import Observation

    status = build_smoke_status(
        sample_id="ID02",
        bam_sha256="a" * 64,
        region="NC_037350.1:25500000-28900000",
        secondary_max=20,
        observations=[
            Observation("secondary_saturation", "passed", "tail decays", {}),
            Observation("sv_smoke", "passed", "all segments retained", {}),
            Observation("assembly_smoke", "passed", "faithful", {}),
        ],
    )

    assert status["status"] == "passed"
    assert {c["name"] for c in status["checks"]} == {
        "secondary_saturation",
        "sv_smoke",
        "assembly_smoke",
    }
    assert status["bam_sha256"] == "a" * 64
    json.dumps(status)


def test_one_failed_check_fails_the_whole_smoke():
    from scripts.codex_reanalysis.pilot_smoke import Observation

    status = build_smoke_status(
        sample_id="s",
        bam_sha256="b" * 64,
        region="r",
        secondary_max=20,
        observations=[
            Observation("secondary_saturation", "passed", "", {}),
            Observation("sv_smoke", "failed", "", {}),
            Observation("assembly_smoke", "passed", "", {}),
        ],
    )

    assert status["status"] == "failed"


def test_a_check_that_did_not_run_leaves_the_smoke_incomplete():
    from scripts.codex_reanalysis.pilot_smoke import Observation

    status = build_smoke_status(
        sample_id="s",
        bam_sha256="c" * 64,
        region="r",
        secondary_max=20,
        observations=[
            Observation("secondary_saturation", "passed", "", {}),
            Observation("sv_smoke", "not_run", "", {}),
            Observation("assembly_smoke", "passed", "", {}),
        ],
    )

    assert status["status"] == "incomplete"


def test_a_missing_tool_becomes_a_smoke_failure_not_an_escaping_oserror():
    """An OSError would escape before any status document is written, leaving
    the operator with an exit code and no record of what blocked the pilot."""
    from scripts.codex_reanalysis.pilot_smoke import _run

    with pytest.raises(SmokeError, match="could not be run"):
        _run(["definitely-not-a-real-binary-xyz"], label="flye")


def test_pilot_smoke_runs_as_a_direct_command():
    """The remote runner invokes this as a script, not as a module."""
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(
        [sys.executable, str(root / "scripts" / "codex_reanalysis" / "pilot_smoke.py"), "--help"],
        cwd=root,
        capture_output=True,
        text=True,
        # The help text carries the docstring, which is not representable in the
        # Windows default console encoding.
        encoding="utf-8",
        check=False,
    )

    assert result.returncode == 0, result.stderr
    for option in ("--bam", "--mapping-bed", "--analytical-bed", "--region", "--secondary-max"):
        assert option in result.stdout


def test_a_status_missing_a_required_check_is_rejected():
    from scripts.codex_reanalysis.pilot_smoke import Observation

    with pytest.raises(SmokeError, match="missing"):
        build_smoke_status(
            sample_id="s",
            bam_sha256="d" * 64,
            region="r",
            secondary_max=20,
            observations=[Observation("secondary_saturation", "passed", "", {})],
        )
