from __future__ import annotations

import pytest

from scripts.codex_reanalysis.assemble_interval import (
    CLASS_I_GENOME_SIZE,
    CLASS_I_REGION,
    build_extract_view_argv,
    build_flye_argv,
    build_ipd_map_argv,
)
from scripts.codex_reanalysis.assembly_qc import (
    COPY_NUMBER_CAVEAT,
    render_report,
    summarize_sample,
)
from scripts.codex_reanalysis.probe_assignment import (
    DUPLICATED_ON_CONTIG,
    HAPLOTYPE_SPLIT,
    HIGH_CONFIDENCE,
    NOT_RECOVERED,
    RECOVERED,
    UNRESOLVED,
    assign_probe,
    cluster_alignments,
    classify_locus_structure,
    parse_paf,
)
from scripts.codex_reanalysis.summarize_mapping_qc import CohortError, validate_cohort_ids


PAF_FIXTURE = """\
BoLA-2*012:01|BoLA-2	1098	0	1098	+	contig_1	50000	1000	2200	1090	1100	60	tp:A:P	dv:f:0.01
BoLA-2*026:03|BoLA-2	1098	0	1090	+	contig_1	50000	1500	2500	1000	1100	60	tp:A:S	dv:f:0.02
BoLA-3*004:03|BoLA-3	1098	0	1098	+	contig_2	40000	500	1600	1080	1100	60	tp:A:P	dv:f:0.005
BoLA-1*007:01|BoLA-1	1098	0	800	+	contig_3	30000	100	900	700	800	5	tp:A:P	dv:f:0.20
BoLA-6*013:01|BoLA-6	1098	0	1098	+	contig_4	80000	100	50000	900	1100	60	tp:A:P	dv:f:0.01
"""


def test_probe_assignment_requires_identity_and_coverage():
    assert assign_probe(identity=0.99, coverage=0.95, mapq=60) == HIGH_CONFIDENCE


def test_low_identity_probe_is_not_a_locus_call():
    assert assign_probe(identity=0.85, coverage=0.95, mapq=60) == UNRESOLVED


def test_low_coverage_or_mapq_is_not_a_locus_call():
    assert assign_probe(identity=0.99, coverage=0.50, mapq=60) == UNRESOLVED
    assert assign_probe(identity=0.99, coverage=0.95, mapq=10) == UNRESOLVED


def test_parse_paf_uses_dv_identity_and_query_coverage():
    hits = parse_paf(PAF_FIXTURE)
    first = hits[0]
    assert first.allele == "BoLA-2*012:01"
    assert first.locus == "BoLA-2"
    assert first.identity == 0.99
    assert first.coverage == 1.0
    assert first.mapq == 60
    assert first.is_primary is True


def test_low_identity_and_overlong_span_are_excluded_from_clusters():
    copies = cluster_alignments(parse_paf(PAF_FIXTURE))
    loci = {copy.best.locus for copy in copies}
    assert "BoLA-1" not in loci
    assert "BoLA-6" not in loci
    assert loci == {"BoLA-2", "BoLA-3"}


def test_nearby_hits_on_one_contig_collapse_to_one_copy():
    copies = cluster_alignments(parse_paf(PAF_FIXTURE))
    bola2 = [copy for copy in copies if copy.best.locus == "BoLA-2"]
    assert len(bola2) == 1
    assert bola2[0].best.allele == "BoLA-2*012:01"
    assert bola2[0].n_distinct_alleles == 2


def test_two_panel_alleles_in_one_cluster_are_not_collapse():
    """Many IPD alleles hit one gene copy; that is panel competition, not haplotype collapse."""
    assert classify_locus_structure(n_contigs=1, n_clusters=1, n_distinct_alleles_in_largest_cluster=12) == RECOVERED


def test_overlapping_paralog_hits_keep_the_best_locus_only():
    paf = (
        "BoLA-2*026:03|BoLA-2\t1098\t0\t1098\t+\tcontig_7\t50000\t28507\t31786\t1090\t1100\t60\ttp:A:P\tdv:f:0.001\n"
        "BoLA-1*049:01|BoLA-1\t1098\t0\t1098\t+\tcontig_7\t50000\t28509\t31786\t1080\t1100\t60\ttp:A:P\tdv:f:0.01\n"
    )
    copies = cluster_alignments(parse_paf(paf))
    assert len(copies) == 1
    assert copies[0].best.locus == "BoLA-2"


def test_loci_farther_than_gap_on_one_contig_stay_separate():
    paf = (
        "BoLA-2*012:01|BoLA-2\t1098\t0\t1098\t+\tcontig_1\t80000\t1000\t4000\t1080\t1100\t60\ttp:A:P\tdv:f:0.01\n"
        "BoLA-3*004:03|BoLA-3\t1098\t0\t1098\t+\tcontig_1\t80000\t25000\t28000\t1080\t1100\t60\ttp:A:P\tdv:f:0.01\n"
    )
    copies = cluster_alignments(parse_paf(paf))
    assert sorted(copy.best.locus for copy in copies) == ["BoLA-2", "BoLA-3"]


def test_same_locus_on_two_contigs_is_haplotype_split_not_copy_number():
    assert classify_locus_structure(n_contigs=2, n_clusters=2, n_distinct_alleles_in_largest_cluster=1) == HAPLOTYPE_SPLIT


def test_two_clusters_on_one_contig_are_duplicated_not_diploid_cn():
    assert classify_locus_structure(n_contigs=1, n_clusters=2, n_distinct_alleles_in_largest_cluster=1) == DUPLICATED_ON_CONTIG


def test_missing_locus_is_not_recovered():
    assert classify_locus_structure(n_contigs=0, n_clusters=0, n_distinct_alleles_in_largest_cluster=0) == NOT_RECOVERED


def test_single_clean_cluster_is_recovered():
    assert classify_locus_structure(n_contigs=1, n_clusters=1, n_distinct_alleles_in_largest_cluster=1) == RECOVERED


def test_extract_command_uses_target_bam_region_and_drops_secondary():
    argv = build_extract_view_argv("s.targets.bam", CLASS_I_REGION)
    assert argv[:4] == ["samtools", "view", "-b", "-F"]
    assert argv[4] == "0x900"
    assert argv[-1] == CLASS_I_REGION
    assert "Deliver_ADD" not in " ".join(argv)


def test_flye_command_keeps_haplotypes_and_pins_class_i_genome_size():
    argv = build_flye_argv("reads.fq", "out/flye")
    assert "--nano-raw" in argv
    assert "--keep-haplotypes" in argv
    assert argv[argv.index("--min-overlap") + 1] == "3000"
    assert argv[argv.index("--genome-size") + 1] == CLASS_I_GENOME_SIZE


def test_ipd_map_command_is_splice_aware_with_secondary_kept():
    argv = build_ipd_map_argv("assembly.fasta", "panel.fasta")
    assert argv[argv.index("-x") + 1] == "splice"
    assert "--secondary=yes" in argv
    assert argv[-2:] == ["assembly.fasta", "panel.fasta"]


def test_recovery_table_requires_the_full_cohort():
    rows = [{"sample_id": "ID01", "status": "passed"}]
    with pytest.raises(CohortError, match="missing"):
        validate_cohort_ids(["ID01", "ID02"], [row["sample_id"] for row in rows])


def test_report_states_cluster_counts_are_not_validated_copy_numbers():
    copies = cluster_alignments(parse_paf(PAF_FIXTURE))
    summary = summarize_sample(
        sample_id="ID02",
        copies=copies,
        n_reads=478,
        n_bases=7_000_000,
        contig_lengths=[50000, 40000],
        status="passed",
        hap1="*015:01",
        hap2="*016:01",
    )
    report = render_report([summary], interval=CLASS_I_REGION)
    assert COPY_NUMBER_CAVEAT in report
    assert "not a diploid copy number" in report.lower()
    assert "4-digit" in report.lower()
    assert summary["n_classical_clusters"] == 2
    assert "copy_number" not in summary
