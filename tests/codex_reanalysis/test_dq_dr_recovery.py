from __future__ import annotations

from scripts.codex_reanalysis.dq_dr_recovery import (
    AMBIGUOUS,
    DIPLOID_BASELINE,
    EXACT_TWO,
    EXPANDED,
    NO_RECOVERY,
    ONE_ALLELE,
    UNDER_ASSEMBLED,
    classify_dq_interpretation,
    classify_drb3_pcr_concordance,
    count_family_clusters,
    normalize_drb3_allele,
    parse_classii_query,
    parse_ipd_classii_paf,
)
from scripts.codex_reanalysis.probe_assignment import cluster_alignments


CLASSII_PAF = """\
BoLA-DQA*026:03 BoLA03028	1098	0	1098	+	contig_1	120000	1000	3200	1090	1100	60	tp:A:P	dv:f:0.01
BoLA-DQB*027:02 BoLA03087	1098	0	1098	+	contig_1	120000	52000	55200	1080	1100	60	tp:A:P	dv:f:0.005
BoLA-DRB3*016:01 BoLA03124	1098	0	1098	+	contig_2	80000	5000	8200	1090	1100	60	tp:A:P	dv:f:0.002
BoLA-DRB3*015:01 BoLA03122	1098	0	1090	+	contig_3	90000	5100	8300	1000	1100	60	tp:A:P	dv:f:0.03
"""


def test_parse_classii_query_extracts_family():
    allele, family, locus = parse_classii_query("BoLA-DQB*027:02 BoLA03087")
    assert allele == "BoLA-DQB*027:02"
    assert family == "DQB"
    assert locus == "DQB"


def test_normalize_drb3_allele_accepts_pcr_and_ipd_labels():
    assert normalize_drb3_allele("*016:01") == "*016:01"
    assert normalize_drb3_allele("BoLA-DRB3*016:01") == "*016:01"


def test_count_family_clusters_from_ipd_paf():
    copies = cluster_alignments(parse_ipd_classii_paf(CLASSII_PAF))
    counts = count_family_clusters(copies)
    assert counts["DQA"] == 1
    assert counts["DQB"] == 1
    assert counts["DRB3"] == 2


def test_diploid_baseline_requires_two_dqa_and_two_dqb_clusters():
    assert classify_dq_interpretation(2, 2, flye_status="passed") == DIPLOID_BASELINE
    assert classify_dq_interpretation(1, 1, flye_status="passed") == UNDER_ASSEMBLED
    assert classify_dq_interpretation(2, 4, flye_status="passed") == EXPANDED


def test_drb3_exact_two_allele_recovery():
    copies = [c for c in cluster_alignments(parse_ipd_classii_paf(CLASSII_PAF)) if c.best.locus == "DRB3"]
    assert classify_drb3_pcr_concordance("*015:01", "*016:01", copies) == EXACT_TWO


def test_missing_drb3_is_no_recovery():
    assert classify_drb3_pcr_concordance("*015:01", "*016:01", []) == NO_RECOVERY


def test_low_identity_drb3_is_not_recovered():
    paf = "BoLA-DRB3*016:01 BoLA03124\t1098\t0\t1098\t+\tcontig_1\t1000\t5000\t8200\t900\t1100\t60\ttp:A:P\tdv:f:0.20\n"
    copies = cluster_alignments(parse_ipd_classii_paf(paf))
    assert classify_drb3_pcr_concordance("*016:01", "*016:01", copies) == NO_RECOVERY
