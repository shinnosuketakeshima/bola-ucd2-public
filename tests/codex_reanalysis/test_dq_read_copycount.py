from __future__ import annotations

import pytest

from scripts.codex_reanalysis.dq_read_copycount import (
    Hit,
    count_family_copies,
    merge_intervals,
    parse_paf_hits,
    parse_panel_family,
    summarise_sample,
)


def test_parse_panel_family_accepts_bola_classii_only():
    assert parse_panel_family(">BoLA-DQA*001:01 BoLA02987") == "DQA"
    assert parse_panel_family(">BoLA-DQB*010:01:01 BoLA03100") == "DQB"
    assert parse_panel_family(">BoLA-DRB3*016:01 BoLA03200") == "DRB3"
    assert parse_panel_family(">BoLA-2*012:01 BoLA02982") is None
    assert parse_panel_family(">Bubu-DQA*001:01 x") is None


def test_merge_intervals_gene_span_collapses_exons():
    # Two exons 2 kb apart merge under a 15 kb gene span.
    merged = merge_intervals([(1000, 1500), (3000, 3500), (50000, 51000)], gap=15_000)
    assert merged == [(1000, 3500), (50000, 51000)]


def test_count_family_copies_separates_families_and_distant_paralogs():
    hits = [
        Hit("r1", "DQA", "BoLA-DQA*001:01", 0, 500, 0.9, 500),
        Hit("r1", "DQA", "BoLA-DQA*002:01", 2000, 2500, 0.85, 500),  # same gene (exons)
        Hit("r1", "DQA", "BoLA-DQA*010:01", 60_000, 60_500, 0.88, 500),  # second copy
        Hit("r1", "DQB", "BoLA-DQB*001:01", 10_000, 10_400, 0.9, 400),
    ]
    copies = count_family_copies(hits, gene_span=15_000)
    assert copies == {"DQA": 2, "DQB": 1}


def test_parse_paf_hits_filters_by_identity_and_family():
    paf = "\n".join(
        [
            # read, qlen, qs, qe, strand, target, tlen, ts, te, nmatch, aln, mapq, tags
            "read1\t20000\t100\t600\t+\tBoLA-DQA*001:01\t800\t0\t500\t480\t500\t60\tdv:f:0.05",
            "read1\t20000\t200\t700\t+\tBoLA-2*012:01\t1000\t0\t500\t480\t500\t60\tdv:f:0.05",  # class I ignored
            "read1\t20000\t8000\t8200\t+\tBoLA-DQB*001:01\t700\t0\t200\t100\t200\t60\tdv:f:0.40",  # low id
        ]
    )
    by_read = parse_paf_hits(paf, min_aln=150, min_identity=0.70)
    assert list(by_read) == ["read1"]
    assert len(by_read["read1"]) == 1
    assert by_read["read1"][0].family == "DQA"


def test_summarise_sample_uses_long_reads_as_pct_denominator():
    by_read = {
        "r1": [Hit("r1", "DQA", "a", 0, 100, 0.9, 100), Hit("r1", "DQA", "b", 50_000, 50_100, 0.9, 100)],
        "r2": [Hit("r2", "DQB", "c", 0, 100, 0.9, 100)],
    }
    # 100 long reads analysed; only 1 has DQA≥2 → 1.0%
    summary = summarise_sample(by_read, n_long_reads=100, gene_span=15_000)
    assert summary["long_reads"] == 100
    assert summary["reads_with_dqa_or_dqb_hit"] == 2
    assert summary["pct_DQA_ge2"] == 1.0
    assert summary["max_DQA_cis"] == 2
    assert summary["pct_DQB_ge2"] == 0.0
    assert summary["max_DQB_cis"] == 1
