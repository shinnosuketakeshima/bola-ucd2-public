from __future__ import annotations

from collections import Counter

from scripts.codex_reanalysis.drb3_wgs_typing import (
    AlleleHit,
    format_top_alleles,
    known_presence_label,
    parse_drb3_header,
    parse_paf_best_per_read,
    summarise_read_typing,
    write_drb3_panel,
)


def test_parse_drb3_header_accepts_bola_drb3_only():
    assert parse_drb3_header(">BoLA-DRB3*016:01 BoLA03200") == "*016:01"
    assert parse_drb3_header(">BoLA-DRB3*010:01:01 x") == "*010:01:01"
    assert parse_drb3_header(">BoLA-DQA*001:01 x") is None
    assert parse_drb3_header(">BoLA-2*012:01 x") is None


def test_write_drb3_panel_filters(tmp_path):
    full = tmp_path / "panel.fa"
    full.write_text(
        ">BoLA-DRB3*016:01 a\nACGT\n>BoLA-DQA*001:01 b\nTTTT\n>BoLA-DRB3*010:01 c\nGGGG\n",
        encoding="utf-8",
    )
    out = tmp_path / "drb3.fa"
    n = write_drb3_panel(full, out)
    text = out.read_text(encoding="utf-8")
    assert n == 2
    assert "BoLA-DRB3*016:01" in text
    assert "BoLA-DRB3*010:01" in text
    assert "BoLA-DQA" not in text


def test_parse_paf_best_per_read_picks_highest_identity():
    paf = "\n".join(
        [
            "r1\t10000\t0\t500\t+\tBoLA-DRB3*016:01\t800\t0\t500\t480\t500\t60\tdv:f:0.05",
            "r1\t10000\t0\t500\t+\tBoLA-DRB3*010:01\t800\t0\t500\t400\t500\t60\tdv:f:0.20",
            "r1\t10000\t0\t100\t+\tBoLA-DRB3*027:03\t800\t0\t100\t90\t100\t60\tdv:f:0.01",  # short aln
            "r2\t10000\t0\t500\t+\tBoLA-DQA*001:01\t800\t0\t500\t480\t500\t60\tdv:f:0.01",  # wrong family
        ]
    )
    best = parse_paf_best_per_read(paf, min_aln=200, min_identity=0.70)
    assert list(best) == ["r1"]
    assert best["r1"].allele == "*016:01"
    assert best["r1"].identity == 0.95


def test_summarise_flags_central_bias_for_homozygous_010():
    best = {
        f"r{i}": AlleleHit(f"r{i}", "*016:01", 0.95, 500, 500) for i in range(10)
    }
    best["r10"] = AlleleHit("r10", "*010:01", 0.90, 400, 400)
    summary = summarise_read_typing(best, hap1="*010:01", hap2="*010:01", n_reads=20)
    assert summary["top1_allele"] == "*016:01"
    assert summary["top1_is_central_bias_artifact"] is True
    assert summary["top1_matches_expected"] is False
    assert "010" in str(summary["known_alleles_present"])


def test_known_presence_and_format_helpers():
    assert format_top_alleles(Counter({"*016:01": 31, "*010:01": 2})) == "*016:01(31),*010:01(2)"
    assert known_presence_label({"*010:01"}, {"*010:01", "*016:01"}, top_alleles=["*016:01", "*010:01"]) == (
        "*010:01 present"
    )
    assert known_presence_label({"*005:02", "*005:03"}, {"*005:02"}, top_alleles=["*016:01", "*005:02"]) == (
        "*005:02 only(top3)"
    )
