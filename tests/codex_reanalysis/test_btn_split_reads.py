import pytest

from scripts.codex_reanalysis.btn_split_reads import parse_region, partner_hits, summarize_records

PARTNERS = [("NC_037350.1", 26240001, 26270000), ("NC_037350.1", 25940001, 25950000)]


def record(flag, mapq, sa=""):
    tags = f"\tSA:Z:{sa}" if sa else ""
    return f"r\t{flag}\tNC_037350.1\t25870000\t{mapq}\t100M\t*\t0\t0\tA\t*{tags}"


def test_partner_needs_window_and_mapq():
    tag = "NC_037350.1,26250000,+,50S50M,60,0;NC_037350.1,25945000,-,50S50M,5,0;"
    assert partner_hits(tag, PARTNERS, 20) == {0}


def test_summary_counts_each_read_once_per_partner():
    lines = [
        record(0, 60, "NC_037350.1,26250000,+,5S5M,60,0;NC_037350.1,26260000,-,5S5M,60,0;"),
        record(0, 60),
        record(0, 3, "NC_037350.1,26250000,+,5S5M,60,0;"),  # low MAPQ primary is not a partner carrier
        record(256, 0),
        record(2048, 60, "NC_037350.1,25870000,+,5S5M,60,0;"),
    ]
    summary = summarize_records(lines, PARTNERS, 20)
    assert summary["primary_reads"] == 3
    assert summary["low_mapq_fraction"] == "0.33"
    assert summary["secondary_alignments"] == 1
    assert summary["primary_mapq_reads"] == 2
    assert summary["partner_counts"] == [1, 0]


def test_region_must_be_well_formed():
    assert parse_region("NC_037350.1:1-10") == ("NC_037350.1", 1, 10)
    with pytest.raises(ValueError):
        parse_region("NC_037350.1:10-1")
