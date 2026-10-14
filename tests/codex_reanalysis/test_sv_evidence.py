from __future__ import annotations

from pathlib import Path

import pytest

from scripts.codex_reanalysis.merge_sv_evidence import (
    EXPLORATORY,
    RECURRENT_REFERENCE_DIFFERENCE,
    SUPPORTED,
    SvParseError,
    build_call_status,
    build_sv_table,
    classify_sv,
    group_by_breakpoint,
    is_recurrent_reference_difference,
    parse_sniffles_vcf,
    render_report,
)

# Header and a handful of data lines captured verbatim from a real Sniffles2
# 2.8.0 run (`sniffles --input ID02.targets.bam`) against this
# project's target-only BAM, so the parser is checked against the actual
# field names/format this tool emits rather than the spec alone.
SNIFFLES_FIXTURE = """\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
##INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of reads supporting a DEL SV through supplementary alignments (split events)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
NC_037336.1\t104055383\tSniffles2.DEL.2S8\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;SVLEN=-201;END=104055584;SUPPORT=14;SUPPORT_SA=9;STRAND=+-;VAF=0.560\tGT:GQ:DR:DV:PS\t0/1:60:11:14:.
NC_037340.1\t47147887\tSniffles2.INS.1SC\tN\tTACTAGCTGCATCACTTATATCC\t60\tPASS\tPRECISE;SVTYPE=INS;SVLEN=102;END=47147887;SUPPORT=4;VAF=0.800\tGT:GQ:DR:DV:PS\t1/1:1:1:4:.
NC_037350.1\t26707388\tSniffles2.BND.3E0S16\tN\t[NC_037350.1:26813151[N\t60\tPASS\tPRECISE;SVTYPE=BND;SUPPORT=14;CHR2=NC_037350.1;VAF=0.438\tGT:GQ:DR:DV:PS\t0/1:60:18:14:.
NC_037350.1\t25686861\tSniffles2.INV.39ES16\tN\t<INV>\t60\tPASS\tPRECISE;SVTYPE=INV;SVLEN=519;END=25687380;SUPPORT=20;STRAND=-;VAF=0.417\tGT:GQ:DR:DV:PS\t0/1:60:28:20:.
NC_037350.1\t71159075\tSniffles2.INS.8S13\tN\tTGCTACAGTTTGTGGG\t60\tPASS\tIMPRECISE;SVTYPE=INS;SVLEN=441;END=71159075;SUPPORT=23;VAF=0.885\tGT:GQ:DR:DV:PS\t1/1:34:3:23:.
"""


def test_caller_only_event_is_exploratory():
    assert classify_sv(callers=1, split_read_support=0, assembly_support=False) == EXPLORATORY


def test_independent_breakpoint_support_is_accepted_with_scope():
    assert classify_sv(callers=2, split_read_support=12, assembly_support=True) == SUPPORTED


def test_second_caller_alone_is_supported():
    assert classify_sv(callers=2, split_read_support=0, assembly_support=False) == SUPPORTED


def test_split_read_support_alone_is_supported():
    assert classify_sv(callers=1, split_read_support=1, assembly_support=False) == SUPPORTED


def test_assembly_support_alone_is_supported():
    assert classify_sv(callers=1, split_read_support=0, assembly_support=True) == SUPPORTED


def test_split_read_support_below_minimum_is_exploratory():
    assert (
        classify_sv(callers=1, split_read_support=2, assembly_support=False, minimum_split_read_support=5)
        == EXPLORATORY
    )


def test_recurrent_reference_difference_majority_threshold():
    assert is_recurrent_reference_difference(sample_count=8, cohort_size=13) is True
    assert is_recurrent_reference_difference(sample_count=6, cohort_size=13) is False
    assert is_recurrent_reference_difference(sample_count=7, cohort_size=13) is True


def test_recurrent_reference_difference_rejects_nonpositive_cohort():
    with pytest.raises(ValueError, match="cohort_size"):
        is_recurrent_reference_difference(sample_count=1, cohort_size=0)


def test_parse_sniffles_vcf_reads_del_with_split_read_support():
    records = parse_sniffles_vcf(SNIFFLES_FIXTURE)
    deletion = next(record for record in records if record["sv_type"] == "DEL")
    assert deletion["chrom"] == "NC_037336.1"
    assert deletion["pos"] == 104055383
    assert deletion["end"] == 104055584
    assert deletion["sv_len"] == -201
    assert deletion["support"] == 14
    assert deletion["support_sa"] == 9
    assert deletion["precise"] is True
    assert deletion["genotype"] == "0/1"
    assert deletion["ref_reads"] == 11
    assert deletion["alt_reads"] == 14


def test_parse_sniffles_vcf_reads_bnd_mate_from_alt_bracket_notation():
    records = parse_sniffles_vcf(SNIFFLES_FIXTURE)
    bnd = next(record for record in records if record["sv_type"] == "BND")
    assert bnd["chr2"] == "NC_037350.1"
    assert bnd["mate_pos"] == 26813151


def test_parse_sniffles_vcf_defaults_missing_support_sa_to_zero():
    records = parse_sniffles_vcf(SNIFFLES_FIXTURE)
    insertion = next(record for record in records if record["id"] == "Sniffles2.INS.1SC")
    assert insertion["support_sa"] == 0


def test_parse_sniffles_vcf_reads_imprecise_flag():
    records = parse_sniffles_vcf(SNIFFLES_FIXTURE)
    imprecise = next(record for record in records if record["id"] == "Sniffles2.INS.8S13")
    assert imprecise["precise"] is False


def test_parse_sniffles_vcf_rejects_short_data_line():
    with pytest.raises(SvParseError):
        parse_sniffles_vcf("#CHROM\tPOS\nchr1\t1\n")


def test_group_by_breakpoint_clusters_nearby_positions_across_samples():
    calls = [
        {"sample_id": "a", "chrom": "chr1", "sv_type": "DEL", "pos": 1000, "support": 5, "support_sa": 2, "precise": True},
        {"sample_id": "b", "chrom": "chr1", "sv_type": "DEL", "pos": 1010, "support": 6, "support_sa": 3, "precise": True},
        {"sample_id": "c", "chrom": "chr1", "sv_type": "DEL", "pos": 5000, "support": 4, "support_sa": 0, "precise": False},
    ]
    groups = group_by_breakpoint(calls, position_tolerance=50)
    assert len(groups) == 2
    clustered = next(group for group in groups if group["sample_count"] == 2)
    assert clustered["sample_ids"] == ["a", "b"]
    assert clustered["total_support_sa"] == 5
    assert clustered["any_imprecise"] is False
    distant = next(group for group in groups if group["sample_count"] == 1)
    assert distant["any_imprecise"] is True


def test_build_sv_table_flags_recurrent_events_before_classification():
    groups = [
        {
            "chrom": "chr1",
            "sv_type": "DEL",
            "position_min": 1000,
            "position_max": 1010,
            "sample_count": 10,
            "total_support_sa": 20,
        }
    ]
    rows = build_sv_table(groups, cohort_size=13)
    assert rows[0]["status"] == RECURRENT_REFERENCE_DIFFERENCE


def test_build_sv_table_classifies_non_recurrent_single_caller_event_as_exploratory():
    groups = [
        {
            "chrom": "chr1",
            "sv_type": "DEL",
            "position_min": 1000,
            "position_max": 1010,
            "sample_count": 1,
            "total_support_sa": 0,
        }
    ]
    rows = build_sv_table(groups, cohort_size=13)
    assert rows[0]["status"] == EXPLORATORY
    assert rows[0]["callers"] == 1
    assert rows[0]["assembly_support"] is False


def test_sniffles_support_sa_does_not_self_promote_to_supported():
    groups = [
        {
            "chrom": "chr1",
            "sv_type": "DEL",
            "position_min": 1000,
            "position_max": 1010,
            "sample_count": 1,
            "total_support_sa": 20,
        }
    ]
    rows = build_sv_table(groups, cohort_size=13)
    assert rows[0]["status"] == EXPLORATORY
    assert rows[0]["sniffles_support_sa"] == 20
    assert rows[0]["split_read_support"] == 0


def test_build_sv_table_rejects_mismatched_evidence_list_length():
    groups = [
        {"chrom": "chr1", "sv_type": "DEL", "position_min": 1, "position_max": 1, "sample_count": 1, "total_support_sa": 0}
    ]
    with pytest.raises(ValueError, match="evidence"):
        build_sv_table(groups, cohort_size=13, callers_by_group=[1, 1])


def test_build_call_status_counts_sv_types_from_a_real_vcf(tmp_path: Path):
    vcf_path = tmp_path / "sample.sniffles2.vcf"
    vcf_path.write_text(SNIFFLES_FIXTURE, encoding="utf-8")
    status = build_call_status("sample", vcf_path, sniffles_version="2.8.0")
    assert status["sample_id"] == "sample"
    assert status["sv_count"] == 5
    assert status["sv_type_counts"] == {"BND": 1, "DEL": 1, "INS": 2, "INV": 1}
    assert status["status"] == "passed"


def test_build_call_status_fails_when_vcf_is_absent(tmp_path: Path):
    status = build_call_status("sample", tmp_path / "missing.vcf", sniffles_version="2.8.0")
    assert status["status"] == "failed"
    assert status["sv_count"] == 0


def test_render_report_only_lists_supported_events_in_the_supported_table():
    rows = [
        {
            "chrom": "chr1",
            "sv_type": "DEL",
            "position_min": 1000,
            "position_max": 1010,
            "sample_count": 2,
            "cohort_size": 13,
            "callers": 1,
            "split_read_support": 5,
            "assembly_support": False,
            "status": SUPPORTED,
        },
        {
            "chrom": "chr2",
            "sv_type": "INS",
            "position_min": 2000,
            "position_max": 2000,
            "sample_count": 1,
            "cohort_size": 13,
            "callers": 1,
            "split_read_support": 0,
            "assembly_support": False,
            "status": EXPLORATORY,
        },
    ]
    report = render_report(rows, sample_ids=[f"s{i}" for i in range(13)])
    assert "chr1" in report
    assert "chr2" not in report
    assert "Only `supported` events are eligible" in report
