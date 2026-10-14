from __future__ import annotations

from scripts.codex_reanalysis.call_small_variants import (
    build_bcftools_mpileup_command,
    build_samtools_mpileup_command,
    classify_dosage,
    genotype_label,
    parse_bcftools_ad_vcf,
    parse_samtools_pileup_line,
)


def test_build_bcftools_mpileup_command_requests_ad_and_single_position():
    cmd = build_bcftools_mpileup_command(
        "s.bam", "ref.fa", "NC_1", 100, bcftools="bcftools", mapq=20, baseq=20
    )
    assert cmd == [
        "bcftools",
        "mpileup",
        "-f",
        "ref.fa",
        "-r",
        "NC_1:100-100",
        "-a",
        "AD,DP",
        "--no-BAQ",
        "-q",
        "20",
        "-Q",
        "20",
        "-Ov",
        "s.bam",
    ]


def test_build_samtools_mpileup_command_requests_single_position():
    cmd = build_samtools_mpileup_command("s.bam", "ref.fa", "NC_1", 100, samtools="samtools", mapq=20, baseq=20)
    assert cmd == ["samtools", "mpileup", "-f", "ref.fa", "-r", "NC_1:100-100", "-q", "20", "-Q", "20", "s.bam"]


def test_parse_bcftools_ad_vcf_reads_ref_and_alt_counts():
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "NC_037350.1\t27719770\t.\tT\tC,<*>\t0\t.\tDP=27\tPL:DP:AD\t45,0,180,81,192,241:16:12,4,0\n"
    )
    result = parse_bcftools_ad_vcf(text)
    assert result == {"ref_base": "T", "alt_base": "C", "ref_count": 12, "alt_count": 4, "depth": 16}


def test_parse_bcftools_ad_vcf_handles_no_alt_observed():
    text = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        "NC_037350.1\t100\t.\tA\t<*>\t0\t.\tDP=10\tPL:DP:AD\t0,30,300:10:10,0\n"
    )
    result = parse_bcftools_ad_vcf(text)
    assert result == {"ref_base": "A", "alt_base": None, "ref_count": 10, "alt_count": 0, "depth": 10}


def test_parse_bcftools_ad_vcf_handles_empty_output():
    assert parse_bcftools_ad_vcf("") == {
        "ref_base": None,
        "alt_base": None,
        "ref_count": 0,
        "alt_count": 0,
        "depth": 0,
    }


def test_parse_samtools_pileup_line_counts_substitutions_only():
    # "*" (an overlapping deletion at this exact position) is excluded from
    # both ref and alt counts, so 11 pileup characters yield depth 10.
    text = "NC_037350.1\t27719770\tT\t11\t*.C,..,CC..\t:GC96=:8<=>\n"
    result = parse_samtools_pileup_line(text)
    assert result == {"ref_base": "T", "alt_base": "C", "ref_count": 7, "alt_count": 3, "depth": 10}


def test_parse_samtools_pileup_line_skips_indel_annotations():
    # "." then a 2bp insertion "+2AT" attached to it, then two more ref calls.
    text = "chr1\t5\tA\t3\t.+2AT,,\tIII\n"
    result = parse_samtools_pileup_line(text)
    assert result == {"ref_base": "A", "alt_base": None, "ref_count": 3, "alt_count": 0, "depth": 3}


def test_parse_samtools_pileup_line_skips_read_start_and_end_markers():
    text = "chr1\t5\tA\t2\t^].,$\tII\n"
    result = parse_samtools_pileup_line(text)
    assert result == {"ref_base": "A", "alt_base": None, "ref_count": 2, "alt_count": 0, "depth": 2}


def test_classify_dosage_ref_below_15_percent_alt():
    assert classify_dosage(ref_count=90, alt_count=10) == "ref"


def test_classify_dosage_het_between_15_and_85_percent():
    assert classify_dosage(ref_count=60, alt_count=40) == "het"


def test_classify_dosage_alt_above_85_percent():
    assert classify_dosage(ref_count=5, alt_count=95) == "alt"


def test_classify_dosage_no_data_when_no_reads():
    assert classify_dosage(ref_count=0, alt_count=0) == "no_data"


def test_genotype_label_homozygous_ref():
    assert genotype_label("T", "C", "ref") == "T/T"


def test_genotype_label_heterozygous_is_allele_sorted():
    assert genotype_label("T", "C", "het") == "C/T"


def test_genotype_label_homozygous_alt():
    assert genotype_label("T", "C", "alt") == "C/C"


def test_genotype_label_no_data():
    assert genotype_label(None, None, "no_data") == "no_data"
