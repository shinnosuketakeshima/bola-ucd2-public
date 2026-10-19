from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.codex_reanalysis.call_coding_variants import (
    BREED_DIVERGENT,
    CDS_EDGE,
    CodingVariantError,
    INDEL,
    MISSENSE,
    NONCODING_IN_SPAN,
    RARE,
    SEGREGATING,
    STOPGAIN,
    SYNONYMOUS,
    build_gene_sequence,
    build_mpileup_call_command,
    build_variant_table,
    classify_allele_frequency,
    classify_effect,
    complement,
    gene_for_position,
    load_cds_models_from_bed,
    load_cds_models_from_json,
    load_cntn3_cds_models,
    parse_vcf_calls,
    render_report,
    verify_contig_identical,
    write_variant_table_tsv,
)


def test_complement_is_case_preserving():
    assert complement("ACGT") == "TGCA"
    assert complement("acgt") == "tgca"


def test_classify_allele_frequency_thresholds():
    assert classify_allele_frequency(0.90) == BREED_DIVERGENT
    assert classify_allele_frequency(0.85) == BREED_DIVERGENT
    assert classify_allele_frequency(0.50) == SEGREGATING
    assert classify_allele_frequency(0.15) == SEGREGATING
    assert classify_allele_frequency(0.10) == RARE


def test_build_gene_sequence_plus_strand_single_exon():
    # region "ATGGCC..." starting at genomic pos 100; single exon 100-105 (ATGGCC).
    region_seq = "ATGGCC"
    coding_seq, gene_index = build_gene_sequence(region_seq, 100, [(100, 105)], "+")
    assert coding_seq == "ATGGCC"
    assert gene_index == {100: 0, 101: 1, 102: 2, 103: 3, 104: 4, 105: 5}


def test_build_gene_sequence_minus_strand_reverses_and_complements():
    # region "GGCCAT" at genomic pos 100-105; minus strand -> reverse-complement "ATGGCC".
    region_seq = "GGCCAT"
    coding_seq, gene_index = build_gene_sequence(region_seq, 100, [(100, 105)], "-")
    assert coding_seq == "ATGGCC"
    # highest genomic position maps to coding index 0 on minus strand.
    assert gene_index[105] == 0
    assert gene_index[100] == 5


def test_classify_effect_synonymous_third_position_wobble():
    # ATG GCC -> codon1=ATG(M), codon2=GCC(A). GCC->GCT is synonymous (both Ala).
    gene_seq = "ATGGCC"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 105)], "+")
    effect, detail = classify_effect(gene_seq, gene_index, 105, "C", "T", "+")
    assert effect == SYNONYMOUS
    assert "GCC" in detail and "GCT" in detail


def test_classify_effect_missense():
    # ATG GCC, position 104 (2nd base of codon2, G->A): GCC->GAC = Ala->Asp, missense.
    gene_seq = "ATGGCC"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 105)], "+")
    effect, detail = classify_effect(gene_seq, gene_index, 104, "C", "A", "+")
    assert effect == MISSENSE


def test_classify_effect_stopgain():
    # TGG (Trp) -> TAG (stop) at 2nd position.
    gene_seq = "TGG"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 102)], "+")
    effect, _ = classify_effect(gene_seq, gene_index, 101, "G", "A", "+")
    assert effect == STOPGAIN


def test_classify_effect_indel_never_scored():
    gene_seq = "ATGGCC"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 105)], "+")
    assert classify_effect(gene_seq, gene_index, 100, "AT", "A", "+")[0] == INDEL


def test_classify_effect_position_outside_cds_is_noncoding():
    gene_seq = "ATGGCC"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 105)], "+")
    assert classify_effect(gene_seq, gene_index, 999, "A", "T", "+")[0] == NONCODING_IN_SPAN


def test_classify_effect_partial_trailing_codon_is_cds_edge():
    # 4-base CDS: codon0 complete (positions 0-2), position 3 starts codon1 but only 1 base exists.
    gene_seq = "ATGG"
    _, gene_index = build_gene_sequence(gene_seq, 100, [(100, 103)], "+")
    effect, _ = classify_effect(gene_seq, gene_index, 103, "G", "A", "+")
    assert effect == CDS_EDGE


def test_parse_vcf_calls_reads_first_alt_and_normalizes_phased_gt():
    vcf_body = "\n".join(
        [
            "#header line ignored",
            "\t".join(["chr1", "100", ".", "A", "T,C", ".", "PASS", "INFO", "GT:AD", "0|1:5,3"]),
        ]
    )
    calls = parse_vcf_calls(vcf_body)
    assert calls == {100: ("A", "T", "0/1")}


def test_parse_vcf_calls_rejects_short_lines():
    with pytest.raises(CodingVariantError):
        parse_vcf_calls("chr1\t100\t.\tA\tT\n")


def test_gene_for_position_finds_containing_exon():
    models = {"GENE1": {"contig": "c", "strand": "+", "exons": [(100, 110), (200, 210)]}}
    assert gene_for_position(models, 105) == "GENE1"
    assert gene_for_position(models, 205) == "GENE1"
    assert gene_for_position(models, 150) is None


def test_build_variant_table_computes_allele_frequency_and_effect():
    # Per-sample bcftools-call VCFs only carry variant sites, so a sample absent from a
    # position's calls (s3 here) is not "called ref" -- it simply contributes no record,
    # matching scripts/btn_coding_aggregate.py's original n_called semantics exactly.
    models = {"GENE1": {"contig": "c", "strand": "+", "exons": [(100, 105)]}}
    sequences_by_gene = {"GENE1": build_gene_sequence("ATGGCC", 100, [(100, 105)], "+")}
    calls_by_sample = {
        "s1": {105: ("C", "T", "1/1")},
        "s2": {105: ("C", "T", "0/1")},
        "s3": {},
    }
    rows = build_variant_table(models, sequences_by_gene, calls_by_sample)
    assert len(rows) == 1
    row = rows[0]
    assert row["gene"] == "GENE1"
    assert row["pos"] == 105
    assert row["n_carriers"] == 2
    assert row["n_hom_alt"] == 1
    assert row["n_called"] == 2
    # (2*1 hom + 1 het) / (2*2 called) = 3/4 = 0.75 -> segregating
    assert row["alt_af"] == 0.75
    assert row["allele_frequency_class"] == SEGREGATING
    assert row["effect"] == SYNONYMOUS


def test_build_variant_table_ignores_positions_outside_any_gene():
    models = {"GENE1": {"contig": "c", "strand": "+", "exons": [(100, 105)]}}
    sequences_by_gene = {"GENE1": build_gene_sequence("ATGGCC", 100, [(100, 105)], "+")}
    calls_by_sample = {"s1": {9999: ("A", "T", "1/1")}}
    rows = build_variant_table(models, sequences_by_gene, calls_by_sample)
    assert rows == []


def test_load_cds_models_from_json_requires_declared_strand(tmp_path):
    payload = {"GENE1": {"chrom": "c", "cds": [[100, 105]]}}
    path = tmp_path / "models.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(CodingVariantError, match="strand"):
        load_cds_models_from_json(path, strand_by_gene={})


def test_load_cds_models_from_json_parses_exons_and_strand(tmp_path):
    payload = {"GENE1": {"chrom": "NC_1", "cds": [[100, 105], [200, 210]]}}
    path = tmp_path / "models.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    models = load_cds_models_from_json(path, strand_by_gene={"GENE1": "-"})
    assert models["GENE1"]["contig"] == "NC_1"
    assert models["GENE1"]["strand"] == "-"
    assert models["GENE1"]["exons"] == [(100, 105), (200, 210)]


def test_load_cds_models_from_bed_converts_zero_based_half_open(tmp_path):
    bed = tmp_path / "cds.bed"
    bed.write_text("NC_1\t99\t105\tGENE1\t+\n", encoding="utf-8")
    models = load_cds_models_from_bed(bed)
    assert models["GENE1"]["exons"] == [(100, 105)]
    assert models["GENE1"]["strand"] == "+"


def test_load_cds_models_from_bed_uses_strand_fallback_when_column_absent(tmp_path):
    bed = tmp_path / "cds.bed"
    bed.write_text("NC_1\t99\t105\tGENE1\n", encoding="utf-8")
    models = load_cds_models_from_bed(bed, strand_by_gene={"GENE1": "+"})
    assert models["GENE1"]["strand"] == "+"


def test_load_cds_models_from_bed_raises_without_strand(tmp_path):
    bed = tmp_path / "cds.bed"
    bed.write_text("NC_1\t99\t105\tGENE1\n", encoding="utf-8")
    with pytest.raises(CodingVariantError, match="strand"):
        load_cds_models_from_bed(bed)


def test_load_cntn3_cds_models_renames_bed_gene_column_to_cntn3(tmp_path):
    """data/meta/cntn3_cds.bed's own gene column reads "CNTN3_CDS", not "CNTN3"."""
    bed = tmp_path / "cntn3_cds.bed"
    bed.write_text("NC_037349.1\t99\t105\tCNTN3_CDS\n", encoding="utf-8")
    models = load_cntn3_cds_models(bed)
    assert "CNTN3_CDS" not in models
    assert models["CNTN3"]["exons"] == [(100, 105)]
    assert models["CNTN3"]["strand"] == "+"


def test_load_cds_models_from_bed_rejects_gene_spanning_multiple_contigs(tmp_path):
    bed = tmp_path / "cds.bed"
    bed.write_text("NC_1\t99\t105\tGENE1\t+\nNC_2\t199\t205\tGENE1\t+\n", encoding="utf-8")
    with pytest.raises(CodingVariantError, match="multiple contigs"):
        load_cds_models_from_bed(bed)


def test_verify_contig_identical_true_for_identical_status(tmp_path):
    tsv = tmp_path / "sequence_comparison.tsv"
    tsv.write_text("new_accession\tstatus\nNC_037350.1\tidentical\n", encoding="utf-8")
    assert verify_contig_identical(tsv, "NC_037350.1") is True


def test_verify_contig_identical_false_for_non_identical_status(tmp_path):
    tsv = tmp_path / "sequence_comparison.tsv"
    tsv.write_text("new_accession\tstatus\nNC_037350.1\tdecision_required\n", encoding="utf-8")
    assert verify_contig_identical(tsv, "NC_037350.1") is False


def test_verify_contig_identical_false_when_contig_absent(tmp_path):
    tsv = tmp_path / "sequence_comparison.tsv"
    tsv.write_text("new_accession\tstatus\nNC_999.1\tidentical\n", encoding="utf-8")
    assert verify_contig_identical(tsv, "NC_037350.1") is False


def test_write_variant_table_tsv_round_trips(tmp_path):
    rows = [
        {
            "gene": "GENE1",
            "pos": 105,
            "ref": "C",
            "alt": "T",
            "allele_frequency_class": SEGREGATING,
            "alt_af": 0.5,
            "n_carriers": 2,
            "n_hom_alt": 1,
            "n_called": 3,
            "sample_count": 3,
            "effect": SYNONYMOUS,
            "detail": "GCC(A)>GCT(A) codon2",
        }
    ]
    output = tmp_path / "coding_variants.tsv"
    write_variant_table_tsv(output, rows)
    text = output.read_text(encoding="utf-8")
    assert "GENE1" in text
    assert "segregating" in text


def test_render_report_summarizes_per_gene():
    rows = [
        {"gene": "GENE1", "allele_frequency_class": SEGREGATING, "effect": MISSENSE},
        {"gene": "GENE1", "allele_frequency_class": BREED_DIVERGENT, "effect": SYNONYMOUS},
    ]
    report = render_report(rows, ["GENE1"])
    assert "GENE1" in report
    assert "2 coding variants" in report


def test_build_mpileup_call_command_pins_ont_baq_and_thresholds():
    command = build_mpileup_call_command("bcftools", Path("ref.fa"), Path("cds.bed"), Path("s.bam"), Path("out.vcf.gz"))
    assert "-B -d 60 -q 20 -Q 10" in command
    assert "call -mv" in command
