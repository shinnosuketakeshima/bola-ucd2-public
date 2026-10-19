"""BTN + CNTN3/ABHD16A/PRR3 coding-variant screening (ARS-UCD2.0).

Reuses the `02_mapping` target BAMs directly -- BTN, CNTN3, ABHD16A, and PRR3 are all
already full-gene targets in `configs/codex_reanalysis.targets.tsv`, so no separate
mapping stage is needed for this analysis.
"""

CODING_ROOT = f"{MARKERS_ROOT}/coding_variants"


rule build_coding_cds_bed:
    input:
        btn_json="data/meta/btn_cds_models.json",
        # v2: rebuilt from a single canonical GFF transcript per gene (see
        # UCD2_MIGRATION_STATUS.md, "the CNTN3 issue turned out to be broader") --
        # the legacy cntn3_cds.bed/bta23_gwas_cds.bed mix exons from more than one
        # annotated isoform and do not translate to a clean single ORF for these
        # three genes.
        cntn3_bed="data/meta/cntn3_cds.v2.bed",
        gwas_bed="data/meta/bta23_gwas_cds.v2.bed",
    output:
        bed=f"{CODING_ROOT}/coding_cds_restrict.bed",
    log:
        f"{OUTPUT_ROOT}/logs/small_variants/build_coding_cds_bed.log"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/build_coding_cds_bed.py \
          --btn-cds-json {input.btn_json:q} \
          --cntn3-cds-bed {input.cntn3_bed:q} \
          --bta23-gwas-cds-bed {input.gwas_bed:q} \
          --output {output.bed:q} \
          > {log:q} 2>&1
        """


rule call_coding_variants:
    input:
        bam=f"{MAPPING_ROOT}/{{sample}}.targets.bam",
        bed=rules.build_coding_cds_bed.output.bed,
        reference=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        mapping_status=f"{MAPPING_ROOT}/{{sample}}.mapping_status.json",
    output:
        vcf=f"{CODING_ROOT}/{{sample}}.coding.vcf.gz",
        index=f"{CODING_ROOT}/{{sample}}.coding.vcf.gz.csi",
        status=f"{CODING_ROOT}/{{sample}}.coding_call_status.json",
        provenance=f"{CODING_ROOT}/provenance/{{sample}}/provenance.json",
    params:
        output_dir=CODING_ROOT,
        bcftools=config["tools"]["bcftools"],
    log:
        f"{OUTPUT_ROOT}/logs/small_variants/call_coding_variants/{{sample}}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/call_coding_variants.py \
          --mode call \
          --sample-id {wildcards.sample:q} \
          --bam {input.bam:q} \
          --reference {input.reference:q} \
          --cds-bed {input.bed:q} \
          --bcftools {params.bcftools:q} \
          --output-dir {params.output_dir:q} \
          > {log:q} 2>&1
        """


rule merge_coding_variants:
    input:
        vcfs=expand(f"{CODING_ROOT}/{{sample}}.coding.vcf.gz", sample=SAMPLES),
        registry=rules.input_audit.output.registry,
        reference=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        comparison=f"{REFERENCE_ROOT}/assembly_comparison/sequence_comparison.tsv",
        btn_json="data/meta/btn_cds_models.json",
        # v2: rebuilt from a single canonical GFF transcript per gene (see
        # UCD2_MIGRATION_STATUS.md, "the CNTN3 issue turned out to be broader") --
        # the legacy cntn3_cds.bed/bta23_gwas_cds.bed mix exons from more than one
        # annotated isoform and do not translate to a clean single ORF for these
        # three genes.
        cntn3_bed="data/meta/cntn3_cds.v2.bed",
        gwas_bed="data/meta/bta23_gwas_cds.v2.bed",
    output:
        table=f"{CODING_ROOT}/coding_variants.tsv",
        status=f"{CODING_ROOT}/coding_variants_status.json",
        provenance=f"{CODING_ROOT}/provenance.json",
        report="docs/codex_reanalysis/CODING_VARIANTS_REPORT.md",
    params:
        root=CODING_ROOT,
        bcftools=config["tools"]["bcftools"],
        samtools=config["tools"]["samtools"],
    log:
        f"{OUTPUT_ROOT}/logs/small_variants/merge_coding_variants.log"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/call_coding_variants.py \
          --mode merge \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --reference {input.reference:q} \
          --btn-cds-json {input.btn_json:q} \
          --cntn3-cds-bed {input.cntn3_bed:q} \
          --bta23-gwas-cds-bed {input.gwas_bed:q} \
          --assembly-comparison {input.comparison:q} \
          --bcftools {params.bcftools:q} \
          --samtools {params.samtools:q} \
          --output {output.report:q} \
          > {log:q} 2>&1
        """
