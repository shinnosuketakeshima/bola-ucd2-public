"""Haplotype-kept local assembly of corrected BoLA intervals.

Standalone (`assembly_all` / `class_i_all`), not wired into `evidence_audit`:
Flye lives in its own prefix so the mapping lock is unchanged, and the
class I redo must finish and be reviewed before it can gate other stages.
"""


rule assemble_class_i:
    input:
        bam=rules.map_targets.output.bam,
        bai=rules.map_targets.output.bai,
        panel="data/raw/references/ipd_classI_panel.fasta",
    output:
        fasta=f"{ASSEMBLY_ROOT}/class_i/{{sample}}/assembly.fasta",
        paf=f"{ASSEMBLY_ROOT}/class_i/{{sample}}/classI.paf",
        status=f"{ASSEMBLY_ROOT}/class_i/{{sample}}/{{sample}}.assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/class_i/{{sample}}/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
        flye=config["tools"].get("flye", "flye"),
        samtools=config["tools"]["samtools"],
        minimap2=config["tools"]["minimap2"],
    threads:
        4
    log:
        f"{OUTPUT_ROOT}/logs/assembly/class_i/{{sample}}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assemble_interval.py \
          --sample-id {wildcards.sample:q} \
          --bam {input.bam:q} \
          --panel {input.panel:q} \
          --output-dir {params.output_dir:q} \
          --interval class_i \
          --flye {params.flye:q} \
          --samtools {params.samtools:q} \
          --minimap2 {params.minimap2:q} \
          --threads {threads} \
          > {log:q} 2>&1
        """


rule assembly_qc_class_i:
    input:
        statuses=expand(
            f"{ASSEMBLY_ROOT}/class_i/{{sample}}/{{sample}}.assembly_status.json",
            sample=SAMPLES,
        ),
        registry=rules.input_audit.output.registry,
    output:
        table=f"{ASSEMBLY_ROOT}/class_i/classI_recovery.tsv",
        copies=f"{ASSEMBLY_ROOT}/class_i/classI_copies.tsv",
        status=f"{ASSEMBLY_ROOT}/class_i/assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/class_i/provenance.json",
        report="docs/codex_reanalysis/ASSEMBLY_REPORT.md",
    params:
        root=f"{ASSEMBLY_ROOT}/class_i",
    log:
        f"{OUTPUT_ROOT}/logs/assembly/class_i/qc.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assembly_qc.py \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --output {output.report:q} \
          > {log:q} 2>&1
        """


rule assemble_dq_dr:
    input:
        bam=rules.map_targets.output.bam,
        bai=rules.map_targets.output.bai,
        panel="data/raw/references/ipd_mhc_bola_panel.fasta",
    output:
        fasta=f"{ASSEMBLY_ROOT}/dq_dr/{{sample}}/assembly.fasta",
        paf=f"{ASSEMBLY_ROOT}/dq_dr/{{sample}}/dq_dr.paf",
        status=f"{ASSEMBLY_ROOT}/dq_dr/{{sample}}/{{sample}}.assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/dq_dr/{{sample}}/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
        flye=config["tools"].get("flye", "flye"),
        samtools=config["tools"]["samtools"],
        minimap2=config["tools"]["minimap2"],
    threads:
        8
    log:
        f"{OUTPUT_ROOT}/logs/assembly/dq_dr/{{sample}}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assemble_interval.py \
          --sample-id {wildcards.sample:q} \
          --bam {input.bam:q} \
          --panel {input.panel:q} \
          --output-dir {params.output_dir:q} \
          --interval dq_dr \
          --flye {params.flye:q} \
          --samtools {params.samtools:q} \
          --minimap2 {params.minimap2:q} \
          --threads {threads} \
          > {log:q} 2>&1
        """


rule assembly_qc_dq_dr:
    input:
        statuses=expand(
            f"{ASSEMBLY_ROOT}/dq_dr/{{sample}}/{{sample}}.assembly_status.json",
            sample=SAMPLES,
        ),
        registry=rules.input_audit.output.registry,
    output:
        table=f"{ASSEMBLY_ROOT}/dq_dr/dq_dr_recovery.tsv",
        copies=f"{ASSEMBLY_ROOT}/dq_dr/dq_dr_copies.tsv",
        status=f"{ASSEMBLY_ROOT}/dq_dr/assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/dq_dr/provenance.json",
        report="docs/codex_reanalysis/ASSEMBLY_DQ_DR_REPORT.md",
    params:
        root=f"{ASSEMBLY_ROOT}/dq_dr",
    log:
        f"{OUTPUT_ROOT}/logs/assembly/dq_dr/qc.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assembly_qc_dq_dr.py \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --output {output.report:q} \
          > {log:q} 2>&1
        """


rule assemble_wide:
    input:
        bam=rules.map_targets.output.bam,
        bai=rules.map_targets.output.bai,
        panel="data/raw/references/ipd_mhc_bola_panel.fasta",
    output:
        fasta=f"{ASSEMBLY_ROOT}/wide/{{sample}}/assembly.fasta",
        paf=f"{ASSEMBLY_ROOT}/wide/{{sample}}/wide.paf",
        status=f"{ASSEMBLY_ROOT}/wide/{{sample}}/{{sample}}.assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/wide/{{sample}}/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
        flye=config["tools"].get("flye", "flye"),
        samtools=config["tools"]["samtools"],
        minimap2=config["tools"]["minimap2"],
    threads:
        8
    log:
        f"{OUTPUT_ROOT}/logs/assembly/wide/{{sample}}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assemble_interval.py \
          --sample-id {wildcards.sample:q} \
          --bam {input.bam:q} \
          --panel {input.panel:q} \
          --output-dir {params.output_dir:q} \
          --interval wide \
          --flye {params.flye:q} \
          --samtools {params.samtools:q} \
          --minimap2 {params.minimap2:q} \
          --threads {threads} \
          > {log:q} 2>&1
        """


rule assembly_qc_wide:
    input:
        statuses=expand(
            f"{ASSEMBLY_ROOT}/wide/{{sample}}/{{sample}}.assembly_status.json",
            sample=SAMPLES,
        ),
        registry=rules.input_audit.output.registry,
        panel="data/raw/references/ipd_mhc_bola_panel.fasta",
    output:
        table=f"{ASSEMBLY_ROOT}/wide/wide_linkage.tsv",
        status=f"{ASSEMBLY_ROOT}/wide/assembly_status.json",
        provenance=f"{ASSEMBLY_ROOT}/wide/provenance.json",
        report="docs/codex_reanalysis/ASSEMBLY_WIDE_REPORT.md",
    params:
        root=f"{ASSEMBLY_ROOT}/wide",
        reference=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        panel="data/raw/references/ipd_mhc_bola_panel.fasta",
    log:
        f"{OUTPUT_ROOT}/logs/assembly/wide/qc.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/assembly_qc_wide.py \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --reference {params.reference:q} \
          --panel {params.panel:q} \
          --output {output.report:q} \
          --table {output.table:q} \
          > {log:q} 2>&1
        """
