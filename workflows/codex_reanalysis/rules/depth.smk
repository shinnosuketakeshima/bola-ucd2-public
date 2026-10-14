"""MAPQ 0/10/20 depth extraction and cohort-normalized ratios."""


rule compute_depth:
    input:
        bam=rules.map_targets.output.bam,
        bai=rules.map_targets.output.bai,
        bed=rules.reference_registry.output.bed,
    output:
        raw=f"{DEPTH_ROOT}/{{sample}}/depth_raw.tsv",
        status=f"{DEPTH_ROOT}/{{sample}}/depth_status.json",
        provenance=f"{DEPTH_ROOT}/{{sample}}/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
        samtools=config["tools"]["samtools"],
    log:
        f"{OUTPUT_ROOT}/logs/depth/{{sample}}.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/compute_depth.py \
          --sample-id {wildcards.sample:q} \
          --bam {input.bam:q} \
          --bed {input.bed:q} \
          --samtools {params.samtools:q} \
          --output {params.output_dir:q} \
          > {log:q} 2>&1
        """


rule normalize_depth:
    input:
        raw=expand(f"{DEPTH_ROOT}/{{sample}}/depth_raw.tsv", sample=SAMPLES),
        registry=rules.input_audit.output.registry,
    output:
        matrix=f"{DEPTH_ROOT}/normalized_depth_matrix.tsv",
        status=f"{DEPTH_ROOT}/normalized_depth_status.json",
        provenance=f"{DEPTH_ROOT}/provenance.json",
        report="docs/codex_reanalysis/DEPTH_CNV_REPORT.md",
    params:
        root=DEPTH_ROOT,
    log:
        f"{OUTPUT_ROOT}/logs/depth/normalize.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/normalize_depth.py \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --output {output.report:q} \
          > {log:q} 2>&1
        """
