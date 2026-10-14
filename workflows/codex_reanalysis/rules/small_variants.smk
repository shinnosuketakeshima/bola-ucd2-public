"""Published-marker coordinate resolution and two-pileup genotyping."""


rule resolve_marker_coordinates:
    input:
        markers="configs/codex_reanalysis.markers.tsv",
        comparison=f"{REFERENCE_ROOT}/assembly_comparison/sequence_comparison.tsv",
        reference_status=f"{REFERENCE_ROOT}/reference_status.json",
    output:
        registry=f"{MARKERS_ROOT}/coordinates/marker_coordinate_registry.tsv",
        status=f"{MARKERS_ROOT}/coordinates/marker_coordinate_status.json",
        provenance=f"{MARKERS_ROOT}/coordinates/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
    log:
        f"{OUTPUT_ROOT}/logs/small_variants/resolve_coordinates.log"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/resolve_marker_coordinates.py \
          --markers {input.markers:q} \
          --assembly-comparison {input.comparison:q} \
          --output {params.output_dir:q} \
          > {log:q} 2>&1
        """


rule audit_published_markers:
    input:
        registry=rules.resolve_marker_coordinates.output.registry,
        status=rules.resolve_marker_coordinates.output.status,
        samples=rules.input_audit.output.registry,
        bams=expand(f"{MAPPING_ROOT}/{{sample}}.targets.bam", sample=SAMPLES),
        mapping_status=expand(f"{MAPPING_ROOT}/{{sample}}.mapping_status.json", sample=SAMPLES),
        reference=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        seal=f"{REFERENCE_ROOT}/mapping_reference_seal.json",
    output:
        matrix=f"{MARKERS_ROOT}/marker_genotypes.tsv",
        status=f"{MARKERS_ROOT}/marker_audit_status.json",
        provenance=f"{MARKERS_ROOT}/provenance.json",
        report="docs/codex_reanalysis/MARKER_REPORT.md",
    params:
        mapping_root=MAPPING_ROOT,
        output_dir=MARKERS_ROOT,
        bcftools=config["tools"]["bcftools"],
        samtools=config["tools"]["samtools"],
        mapq=config["thresholds"]["minimum_mapping_quality"],
        baseq=config["thresholds"]["minimum_base_quality"],
        minimum_depth=config["thresholds"]["minimum_depth"],
    log:
        f"{OUTPUT_ROOT}/logs/small_variants/audit_markers.log"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/audit_published_markers.py \
          --markers {input.registry:q} \
          --registry {input.samples:q} \
          --mapping-root {params.mapping_root:q} \
          --reference {input.reference:q} \
          --reference-seal {input.seal:q} \
          --bcftools {params.bcftools:q} \
          --samtools {params.samtools:q} \
          --mapq {params.mapq} \
          --baseq {params.baseq} \
          --minimum-depth {params.minimum_depth} \
          --output {params.output_dir:q} \
          --report {output.report:q} \
          > {log:q} 2>&1
        """
