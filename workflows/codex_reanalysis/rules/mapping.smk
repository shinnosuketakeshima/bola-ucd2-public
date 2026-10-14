"""Input/reference gates and streaming target-only mapping rules."""


rule input_audit:
    input:
        config="configs/codex_reanalysis.example.yaml",
        samplesheet="data/meta/samplesheet.tsv",
        mapping="data/meta/deliverable_mapping.tsv",
        drb3="data/meta/drb3_alleles.tsv",
    output:
        status=f"{AUDIT_ROOT}/audit_status.json",
        manifest=f"{AUDIT_ROOT}/input_manifest.tsv",
        qc=f"{AUDIT_ROOT}/fastq_qc.tsv",
        registry=f"{AUDIT_ROOT}/sample_registry.tsv",
        provenance=f"{AUDIT_ROOT}/provenance.json",
    params:
        output_root=lambda wildcards, output: str(Path(output.status).parent),
    log:
        f"{OUTPUT_ROOT}/logs/input_audit.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        test -n "${{CODEX_DELIVER_ROOT:-}}"
        python scripts/codex_reanalysis/audit_inputs.py \
          --deliver-root "$CODEX_DELIVER_ROOT" \
          --config {input.config:q} \
          --output {params.output_root:q} \
          > {log:q} 2>&1
        """


rule normalize_inputs:
    input:
        audit_status=rules.input_audit.output.status,
        qc=rules.input_audit.output.qc,
        manifest=rules.input_audit.output.manifest,
    output:
        status=f"{NORMALIZED_ROOT}/normalization_status.json",
        manifest=f"{NORMALIZED_ROOT}/normalized_manifest.tsv",
    params:
        output_root=lambda wildcards, output: str(Path(output.status).parent),
        audit_root=lambda wildcards, input: str(Path(input.audit_status).parent),
        minimum_read_length=MINIMUM_READ_LENGTH,
    log:
        f"{OUTPUT_ROOT}/logs/normalize_inputs.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        test -n "${{CODEX_DELIVER_ROOT:-}}"
        python scripts/codex_reanalysis/normalize_inputs.py \
          --deliver-root "$CODEX_DELIVER_ROOT" \
          --audit-dir {params.audit_root:q} \
          --output {params.output_root:q} \
          --minimum-read-length {params.minimum_read_length} \
          --compressor pigz \
          > {log:q} 2>&1
        """


rule reference_registry:
    input:
        config="configs/codex_reanalysis.example.yaml",
        targets="configs/codex_reanalysis.targets.tsv",
        legacy_fasta="data/raw/references/ARS-UCD1.2.fa",
        legacy_report="data/raw/references/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt",
        legacy_checksums="data/raw/references/GCF_002263795.1_ARS-UCD1.2_md5checksums.txt",
        ipd_manifest="data/raw/references/ipd_mhc_bola_manifest.tsv",
        ipd_panel="data/raw/references/ipd_mhc_bola_panel.fasta",
    output:
        status=f"{REFERENCE_ROOT}/reference_status.json",
        manifest=f"{REFERENCE_ROOT}/reference_manifest.tsv",
        bed=f"{REFERENCE_ROOT}/{ANALYTICAL_BED_NAME}",
        mapping_bed=f"{REFERENCE_ROOT}/{MAPPING_BED_NAME}",
        fasta=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        fai=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna.fai",
        mmi=f"{REFERENCE_ROOT}/GCF_002263795.3_ARS-UCD2.0_genomic.fna.mmi",
    params:
        output_root=lambda wildcards, output: str(Path(output.status).parent),
    log:
        f"{OUTPUT_ROOT}/logs/reference_registry.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/reference_registry.py \
          --config {input.config:q} \
          --output {params.output_root:q} \
          --legacy-reference {input.legacy_fasta:q} \
          --legacy-assembly-report {input.legacy_report:q} \
          --legacy-checksum-manifest {input.legacy_checksums:q} \
          --ipd-manifest {input.ipd_manifest:q} \
          --ipd-panel {input.ipd_panel:q} \
          > {log:q} 2>&1
        """


rule mapping_reference_seal:
    input:
        status=rules.reference_registry.output.status,
        manifest=rules.reference_registry.output.manifest,
        reference=rules.reference_registry.output.fasta,
        targets=rules.reference_registry.output.mapping_bed,
        config="configs/codex_reanalysis.example.yaml",
    output:
        seal=f"{REFERENCE_ROOT}/mapping_reference_seal.json",
    log:
        f"{OUTPUT_ROOT}/logs/mapping_reference_seal.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/stream_target_mapping.py seal-reference \
          --reference {input.reference:q} \
          --targets-bed {input.targets:q} \
          --reference-manifest {input.manifest:q} \
          --config {input.config:q} \
          --output {output.seal:q} \
          > {log:q} 2>&1
        """


rule map_targets:
    input:
        audit_status=rules.input_audit.output.status,
        sample_registry=rules.input_audit.output.registry,
        input_manifest=rules.input_audit.output.manifest,
        normalized_manifest=rules.normalize_inputs.output.manifest,
        reference_status=rules.reference_registry.output.status,
        reference_manifest=rules.reference_registry.output.manifest,
        reference_seal=rules.mapping_reference_seal.output.seal,
        reference=rules.reference_registry.output.fasta,
        targets=rules.reference_registry.output.mapping_bed,
        config="configs/codex_reanalysis.example.yaml",
    output:
        bam=f"{MAPPING_ROOT}/{{sample}}.targets.bam",
        bai=f"{MAPPING_ROOT}/{{sample}}.targets.bam.bai",
        flagstat=f"{MAPPING_ROOT}/{{sample}}.flagstat.txt",
        mapq=f"{MAPPING_ROOT}/{{sample}}.mapq_summary.tsv",
        status=f"{MAPPING_ROOT}/{{sample}}.mapping_status.json",
        commands=f"{MAPPING_ROOT}/{{sample}}.commands.json",
        provenance=f"{MAPPING_ROOT}/provenance/{{sample}}/provenance.json",
    threads:
        config.get("mapping", {}).get("threads", 8)
    params:
        output_root=lambda wildcards, output: str(Path(output.status).parent),
        minimap2=config["tools"]["minimap2"],
        samtools=config["tools"]["samtools"],
        secondary_max=config["mapping"]["secondary_max"],
    log:
        f"{OUTPUT_ROOT}/logs/mapping/{{sample}}.workflow.log"
    conda:
        "../envs/mapping.yaml"
    shell:
        r"""
        test -n "${{CODEX_DELIVER_ROOT:-}}"
        # The input seal is produced on the execution host, not by this
        # workflow, so it is supplied at runtime like the delivery root. When it
        # is present the delivered hash is a fact verified here; when it is
        # absent provenance says so rather than implying a verification.
        if [ -n "${{CODEX_INPUT_SEAL:-}}" ]; then
          set -- --input-seal "$CODEX_INPUT_SEAL"
        else
          set --
        fi
        python scripts/codex_reanalysis/stream_target_mapping.py map \
          --sample-id {wildcards.sample:q} \
          --sample-registry {input.sample_registry:q} \
          --input-manifest {input.input_manifest:q} \
          --normalized-manifest {input.normalized_manifest:q} \
          --deliver-root "$CODEX_DELIVER_ROOT" \
          --reference {input.reference:q} \
          --targets-bed {input.targets:q} \
          --audit-status {input.audit_status:q} \
          --reference-status {input.reference_status:q} \
          --reference-manifest {input.reference_manifest:q} \
          --reference-seal {input.reference_seal:q} \
          --config {input.config:q} \
          --output-dir {params.output_root:q} \
          --threads {threads} \
          --minimap2 {params.minimap2:q} \
          --samtools {params.samtools:q} \
          --secondary-max {params.secondary_max} \
          "$@" \
          > {log:q} 2>&1
        """
