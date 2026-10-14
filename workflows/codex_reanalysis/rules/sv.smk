"""Sniffles2 structural-variant calling and cross-sample evidence merge.

Standalone (`sv_all`), not wired into `evidence_audit`: no `envs/sv.yaml`
build has been exercised through `snakemake --use-conda` yet, so wiring it
into the shared evidence-audit target would make that target unrunnable for
every other stage until it has been.
"""


rule sniffles_call:
    input:
        bam=rules.map_targets.output.bam,
        bai=rules.map_targets.output.bai,
    output:
        vcf=f"{SV_ROOT}/{{sample}}/{{sample}}.sniffles2.vcf",
        status=f"{SV_ROOT}/{{sample}}/{{sample}}.sniffles_status.json",
        provenance=f"{SV_ROOT}/{{sample}}/provenance.json",
    params:
        output_dir=lambda wildcards, output: str(Path(output.status).parent),
        sniffles=config["tools"]["sniffles"],
        sniffles_version=config["tools"]["sniffles_version"],
    threads:
        config.get("mapping", {}).get("threads", 8)
    log:
        f"{OUTPUT_ROOT}/logs/sv/{{sample}}.log"
    conda:
        "../envs/sv.yaml"
    shell:
        r"""
        {params.sniffles:q} --input {input.bam:q} --vcf {output.vcf:q} --threads {threads} \
          --contig NC_037349.1 --contig NC_037350.1 --allow-overwrite \
          > {log:q} 2>&1
        python scripts/codex_reanalysis/merge_sv_evidence.py --mode record-call \
          --sample-id {wildcards.sample:q} \
          --vcf {output.vcf:q} \
          --sniffles-version {params.sniffles_version:q} \
          --output-dir {params.output_dir:q} \
          >> {log:q} 2>&1
        """


rule merge_sv_evidence:
    input:
        vcfs=expand(f"{SV_ROOT}/{{sample}}/{{sample}}.sniffles2.vcf", sample=SAMPLES),
        statuses=expand(f"{SV_ROOT}/{{sample}}/{{sample}}.sniffles_status.json", sample=SAMPLES),
        registry=rules.input_audit.output.registry,
    output:
        table=f"{SV_ROOT}/sv_evidence_table.tsv",
        status=f"{SV_ROOT}/sv_evidence_status.json",
        provenance=f"{SV_ROOT}/provenance.json",
        report="docs/codex_reanalysis/SV_REPORT.md",
    params:
        root=SV_ROOT,
    log:
        f"{OUTPUT_ROOT}/logs/sv/merge.log"
    conda:
        "../envs/sv.yaml"
    shell:
        r"""
        python scripts/codex_reanalysis/merge_sv_evidence.py --mode merge \
          --root {params.root:q} \
          --registry {input.registry:q} \
          --output {output.report:q} \
          > {log:q} 2>&1
        """
