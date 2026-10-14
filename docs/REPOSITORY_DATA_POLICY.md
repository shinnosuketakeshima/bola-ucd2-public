# Repository data and privacy policy

This private repository contains code, small reviewed evidence files, and
manuscript-source material for the ARS-UCD2.0 reanalysis. It must never contain
source sequencing data, source alignment data, reference genomes, internal
anonymization keys, participant-identifying metadata, tool binaries, temporary
files, or unreviewed generated outputs.

## Protected paths and file types

Git ignores `data/raw/`, `data/meta/private/`, `results/`, `logs/`, `tmp/`,
`work/`, `.snakemake/`, local YAML configuration, and tool binaries. The only
potential exceptions under `results/` are reviewed files in
`results/codex_reanalysis_v1/07_evidence_audit/` and
`results/codex_reanalysis_v1/08_manuscript/`; they require the checks below
before staging.

Genomic formats and indexes (including FASTQ, BAM, CRAM, FASTA, VCF, BCF, and
their common compressed/indexed forms) are prohibited everywhere. So is
`docs/ANONYMIZATION_KEY_INTERNAL.md`.

## Staging gate

Before every commit, run:

```powershell
python -m pytest tests/test_repository_policy.py -q
git status --short
git diff --cached --name-only
```

Do not stage a file larger than 20 MiB, a symlink/junction that resolves into
`data/raw/` or the protected delivery root, or any protected format or path. Evidence/manuscript outputs may
be staged only after their size, privacy, and content review is documented.

The analysis receives its protected data root only at runtime, via the
`--deliver-root` command-line argument. Use
`configs/codex_reanalysis.example.yaml` as a portable template; local YAML
files are ignored.
