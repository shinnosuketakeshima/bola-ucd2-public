# Analysis code: BoLA region of Japanese Black cattle (ARS-UCD2.0)

Code for the Oxford Nanopore whole-genome reanalysis of the bovine MHC (BoLA) region in
13 Japanese Black cattle, mapped to RefSeq ARS-UCD2.0 (`GCF_002263795.3`).

## Data availability

The sequencing reads are **not** in this repository. They are deposited at the DDBJ Sequence Read
Archive (DRA032007) under BioProject PRJDB45880, BioSamples SAMD02115956-SAMD02115968,
Experiments DRX1182054-DRX1182066 and Runs DRR1213128-DRR1213140. The records are held private
until publication and are mirrored to NCBI SRA and ENA after release. Samples are identified only by
the manuscript IDs ID01-ID13.

## Contents

- `scripts/codex_reanalysis/`: pipeline stages (input audit, reference registry, streaming target
  mapping, depth/CNV, DRB3 exon-2 typing, evidence matrix, DDBJ worksheets).
- `workflows/codex_reanalysis/`: Snakemake rules and conda environments.
- `run_remote/`: guarded scripts for the remote mapping host.
- `configs/`: configuration template, target and marker lists, sample-sheet schema.
- `docs/codex_reanalysis/ddbj_metadata/`: the public DDBJ BioProject, BioSample, Experiment and Run
  worksheets and the accession register.
- `tests/`: unit tests, including the repository data policy checks.

## Reproducing

The pipeline needs the delivered FASTQ files, which are supplied at run time through `--deliver-root`
(or `CODEX_DELIVER_ROOT` for Snakemake) and are never stored here. In this public copy the original
sample names are replaced by `ID01`-`ID13` and the vendor delivery numbers by `DELIVERY_nn`, so the
sample sheet of the original delivery must be mapped to these identifiers before it is used.

Python 3.10 or newer is required (3.13 is used in CI).

```bash
python -m pip install -r requirements.txt
python -m pytest tests -q
```

## Data policy

Genomic file formats, raw data and local absolute paths are never committed; see
`docs/REPOSITORY_DATA_POLICY.md`.

## Licence

MIT, see `LICENSE`. If you use the code, please cite it as in `CITATION.cff`.
