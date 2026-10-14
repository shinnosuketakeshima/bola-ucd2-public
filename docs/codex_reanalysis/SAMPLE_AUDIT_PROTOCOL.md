# Delivered FASTQ input audit protocol

This audit is the gate before remapping.  It receives the delivery root only
as an explicit command-line argument and does not change, copy, link, or log
source files.  Its derived tables are written under the ignored
`results/codex_reanalysis_v1/00_input_audit/` directory.

For every primary-cohort FASTQ, the audit streams the gzip file, confirms
four-line FASTQ structure, checks sequence/quality length equality, calculates
non-sensitive length and base-quality summaries, and records a SHA-256 hash of
the compressed input.  It also detects same-sample and cross-sample duplicate
read IDs using a temporary database of one-way digests; no read ID is retained
after the audit stops.

The registry requires exactly one FASTQ for each of the thirteen manuscript
samples.  `DELIVERY_07` and `DELIVERY_10` are explicit absent-cohort rows, not
failures.  PCR-SBT DRB3 diplotypes are preserved only as external metadata and
are not treated as sequence-derived evidence.

Files with `.corrupt.bak` or `.rebuild.bak.*` names are excluded.  A malformed
gzip stream, a truncated FASTQ record, an inconsistent sample mapping, or a
missing expected analytical FASTQ is a blocking error.  A successful audit has
13 analytical rows plus the two absent-cohort rows.

Run with:

```powershell
python scripts/codex_reanalysis/audit_inputs.py --deliver-root H:/Deliver_ADD --config configs/codex_reanalysis.example.yaml --output results/codex_reanalysis_v1/00_input_audit
```
