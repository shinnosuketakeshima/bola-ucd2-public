# ARS-UCD2.0 reference provenance protocol

The primary reference is RefSeq ARS-UCD2.0, assembly accession
`GCF_002263795.3`. The annotation is pinned to NCBI release
`GCF_002263795.3-RS_2024_12` (NCBI Eukaryotic Genome Annotation Pipeline
10.3, submitted 2024-12-26). Resources are retrieved only from the release-
specific NCBI FTP directory declared in `configs/codex_reanalysis.example.yaml`.

The registry downloads `md5checksums.txt` first and verifies the compressed
genomic FASTA, genomic GFF3, assembly report, and annotation report against
the official MD5 values. FASTA/GFF are atomically regenerated from those
verified inputs on every build, and FAI/MMI indexes are rebuilt into temporary
paths before publication. The manifest records the compressed resources,
expanded resources, indexes, official MD5 where available, locally calculated
SHA-256, byte size, URL or derivation, and retrieval date. Provenance also
records configuration checksum, software versions, Git commit, and dirty-tree
status. A mismatch is a hard stop.

All analytical gene intervals are regenerated from `gene` records in the
fixed GFF3 and the coordinate-free tracked target specification. That TSV is
deliberately rejected if it contains chromosome, start, end, or position
columns. GFF 1-based closed coordinates are converted to BED 0-based
half-open coordinates and the extracted interval SHA-256 is recorded.
Ambiguous aliases, duplicate gene matches, missing reference features, and
out-of-range extraction are hard failures.

`DQB2`, `NC10`, and `NC9` are declared as haplotype-specific probe targets.
If no standalone RefSeq gene feature exists, each remains a
`haplotype_probe` and is not written as an invented reference-gene BED
interval. Probe-panel alignment status is recorded independently using a
tracked exact target-to-panel-locus mapping: a generic DQB alignment cannot be
used as DQB2 evidence. The three manuscript SNPs are coordinate-free variant
probe targets until independently retrieved dbSNP sequences are aligned;
CNTN3, ABHD16A, and PRR3 are separately named host-gene targets.

An IPD alignment is accepted as panel evidence only when identity is at least
95%, query coverage is at least 80%, and MAPQ is at least 20. Lower-quality
hits are retained as `partial_alignment` evidence but do not set a target to
`panel_aligned`. Identity, coverage, and MAPQ are taken from the same ranked
alignment row. When multiple targets share one mapped panel locus (for example
the DQA or DQB families), the registry uses `panel_family_aligned`, explicitly
not a target-specific alignment claim.

The BTA22 and BTA23 records from ARS-UCD1.2 and ARS-UCD2.0 are compared by
accession, length, and SHA-256. The old accessions are read from the official,
MD5-verified ARS-UCD1.2 assembly report rather than inferred from `chr` FASTA
headers. Any difference receives
`decision_required`; it must not trigger silent coordinate reuse or liftover.
Sequence alignment evidence is required before a coordinate relationship is
accepted.

The operational build on 2026-07-31 resolved 51 unique RefSeq gene features
and retained six declared coordinate-free targets (DQB2, NC10, NC9, and the
three manuscript SNPs) without invented coordinates. The local ARS-UCD1.2 and
pinned ARS-UCD2.0 BTA22/BTA23 records
were identical in length and SHA-256 and had completed sequence alignments.
These observations are build outputs, not permission to reuse other legacy
coordinates.

Primary sources:

- NCBI annotation report: https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Bos_taurus/GCF_002263795.3-RS_2024_12/
- NCBI release FTP: https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9913/GCF_002263795.3-RS_2024_12/
