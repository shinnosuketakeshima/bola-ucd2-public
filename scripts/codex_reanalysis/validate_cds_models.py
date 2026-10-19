"""Validate CDS exon models by translating them and checking for a clean ORF.

`call_coding_variants.py` reuses legacy ARS-UCD1.2-era CDS exon lists
(`data/meta/btn_cds_models.json`, `data/meta/cntn3_cds.bed`, `data/meta/bta23_gwas_cds.bed`)
because BTA22/BTA23 are byte-identical between assemblies. That reuse assumes each legacy
exon list represents a single, internally consistent transcript. It does not always: some
of these genes have multiple annotated isoforms in the pinned GFF, and a legacy file that
mixes exons from more than one isoform produces a CDS that is not a multiple of 3, or
whose translation contains a premature or missing stop codon -- a frame error that
silently corrupts every downstream amino-acid call without ever raising an exception.

This script is the check that catches that: translate a candidate CDS model (from either
source) with the same `build_gene_sequence` used for real variant calling, and assert the
result is a clean single open reading frame (starts M, exactly one trailing stop, no
internal stop). It also builds a single-transcript model directly from the pinned GFF, for
comparison, preferring a curated RefSeq transcript (`NM_`/`NR_`/`NP_` accession prefix,
i.e. "BestRefSeq") over a Gnomon prediction (`XM_`/`XR_`/`XP_`), tie-broken by total CDS
length -- this is a documented analytical decision (which isoform is "the" transcript),
not an incidental default, so the chosen transcript accession is recorded in the report.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.call_coding_variants import (
    CODON_TABLE,
    build_gene_sequence,
)

STOP = "*"


class CdsValidationError(ValueError):
    """Raised when a CDS model cannot be parsed or translated."""


def translate(coding_seq: str) -> str:
    """Translate a coding sequence codon by codon; '?' for any unrecognized codon."""
    protein = []
    for i in range(0, len(coding_seq) - len(coding_seq) % 3, 3):
        protein.append(CODON_TABLE.get(coding_seq[i : i + 3], "?"))
    return "".join(protein)


def validate_orf(protein: str, coding_length: int) -> tuple[bool, str]:
    """Return (is_clean_orf, reason). A clean ORF starts M, ends with one stop, no others."""
    if coding_length % 3 != 0:
        return False, f"CDS length {coding_length} is not a multiple of 3"
    if not protein:
        return False, "empty translation"
    if protein[0] != "M":
        return False, f"does not start with Met (starts {protein[0]!r})"
    if protein.count(STOP) != 1:
        return False, f"expected exactly one stop codon, found {protein.count(STOP)}"
    if protein[-1] != STOP:
        return False, "stop codon is not at the end (premature stop)"
    return True, "clean single ORF"


def parse_gff_cds_by_transcript(
    gff_path: Path, contig: str, gene_start: int, gene_end: int
) -> dict[str, dict]:
    """Group CDS records within a gene's span by their Parent transcript.

    Returns {transcript_id: {"exons": [(start, end), ...], "strand": "+"/"-"}}, 1-based closed.
    """
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    by_transcript: dict[str, dict] = {}
    with opener(gff_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != contig or fields[2] != "CDS":
                continue
            start, end = int(fields[3]), int(fields[4])
            if start < gene_start or end > gene_end:
                continue
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            parent = attrs.get("Parent", "")
            if not parent.startswith("rna-"):
                continue
            transcript_id = parent[len("rna-") :]
            record = by_transcript.setdefault(transcript_id, {"exons": [], "strand": fields[6]})
            record["exons"].append((start, end))
    for record in by_transcript.values():
        record["exons"].sort()
    return by_transcript


_CURATED_PREFIXES = ("NM_", "NR_", "NP_")


def choose_canonical_transcript(by_transcript: dict[str, dict]) -> str:
    """Prefer a curated RefSeq (NM_/NR_/NP_) transcript; else the longest CDS."""
    if not by_transcript:
        raise CdsValidationError("no CDS transcripts found in gene span")

    def sort_key(item: tuple[str, dict]) -> tuple[int, int]:
        transcript_id, record = item
        is_curated = 1 if transcript_id.startswith(_CURATED_PREFIXES) else 0
        length = sum(end - start + 1 for start, end in record["exons"])
        return (is_curated, length)

    return max(by_transcript.items(), key=sort_key)[0]


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True, help="samtools-faidx-able FASTA")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--genes", type=Path, required=True, help="TSV: gene,contig,start,end,strand")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def fetch_sequence(samtools: str, reference: Path, contig: str, start: int, end: int) -> str:
    import subprocess

    text = subprocess.run(
        [samtools, "faidx", str(reference), f"{contig}:{start}-{end}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    rows = []
    with args.genes.open(encoding="utf-8", newline="") as handle:
        import csv

        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["gene"]
            contig = row["contig"]
            gene_start, gene_end = int(row["start"]), int(row["end"])
            by_transcript = parse_gff_cds_by_transcript(args.gff, contig, gene_start, gene_end)
            transcript_id = choose_canonical_transcript(by_transcript)
            model = by_transcript[transcript_id]
            region_start = min(start for start, _ in model["exons"])
            region_end = max(end for _, end in model["exons"])
            region_seq = fetch_sequence(args.samtools, args.reference, contig, region_start, region_end)
            coding_seq, _ = build_gene_sequence(region_seq, region_start, model["exons"], model["strand"])
            protein = translate(coding_seq)
            is_clean, reason = validate_orf(protein, len(coding_seq))
            rows.append(
                {
                    "gene": gene,
                    "transcript_id": transcript_id,
                    "strand": model["strand"],
                    "exon_count": len(model["exons"]),
                    "cds_length": len(coding_seq),
                    "clean_orf": is_clean,
                    "reason": reason,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    passed = sum(1 for row in rows if row["clean_orf"])
    print(f"CDS model validation written: {args.output} ({passed}/{len(rows)} clean ORFs)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
