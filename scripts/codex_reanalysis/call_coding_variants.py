"""Coding-variant screening at BTN loci and at CNTN3/ABHD16A/PRR3 on the ARS-UCD2.0 pipeline.

Ports `scripts/btn_coding_aggregate.py`'s translation/effect/frequency logic (the
`bcftools mpileup -B -d 60 -q 20 -Q 10` + `call -mv` protocol pinned in
`scripts/_btn_coding_call.sh`) onto target BAMs `02_mapping` already produced -- BTN,
CNTN3, ABHD16A, and PRR3 are all already full-gene targets in
`configs/codex_reanalysis.targets.tsv`, so no re-mapping is needed for this stage.

CDS exon spans are reused from the legacy ARS-UCD1.2-era files
(`data/meta/btn_cds_models.json`, `data/meta/cntn3_cds.bed`,
`data/meta/bta23_gwas_cds.bed`) because BTA22 (`NC_037349.1`, CNTN3) and BTA23
(`NC_037350.1`, everything else) are byte-identical between ARS-UCD1.2 and ARS-UCD2.0.
`verify_contig_identical` checks this against
`01_reference/assembly_comparison/sequence_comparison.tsv` before any exon span is
trusted -- the same discipline `resolve_marker_coordinates.py` already applies to marker
positions, not a new assumption introduced here.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

MISSENSE = "MISSENSE"
STOPGAIN = "STOPGAIN"
SYNONYMOUS = "synonymous"
INDEL = "indel"
NONCODING_IN_SPAN = "noncoding_in_span"
CDS_EDGE = "cds_edge"

BREED_DIVERGENT = "breed_divergent"
SEGREGATING = "segregating"
RARE = "rare"

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


class CodingVariantError(ValueError):
    """Raised when coding-variant screening cannot proceed safely."""


def _standard_codon_table() -> dict[str, str]:
    bases = "TCAG"
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    table = {}
    for i, amino_acid in enumerate(amino_acids):
        codon = "".join(bases[(i >> shift) & 3] for shift in (4, 2, 0))
        table[codon] = amino_acid
    return table


CODON_TABLE = _standard_codon_table()


def complement(base: str) -> str:
    return base.translate(_COMPLEMENT)


def verify_contig_identical(comparison_path: Path, contig_accession: str) -> bool:
    """True only if sequence_comparison.tsv records this contig's status as identical."""
    with Path(comparison_path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("new_accession") == contig_accession:
                return row.get("status") == "identical"
    return False


def load_cds_models_from_json(path: Path, strand_by_gene: dict[str, str]) -> dict[str, dict]:
    """Load data/meta/btn_cds_models.json's shape: {gene: {chrom, cds: [[start,end],...]}}."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    models = {}
    for gene, record in payload.items():
        if gene not in strand_by_gene:
            raise CodingVariantError(f"no strand declared for gene {gene!r} from {path}")
        models[gene] = {
            "contig": record["chrom"],
            "strand": strand_by_gene[gene],
            "exons": [tuple(exon) for exon in record["cds"]],
        }
    return models


def load_cds_models_from_bed(path: Path, strand_by_gene: dict[str, str] | None = None) -> dict[str, dict]:
    """Load a 4- or 5-column CDS BED (chrom, start0, end, gene[, strand]) into gene models.

    BED is 0-based half-open; exons are stored 1-based closed to match the JSON loader.
    """
    strand_by_gene = strand_by_gene or {}
    models: dict[str, dict] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            contig, start0, end, gene = fields[0], int(fields[1]), int(fields[2]), fields[3]
            strand = fields[4] if len(fields) >= 5 else strand_by_gene.get(gene)
            if strand not in ("+", "-"):
                raise CodingVariantError(f"no strand for gene {gene!r} in {path} (row: {line!r})")
            model = models.setdefault(gene, {"contig": contig, "strand": strand, "exons": []})
            if model["contig"] != contig:
                raise CodingVariantError(f"gene {gene!r} spans multiple contigs in {path}")
            model["exons"].append((start0 + 1, end))
    return models


def build_gene_sequence(region_seq: str, region_start: int, exons: Sequence[tuple[int, int]], strand: str) -> tuple[str, dict[int, int]]:
    """Return (coding_seq, {genomic_pos: coding_index}) from a fetched genomic region.

    `region_seq` is 1-based-inclusive-mapped: region_seq[0] is the base at `region_start`.
    """
    genomic_positions = []
    for start, end in sorted(exons):
        genomic_positions.extend(range(start, end + 1))
    if strand == "-":
        genomic_positions = genomic_positions[::-1]

    coding_bases = []
    for pos in genomic_positions:
        base = region_seq[pos - region_start]
        coding_bases.append(complement(base) if strand == "-" else base)
    coding_seq = "".join(coding_bases).upper()
    gene_index = {pos: i for i, pos in enumerate(genomic_positions)}
    return coding_seq, gene_index


def classify_effect(
    gene_seq: str, gene_index: dict[int, int], pos: int, ref: str, alt: str, strand: str
) -> tuple[str, str]:
    """Classify one SNV's coding effect. Indels and non-CDS positions are reported, not scored."""
    if len(ref) != 1 or len(alt) != 1:
        return INDEL, ""
    if pos not in gene_index:
        return NONCODING_IN_SPAN, ""
    i = gene_index[pos]
    codon_no, frame = divmod(i, 3)
    codon_start = codon_no * 3
    if codon_start + 3 > len(gene_seq):
        return CDS_EDGE, ""
    ref_codon = gene_seq[codon_start : codon_start + 3]
    alt_base = complement(alt) if strand == "-" else alt
    alt_codon = list(ref_codon)
    alt_codon[frame] = alt_base.upper()
    alt_codon = "".join(alt_codon)
    ref_aa = CODON_TABLE.get(ref_codon, "?")
    alt_aa = CODON_TABLE.get(alt_codon, "?")
    detail = f"{ref_codon}({ref_aa})>{alt_codon}({alt_aa}) codon{codon_no + 1}"
    if ref_aa == alt_aa:
        return SYNONYMOUS, detail
    if alt_aa == "*":
        return STOPGAIN, detail
    return MISSENSE, detail


def classify_allele_frequency(alt_af: float) -> str:
    """breed_divergent (near-fixed vs. Hereford) / segregating (candidate) / rare (ONT-error-prone)."""
    if alt_af >= 0.85:
        return BREED_DIVERGENT
    if alt_af >= 0.15:
        return SEGREGATING
    return RARE


def parse_vcf_calls(text: str) -> dict[int, tuple[str, str, str]]:
    """Parse a bcftools-call VCF body into {pos: (ref, first_alt, genotype)}."""
    calls: dict[int, tuple[str, str, str]] = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            raise CodingVariantError(f"VCF data line has fewer than 10 columns: {line!r}")
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4].split(",")[0]
        genotype = fields[9].split(":")[0].replace("|", "/")
        calls[pos] = (ref, alt, genotype)
    return calls


def gene_for_position(models: dict[str, dict], pos: int) -> str | None:
    for gene, model in models.items():
        for start, end in model["exons"]:
            if start <= pos <= end:
                return gene
    return None


def build_variant_table(models: dict[str, dict], sequences_by_gene: dict[str, tuple[str, dict[int, int]]], calls_by_sample: dict[str, dict[int, tuple[str, str, str]]]) -> list[dict]:
    """Aggregate per-sample VCF calls into one row per (gene, position, alt) across the cohort."""
    variants: dict[tuple[int, str], dict] = {}
    for sample_id, calls in calls_by_sample.items():
        for pos, (ref, alt, genotype) in calls.items():
            gene = gene_for_position(models, pos)
            if gene is None:
                continue
            key = (pos, alt)
            record = variants.setdefault(
                key, {"pos": pos, "ref": ref, "alt": alt, "gene": gene, "genotypes": {}}
            )
            record["genotypes"][sample_id] = genotype

    sample_count = len(calls_by_sample)
    rows = []
    for (pos, alt), record in variants.items():
        gene = record["gene"]
        gene_seq, gene_index = sequences_by_gene[gene]
        strand = models[gene]["strand"]
        genotypes = record["genotypes"]
        called = [gt for gt in genotypes.values() if gt not in ("./.", ".")]
        carriers = [gt for gt in genotypes.values() if gt not in ("0/0", "./.", ".")]
        homozygous_alt = [gt for gt in genotypes.values() if gt == "1/1"]
        n_called = max(len(called), 1)
        heterozygous_count = len(carriers) - len(homozygous_alt)
        alt_af = round((2 * len(homozygous_alt) + heterozygous_count) / (2 * n_called), 3)
        effect, detail = classify_effect(gene_seq, gene_index, pos, record["ref"], alt, strand)
        rows.append(
            {
                "gene": gene,
                "pos": pos,
                "ref": record["ref"],
                "alt": alt,
                "allele_frequency_class": classify_allele_frequency(alt_af),
                "alt_af": alt_af,
                "n_carriers": len(carriers),
                "n_hom_alt": len(homozygous_alt),
                "n_called": len(called),
                "sample_count": sample_count,
                "effect": effect,
                "detail": detail,
            }
        )
    return sorted(rows, key=lambda row: (row["gene"], row["pos"]))


def write_variant_table_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "gene",
        "pos",
        "ref",
        "alt",
        "allele_frequency_class",
        "alt_af",
        "n_carriers",
        "n_hom_alt",
        "n_called",
        "sample_count",
        "effect",
        "detail",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_report(rows: Sequence[dict], gene_names: Sequence[str]) -> str:
    lines = [
        "# Coding-Variant Screen (BTN + CNTN3/ABHD16A/PRR3, ARS-UCD2.0)",
        "",
        "`bcftools mpileup -B -d 60 -q 20 -Q 10` (BAQ disabled for ONT data) + `call -mv`, "
        "against the ARS-UCD2.0 target BAM, restricted to each gene's CDS exons. CDS exon "
        "spans are reused from ARS-UCD1.2-era files because BTA22/BTA23 are byte-identical "
        "between the two assemblies (verified per-contig before use, not assumed).",
        "",
        "Because this cohort derives from 2017-era R9/Guppy5 chemistry with a comparatively "
        "high base-calling error rate, `rare` (alt_af < 0.15) calls should be treated as "
        "ONT-error-prone; `breed_divergent` (>= 0.85) reflects a near-fixed difference from "
        "the Hereford reference rather than segregating diversity. `segregating` "
        "(0.15-0.85) MISSENSE/STOPGAIN calls are the Sanger-confirmation candidates.",
        "",
        f"Genes screened: {', '.join(gene_names)}. {len(rows)} coding variants across the cohort.",
        "",
        "| gene | n_variants | breed_divergent | segregating | rare | segregating MISSENSE | segregating STOPGAIN |",
        "|---|---|---|---|---|---|---|",
    ]
    for gene in gene_names:
        gene_rows = [row for row in rows if row["gene"] == gene]
        if not gene_rows:
            continue
        breed_divergent = sum(1 for row in gene_rows if row["allele_frequency_class"] == BREED_DIVERGENT)
        segregating = sum(1 for row in gene_rows if row["allele_frequency_class"] == SEGREGATING)
        rare = sum(1 for row in gene_rows if row["allele_frequency_class"] == RARE)
        segregating_missense = sum(
            1 for row in gene_rows if row["allele_frequency_class"] == SEGREGATING and row["effect"] == MISSENSE
        )
        segregating_stopgain = sum(
            1 for row in gene_rows if row["allele_frequency_class"] == SEGREGATING and row["effect"] == STOPGAIN
        )
        lines.append(
            f"| {gene} | {len(gene_rows)} | {breed_divergent} | {segregating} | {rare} | "
            f"{segregating_missense} | {segregating_stopgain} |"
        )
    lines.append("")
    return "\n".join(lines)


def build_mpileup_call_command(
    bcftools: str, reference: Path, bed: Path, bam: Path, output_vcf: Path
) -> str:
    return (
        f"{bcftools} mpileup -B -d 60 -q 20 -Q 10 -a FORMAT/AD,FORMAT/DP "
        f"-f {reference} -R {bed} -Ou {bam} 2>/dev/null | "
        f"{bcftools} call -mv -Oz -o {output_vcf}"
    )


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("call", "merge"), default="merge")
    parser.add_argument("--sample-id", help="call mode")
    parser.add_argument("--bam", type=Path, help="call mode")
    parser.add_argument("--reference", type=Path, help="both modes")
    parser.add_argument("--cds-bed", type=Path, help="call mode: BED restricting mpileup")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--output-dir", type=Path, help="call mode")
    parser.add_argument("--root", type=Path, help="merge mode: directory of per-sample VCFs")
    parser.add_argument("--registry", type=Path, help="merge mode: sample_registry.tsv")
    parser.add_argument("--btn-cds-json", type=Path, help="merge mode")
    parser.add_argument("--cntn3-cds-bed", type=Path, help="merge mode")
    parser.add_argument("--bta23-gwas-cds-bed", type=Path, help="merge mode")
    parser.add_argument("--assembly-comparison", type=Path, help="merge mode")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--output", type=Path, help="merge mode: report destination")
    return parser.parse_args(argv)


_BTN_STRAND = {
    "BTNL2": "+",
    "LOC104969893_BTN2A1L": "+",
    "LOC107131258_BTN2A1L": "+",
    "LOC615223_BTN3A2": "-",
    "LOC107131749_BTN2A2": "-",
    "LOC786760_BTN3A3": "-",
}
#: data/meta/cntn3_cds.bed's own gene column reads "CNTN3_CDS", not "CNTN3" -- loaded
#: under that name and renamed to "CNTN3" by load_cntn3_cds_models below so downstream
#: reporting (and the manuscript's own gene naming) stays "CNTN3".
_CNTN3_BED_GENE_NAME = "CNTN3_CDS"
_CNTN3_STRAND = {_CNTN3_BED_GENE_NAME: "+"}


def load_cntn3_cds_models(path: Path) -> dict[str, dict]:
    """Load data/meta/cntn3_cds.bed and rename its "CNTN3_CDS" gene column to "CNTN3"."""
    models = load_cds_models_from_bed(path, _CNTN3_STRAND)
    if _CNTN3_BED_GENE_NAME in models:
        models["CNTN3"] = models.pop(_CNTN3_BED_GENE_NAME)
    return models

#: bta23_gwas_cds.bed carries the whole class III region (ABCF1, APOM, BAG6, ...); only
#: ABHD16A and PRR3 are in scope for this stage. Named after the post-rename model keys
#: (load_cntn3_cds_models already renamed "CNTN3_CDS" to "CNTN3"), not the raw BED column.
_GENES_IN_SCOPE = frozenset({*_BTN_STRAND, "CNTN3", "ABHD16A", "PRR3"})


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)

    if args.mode == "call":
        for name, value in (
            ("--sample-id", args.sample_id),
            ("--bam", args.bam),
            ("--reference", args.reference),
            ("--cds-bed", args.cds_bed),
            ("--output-dir", args.output_dir),
        ):
            if value is None:
                raise SystemExit(f"call mode requires {name}")
        args.output_dir.mkdir(parents=True, exist_ok=True)
        output_vcf = args.output_dir / f"{args.sample_id}.coding.vcf.gz"
        command = build_mpileup_call_command(args.bcftools, args.reference, args.cds_bed, args.bam, output_vcf)
        subprocess.run(command, shell=True, check=True)
        subprocess.run([args.bcftools, "index", "-f", str(output_vcf)], check=True)
        status = {"sample_id": args.sample_id, "status": "passed" if output_vcf.is_file() else "failed"}
        (args.output_dir / f"{args.sample_id}.coding_call_status.json").write_text(
            json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_provenance(
            args.output_dir / "provenance" / args.sample_id,
            {"command": command, "input_hashes": {"bam": sha256_file(args.bam)} if args.bam.is_file() else {}},
        )
        print(f"Coding-variant call recorded for {args.sample_id}: {output_vcf}")
        return 0

    for name, value in (
        ("--root", args.root),
        ("--registry", args.registry),
        ("--reference", args.reference),
        ("--btn-cds-json", args.btn_cds_json),
        ("--cntn3-cds-bed", args.cntn3_cds_bed),
        ("--bta23-gwas-cds-bed", args.bta23_gwas_cds_bed),
        ("--assembly-comparison", args.assembly_comparison),
        ("--output", args.output),
    ):
        if value is None:
            raise SystemExit(f"merge mode requires {name}")

    models: dict[str, dict] = {}
    models.update(load_cds_models_from_json(args.btn_cds_json, _BTN_STRAND))
    models.update(load_cntn3_cds_models(args.cntn3_cds_bed))
    models.update(load_cds_models_from_bed(args.bta23_gwas_cds_bed))
    models = {gene: model for gene, model in models.items() if gene in _GENES_IN_SCOPE}

    verified_contigs = {model["contig"] for model in models.values()}
    for contig in verified_contigs:
        if not verify_contig_identical(args.assembly_comparison, contig):
            raise CodingVariantError(
                f"contig {contig!r} is not verified identical between ARS-UCD1.2 and ARS-UCD2.0; "
                "cannot reuse its legacy CDS exon spans"
            )

    with args.registry.open(encoding="utf-8", newline="") as handle:
        sample_ids = [
            row["sample_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("analysis_status", "").strip() == "analytical"
        ]

    calls_by_sample: dict[str, dict[int, tuple[str, str, str]]] = {}
    input_hashes = {}
    for sample_id in sample_ids:
        vcf_path = args.root / f"{sample_id}.coding.vcf.gz"
        if not vcf_path.is_file():
            continue
        text = subprocess.run(
            [args.bcftools, "view", "-H", str(vcf_path)], check=True, capture_output=True, text=True
        ).stdout
        calls_by_sample[sample_id] = parse_vcf_calls(text)
        input_hashes[sample_id] = sha256_file(vcf_path)

    sequences_by_gene: dict[str, tuple[str, dict[int, int]]] = {}
    for gene, model in models.items():
        region_start = min(start for start, _ in model["exons"])
        region_end = max(end for _, end in model["exons"])
        region_seq = subprocess.run(
            [args.samtools, "faidx", str(args.reference), f"{model['contig']}:{region_start}-{region_end}"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        sequence = "".join(line.strip() for line in region_seq.splitlines() if not line.startswith(">"))
        sequences_by_gene[gene] = build_gene_sequence(sequence, region_start, model["exons"], model["strand"])

    rows = build_variant_table(models, sequences_by_gene, calls_by_sample)

    args.root.mkdir(parents=True, exist_ok=True)
    write_variant_table_tsv(args.root / "coding_variants.tsv", rows)

    status_payload = {
        "sample_count": len(calls_by_sample),
        "gene_count": len(models),
        "variant_count": len(rows),
        "status": "passed",
    }
    (args.root / "coding_variants_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(args.root, {"input_hashes": input_hashes})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(rows, sorted(models)), encoding="utf-8")
    print(f"Coding-variant report written: {args.output} ({len(rows)} variants, {len(models)} genes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
