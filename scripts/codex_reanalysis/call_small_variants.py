"""Genotype one marker position from a target BAM by two independent pileups.

No Clair3 model verified against this cohort's Guppy5/R9.4.1 chemistry is
available on this host (see
``docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md`` §5.5),
so every genotype produced here is exploratory: it is a BCFtools
allele-depth pileup cross-checked against an independently coded raw-pileup
base counter, not a caller/model result. The 15%/85% allele-fraction dosage
rule below is the same one validated at 98.5% agreement against the
official caller for 65 sample x position combinations in
``docs/CLASSIII_016_TAG_SNP_EXPLORATION.md``, reused here rather than
invented fresh.
"""

from __future__ import annotations

import subprocess
from typing import Sequence


class PileupParseError(ValueError):
    """Raised when a pileup line cannot be parsed at the requested position."""


def build_bcftools_mpileup_command(
    bam: str,
    ref_fasta: str,
    contig: str,
    position: int,
    bcftools: str = "bcftools",
    mapq: int = 20,
    baseq: int = 20,
) -> list[str]:
    region = f"{contig}:{position}-{position}"
    return [
        bcftools,
        "mpileup",
        "-f",
        ref_fasta,
        "-r",
        region,
        "-a",
        "AD,DP",
        "--no-BAQ",
        "-q",
        str(mapq),
        "-Q",
        str(baseq),
        "-Ov",
        bam,
    ]


def build_samtools_mpileup_command(
    bam: str,
    ref_fasta: str,
    contig: str,
    position: int,
    samtools: str = "samtools",
    mapq: int = 20,
    baseq: int = 20,
) -> list[str]:
    region = f"{contig}:{position}-{position}"
    return [samtools, "mpileup", "-f", ref_fasta, "-r", region, "-q", str(mapq), "-Q", str(baseq), bam]


def parse_bcftools_ad_vcf(text: str) -> dict:
    """Parse the one data line bcftools mpileup -a AD emits for a single-position region."""
    data_line = next((line for line in text.splitlines() if line and not line.startswith("#")), None)
    if data_line is None:
        return {"ref_base": None, "alt_base": None, "ref_count": 0, "alt_count": 0, "depth": 0}
    fields = data_line.split("\t")
    ref_base, alt_field, format_keys, sample_value = fields[3], fields[4], fields[8], fields[9]
    alleles = [ref_base] + [allele for allele in alt_field.split(",") if allele not in {".", "<*>"}]
    key_index = {key: index for index, key in enumerate(format_keys.split(":"))}
    if "AD" not in key_index:
        raise PileupParseError("bcftools mpileup output is missing an AD field")
    ad_values = [int(value) for value in sample_value.split(":")[key_index["AD"]].split(",")]
    ad_values = ad_values[: len(alleles)]
    ref_count = ad_values[0] if ad_values else 0
    alt_counts = ad_values[1:]
    alt_index = max(range(len(alt_counts)), default=None, key=lambda i: alt_counts[i]) if alt_counts else None
    alt_base = alleles[1 + alt_index] if alt_index is not None and alt_counts and alt_counts[alt_index] > 0 else None
    return {
        "ref_base": ref_base,
        "alt_base": alt_base,
        "ref_count": ref_count,
        "alt_count": sum(alt_counts),
        "depth": ref_count + sum(alt_counts),
    }


def parse_samtools_pileup_line(text: str) -> dict:
    """Parse one classic ``samtools mpileup`` line, counting substitutions only.

    Insertions, deletions and skipped/reference-skip bases are excluded from
    both the ref and alt counts, matching the indel-avoidance policy already
    validated in ``docs/CLASSIII_016_TAG_SNP_EXPLORATION.md`` for this noisy
    ONT R9 chemistry.
    """
    line = next((candidate for candidate in text.splitlines() if candidate.strip()), None)
    if line is None:
        return {"ref_base": None, "alt_base": None, "ref_count": 0, "alt_count": 0, "depth": 0}
    fields = line.split("\t")
    ref_base = fields[2].upper()
    bases = fields[4]

    ref_count = 0
    alt_counts: dict[str, int] = {}
    i = 0
    while i < len(bases):
        char = bases[i]
        if char == "^":
            i += 2
            continue
        if char == "$":
            i += 1
            continue
        if char in "+-":
            i += 1
            digits = ""
            while i < len(bases) and bases[i].isdigit():
                digits += bases[i]
                i += 1
            i += int(digits)
            continue
        if char in ".,":
            ref_count += 1
        elif char == "*":
            pass
        elif char.upper() in "ACGT":
            alt_counts[char.upper()] = alt_counts.get(char.upper(), 0) + 1
        i += 1

    alt_base = max(alt_counts, key=alt_counts.get) if alt_counts else None
    alt_count = alt_counts[alt_base] if alt_base else 0
    return {
        "ref_base": ref_base,
        "alt_base": alt_base,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": ref_count + alt_count,
    }


def classify_dosage(ref_count: int, alt_count: int) -> str:
    """Apply the project's validated 15%/85% allele-fraction dosage rule."""
    total = ref_count + alt_count
    if total == 0:
        return "no_data"
    alt_fraction = alt_count / total
    if alt_fraction < 0.15:
        return "ref"
    if alt_fraction > 0.85:
        return "alt"
    return "het"


def genotype_label(ref_base: str | None, alt_base: str | None, dosage_class: str) -> str:
    if dosage_class == "no_data" or ref_base is None:
        return "no_data"
    if dosage_class == "ref":
        return f"{ref_base}/{ref_base}"
    if dosage_class == "alt" and alt_base is not None:
        return f"{alt_base}/{alt_base}"
    if dosage_class == "het" and alt_base is not None:
        return "/".join(sorted([ref_base, alt_base]))
    return "no_data"


def run_pileup(command: Sequence[str]) -> str:
    result = subprocess.run(list(command), capture_output=True, check=True, text=True)
    return result.stdout


def genotype_marker(
    bam: str,
    ref_fasta: str,
    contig: str,
    position: int,
    bcftools: str = "bcftools",
    samtools: str = "samtools",
    mapq: int = 20,
    baseq: int = 20,
) -> dict:
    """Genotype one marker with both independent pileup methods."""
    bcftools_text = run_pileup(
        build_bcftools_mpileup_command(bam, ref_fasta, contig, position, bcftools, mapq, baseq)
    )
    primary = parse_bcftools_ad_vcf(bcftools_text)
    primary_dosage = classify_dosage(primary["ref_count"], primary["alt_count"])
    primary_genotype = genotype_label(primary["ref_base"], primary["alt_base"], primary_dosage)

    samtools_text = run_pileup(
        build_samtools_mpileup_command(bam, ref_fasta, contig, position, samtools, mapq, baseq)
    )
    pileup = parse_samtools_pileup_line(samtools_text)
    pileup_dosage = classify_dosage(pileup["ref_count"], pileup["alt_count"])
    pileup_genotype = genotype_label(pileup["ref_base"], pileup["alt_base"], pileup_dosage)

    return {
        "contig": contig,
        "position": position,
        "primary_genotype": primary_genotype,
        "primary_depth": primary["depth"],
        "primary_ref_count": primary["ref_count"],
        "primary_alt_count": primary["alt_count"],
        "pileup_genotype": pileup_genotype,
        "pileup_depth": pileup["depth"],
        "pileup_ref_count": pileup["ref_count"],
        "pileup_alt_count": pileup["alt_count"],
    }
