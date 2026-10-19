"""Re-genotype the CNTN3/ABHD16A/PRR3 missense sites in every animal.

The discovery calls (`call_coding_variants.py`, `call -mv`) emit a record only
where an alternative allele was called, so a missing record is not a reference
genotype. This stage revisits every candidate missense site in every animal
and assigns a genotype only where the evidence supports one:

* allele depths come from `bcftools mpileup` with the pinned coding protocol
  (`-B -d 60 -q 20`, base quality `--min-base-quality`, default 10);
* the allele-fraction genotype uses the project's 15%/85% rule on
  ref + alt reads (`call_small_variants.classify_dosage`);
* the model genotype comes from `bcftools call -m -C alleles` at the same
  site in the same animal (single-sample model, as in discovery);
* a genotype is `called` only when ref + alt depth >= `--min-depth`, both
  genotypes agree, and a non-reference call has GQ >= `--min-gq`.
  `bcftools call` reports no GQ for homozygous-reference records at forced
  sites, so 0/0 calls rest on depth, allele fraction and agreement.

Everything else is `unresolved` with its reason and counts towards no
denominator. Frequencies are computed over called alleles only.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.call_small_variants import classify_dosage
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

GENES = {"CNTN3": "NC_037349.1", "ABHD16A": "NC_037350.1", "PRR3": "NC_037350.1"}
DOSAGE_GT = {"ref": "0/0", "het": "0/1", "alt": "1/1"}
HIGH, INTERMEDIATE, LOW = "high", "intermediate", "low"


class CandidateSiteError(ValueError):
    """Raised when inputs cannot support a genotype table."""


def load_candidate_sites(coding_variants: Path) -> list[dict]:
    sites = []
    with coding_variants.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene"] in GENES and row["effect"] == "MISSENSE":
                sites.append(
                    {
                        "gene": row["gene"],
                        "contig": GENES[row["gene"]],
                        "position": int(row["pos"]),
                        "ref": row["ref"],
                        "alt": row["alt"],
                        "detail": row["detail"],
                    }
                )
    if not sites:
        raise CandidateSiteError(f"no CNTN3/ABHD16A/PRR3 missense sites in {coding_variants}")
    return sorted(sites, key=lambda site: (site["contig"], site["position"]))


def frequency_class(alt_af: float) -> str:
    if alt_af >= 0.85:
        return HIGH
    if alt_af >= 0.15:
        return INTERMEDIATE
    return LOW


def judge_genotype(
    ref_count: int,
    alt_count: int,
    caller_gt: str,
    caller_gq: int | None,
    min_depth: int,
    min_gq: int,
) -> tuple[str, str, str]:
    """Return (status, genotype, reason) for one animal at one site."""
    depth = ref_count + alt_count
    if depth == 0:
        return "unresolved", "./.", "no_data"
    fraction_gt = DOSAGE_GT[classify_dosage(ref_count, alt_count)]
    if depth < min_depth:
        return "unresolved", "./.", "depth_below_minimum"
    if caller_gt != fraction_gt:
        return "unresolved", "./.", "caller_disagrees_with_allele_fraction"
    if fraction_gt != "0/0" and (caller_gq is None or caller_gq < min_gq):
        return "unresolved", "./.", "genotype_quality_below_minimum"
    return "called", fraction_gt, ""


def summarize_site(site: dict, calls: list[dict]) -> dict:
    called = [call for call in calls if call["status"] == "called"]
    het = sum(1 for call in called if call["genotype"] == "0/1")
    hom_alt = sum(1 for call in called if call["genotype"] == "1/1")
    called_alleles = 2 * len(called)
    alt_alleles = het + 2 * hom_alt
    alt_af = alt_alleles / called_alleles if called_alleles else None
    return {
        **{key: site[key] for key in ("gene", "contig", "position", "ref", "alt", "detail")},
        "animals_called": len(called),
        "animals_unresolved": len(calls) - len(called),
        "carriers": het + hom_alt,
        "het": het,
        "hom_alt": hom_alt,
        "alt_alleles": alt_alleles,
        "called_alleles": called_alleles,
        "alt_af": "" if alt_af is None else f"{alt_af:.3f}",
        "frequency_class": "" if alt_af is None else frequency_class(alt_af),
        "depth_min": min(call["depth"] for call in calls),
        "depth_max": max(call["depth"] for call in calls),
    }


def _parse_query_value(value: str) -> int | None:
    return None if value in {"", "."} else int(value)


def genotype_sample(
    bam: Path,
    reference: Path,
    sites_file: Path,
    sites: list[dict],
    bcftools: str,
    min_base_quality: int,
    workdir: Path,
) -> dict[tuple[str, int], dict]:
    """Return allele counts and the caller genotype for every site in one BAM."""
    pileup = workdir / f"{bam.stem}.pileup.bcf"
    subprocess.run(
        [
            bcftools, "mpileup", "-B", "-d", "60", "-q", "20", "-Q", str(min_base_quality),
            "-a", "FORMAT/AD,FORMAT/DP", "-f", str(reference), "-T", str(sites_file),
            "-Ob", "-o", str(pileup), str(bam),
        ],
        check=True,
        capture_output=True,
    )
    counts_text = subprocess.run(
        [bcftools, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n", str(pileup)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    called = workdir / f"{bam.stem}.called.bcf"
    subprocess.run(
        [
            bcftools, "call", "-m", "-C", "alleles", "-T", str(sites_file), "-a", "GQ",
            "-Ob", "-o", str(called), str(pileup),
        ],
        check=True,
        capture_output=True,
    )
    call_text = subprocess.run(
        [bcftools, "query", "-f", "%CHROM\t%POS\t%REF[\t%GT\t%GQ]\n", str(called)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    wanted = {(site["contig"], site["position"]): site for site in sites}
    result: dict[tuple[str, int], dict] = {}
    for line in counts_text.splitlines():
        contig, position, ref, alts, ad = line.split("\t")
        key = (contig, int(position))
        site = wanted.get(key)
        if site is None or ref != site["ref"]:
            continue  # indel record at the same position
        alleles = alts.split(",")
        depths = [_parse_query_value(value) or 0 for value in ad.split(",")]
        alt_count = depths[1 + alleles.index(site["alt"])] if site["alt"] in alleles else 0
        result[key] = {"ref_count": depths[0], "alt_count": alt_count, "caller_gt": "./.", "caller_gq": None}
    for line in call_text.splitlines():
        contig, position, ref, gt, gq = line.split("\t")
        key = (contig, int(position))
        if key in result and ref == wanted[key]["ref"]:
            result[key]["caller_gt"] = gt.replace("|", "/")
            result[key]["caller_gq"] = _parse_query_value(gq)
    for key in wanted:
        result.setdefault(key, {"ref_count": 0, "alt_count": 0, "caller_gt": "./.", "caller_gq": None})
    return result


def write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coding-variants", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True, help="00_input_audit/sample_registry.tsv")
    parser.add_argument("--bam-dir", type=Path, required=True, help="02_mapping (<sample>.targets.bam)")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--min-depth", type=int, default=10)
    parser.add_argument("--min-gq", type=int, default=20)
    parser.add_argument("--min-base-quality", type=int, default=10)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    sites = load_candidate_sites(args.coding_variants)
    with args.registry.open(encoding="utf-8", newline="") as handle:
        samples = [
            row["sample_id"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row["analysis_status"] == "analytical"
        ]
    if len(samples) != 13:
        raise CandidateSiteError(f"expected 13 analytical samples, found {len(samples)}")
    bams = {sample: args.bam_dir / f"{sample}.targets.bam" for sample in samples}
    for sample, bam in bams.items():
        if not bam.is_file():
            raise CandidateSiteError(f"missing target BAM for {sample}: {bam}")

    args.output.mkdir(parents=True, exist_ok=True)
    genotype_rows: list[dict] = []
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        sites_file = workdir / "sites.tsv"
        sites_file.write_text(
            "".join(f"{s['contig']}\t{s['position']}\t{s['ref']},{s['alt']}\n" for s in sites),
            encoding="utf-8",
        )
        for sample in samples:
            evidence = genotype_sample(
                bams[sample], args.reference, sites_file, sites, args.bcftools, args.min_base_quality, workdir
            )
            for site in sites:
                found = evidence[(site["contig"], site["position"])]
                status, genotype, reason = judge_genotype(
                    found["ref_count"], found["alt_count"], found["caller_gt"], found["caller_gq"],
                    args.min_depth, args.min_gq,
                )
                depth = found["ref_count"] + found["alt_count"]
                genotype_rows.append(
                    {
                        "gene": site["gene"],
                        "contig": site["contig"],
                        "position": site["position"],
                        "ref": site["ref"],
                        "alt": site["alt"],
                        "sample_id": sample,
                        "depth": depth,
                        "ref_count": found["ref_count"],
                        "alt_count": found["alt_count"],
                        "alt_fraction": f"{found['alt_count'] / depth:.3f}" if depth else "",
                        "caller_gt": found["caller_gt"],
                        "caller_gq": "" if found["caller_gq"] is None else found["caller_gq"],
                        "status": status,
                        "genotype": genotype,
                        "reason": reason,
                    }
                )

    summary_rows = [
        summarize_site(
            site,
            [
                row for row in genotype_rows
                if (row["contig"], row["position"]) == (site["contig"], site["position"])
            ],
        )
        for site in sites
    ]
    write_tsv(args.output / "candidate_site_genotypes.tsv", genotype_rows)
    write_tsv(args.output / "candidate_site_summary.tsv", summary_rows)
    parameters = {
        "min_depth": args.min_depth,
        "min_gq": args.min_gq,
        "min_base_quality": args.min_base_quality,
        "min_mapping_quality": 20,
        "max_depth_per_file": 60,
        "allele_fraction_rule": "ref <0.15, het 0.15-0.85, alt >0.85 over ref+alt reads",
        "frequency_classes": "high >=0.85, intermediate 0.15-0.85, low <0.15 over called alleles",
    }
    status = {
        "status": "passed",
        "samples": len(samples),
        "sites": len(sites),
        "called_genotypes": sum(1 for row in genotype_rows if row["status"] == "called"),
        "unresolved_genotypes": sum(1 for row in genotype_rows if row["status"] != "called"),
        "parameters": parameters,
    }
    (args.output / "candidate_sites_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(
        args.output,
        {
            "parameters": parameters,
            "input_hashes": {
                "coding_variants": sha256_file(args.coding_variants),
                **{f"bam:{sample}": sha256_file(bam) for sample, bam in bams.items()},
            },
        },
    )
    print(f"Genotyped {len(sites)} sites in {len(samples)} animals: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
