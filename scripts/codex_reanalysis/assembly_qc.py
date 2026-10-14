"""Cohort QC and report for haplotype-kept local assemblies.

Cluster counts, locus presence, and best-allele lineage are recovery
summaries. They are not validated diploid copy numbers and not official
IPD 4-digit SBT calls.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.assemble_interval import CLASS_I_REGION
from scripts.codex_reanalysis.probe_assignment import (
    CopyCluster,
    cluster_alignments,
    is_classical,
    locus_structure,
    parse_paf,
)
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file
from scripts.codex_reanalysis.summarize_mapping_qc import validate_cohort_ids


COPY_NUMBER_CAVEAT = (
    "Classical cluster counts are haplotype-contig recovery summaries, not a "
    "diploid copy number: `--keep-haplotypes` can split one gene across contigs, "
    "and collapse can merge two alleles onto one contig."
)

def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_analytical_samples(registry_path: Path) -> list[dict[str, str]]:
    return [row for row in read_tsv(registry_path) if row.get("analysis_status", "").strip() == "analytical"]


def summarize_sample(
    *,
    sample_id: str,
    copies: Sequence[CopyCluster],
    n_reads: int,
    n_bases: int,
    contig_lengths: Sequence[int],
    status: str,
    hap1: str,
    hap2: str,
    estimated_coverage: float | None = None,
    n50_value: int | None = None,
) -> dict:
    classical = [copy for copy in copies if is_classical(copy.best.locus)]
    nonclassical = [copy for copy in copies if not is_classical(copy.best.locus)]
    loci = sorted({copy.best.locus for copy in classical})
    best_by_locus: dict[str, CopyCluster] = {}
    for copy in classical:
        previous = best_by_locus.get(copy.best.locus)
        if previous is None or (copy.best.nmatch, copy.best.identity) > (
            previous.best.nmatch,
            previous.best.identity,
        ):
            best_by_locus[copy.best.locus] = copy
    structures = {locus: locus_structure(copies, locus) for locus in loci}
    best_alleles = "; ".join(
        f"{cluster.best.allele}@{cluster.best.identity:.2f}"
        for _locus, cluster in sorted(best_by_locus.items())
    )
    lengths = [int(length) for length in contig_lengths]
    return {
        "sample_id": sample_id,
        "pcr_sbt": f"{hap1}/{hap2}",
        "status": status,
        "n_reads": n_reads,
        "n_bases": n_bases,
        "estimated_coverage": estimated_coverage,
        "n_contigs": len(lengths),
        "total_assembly_length": sum(lengths),
        "n50": n50_value if n50_value is not None else (max(lengths) if lengths else 0),
        "longest_contig": max(lengths) if lengths else 0,
        "n_classical_clusters": len(classical),
        "n_nonclassical_clusters": len(nonclassical),
        "n_classical_loci": len(loci),
        "classical_loci": ",".join(loci) if loci else "-",
        "best_alleles": best_alleles,
        "structures": structures,
        "hap1": hap1,
        "hap2": hap2,
    }


def render_report(summaries: Sequence[dict], interval: str = CLASS_I_REGION) -> str:
    lines = [
        "# ARS-UCD2.0 Class I Local-Assembly Report",
        "",
        f"Interval: `{interval}` (NC10 through classical BOLA on ARS-UCD2.0 BTA23). "
        "Reads were extracted from the UCD2.0 target BAM (primary alignments, `-F 0x900`), not from the delivery drive. "
        "Assembler: Flye 2.9.6 `--nano-raw --keep-haplotypes --min-overlap 3000 --genome-size 0.35m`. "
        "IPD class I CDS panel (145 alleles) mapped with `minimap2 -x splice -p 0.2 -N 30 --secondary=yes`.",
        "",
        COPY_NUMBER_CAVEAT
        + " Best-allele labels are the highest-identity IPD match at ≥95% identity, ≥80% query coverage, "
        "MAPQ ≥20; they are not official 4-digit SBT genotypes.",
        "",
        "Wide-interval assembly (DRB3–rs29026690) is reported in "
        "`docs/codex_reanalysis/ASSEMBLY_WIDE_REPORT.md`.",
        "",
        f"Cohort: {len(summaries)} animals.",
        "",
        "## Per-sample recovery",
        "",
        "| sample_id | PCR-SBT DRB3 | status | reads | contigs | N50 | classical clusters | loci | best allele (identity) | structure |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in summaries:
        structures = row.get("structures") or {}
        structure_text = (
            ", ".join(f"{locus}={label}" for locus, label in sorted(structures.items())) or "-"
        )
        lines.append(
            f"| {row['sample_id']} | {row['pcr_sbt']} | {row['status']} | {row['n_reads']} | "
            f"{row['n_contigs']} | {row['n50']} | {row['n_classical_clusters']} | "
            f"{row['classical_loci']} | {row['best_alleles'] or '-'} | {structure_text} |"
        )
    present = {
        locus: sum(
            1
            for row in summaries
            if locus in (row.get("classical_loci") or "").split(",")
        )
        for locus in ("BoLA-1", "BoLA-2", "BoLA-3", "BoLA-6")
    }
    lines.extend(
        [
            "",
            "## Locus presence (high-confidence clusters)",
            "",
            f"- BoLA-2: {present['BoLA-2']}/{len(summaries)}",
            f"- BoLA-3: {present['BoLA-3']}/{len(summaries)}",
            f"- BoLA-1: {present['BoLA-1']}/{len(summaries)}",
            f"- BoLA-6: {present['BoLA-6']}/{len(summaries)}",
            "",
            "Locus presence is a recovery summary under the project IPD evidence rule "
            "(0.95 / 0.80 / MAPQ 20) on primary target-BAM reads, not manuscript copy-number evidence. "
            "Paralogous panel alleles that hit the same contig interval are assigned the single best locus; "
            "they do not count as extra genes. Absences (e.g., ID05 missing BoLA-2/6; "
            "ID08 and ID13 missing BoLA-3) are assembly/threshold outcomes, "
            "not evidence that the gene is deleted.",
            "",
        ]
    )
    return "\n".join(lines)


def write_recovery_tsv(path: Path, summaries: Sequence[dict]) -> None:
    fieldnames = (
        "sample_id",
        "pcr_sbt",
        "status",
        "n_reads",
        "n_bases",
        "estimated_coverage",
        "n_contigs",
        "total_assembly_length",
        "n50",
        "longest_contig",
        "n_classical_clusters",
        "n_nonclassical_clusters",
        "n_classical_loci",
        "classical_loci",
        "best_alleles",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)


def write_copies_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "sample_id",
        "contig",
        "ref_start",
        "ref_end",
        "class",
        "locus",
        "best_allele",
        "identity",
        "coverage",
        "mapq",
        "n_distinct_alleles",
        "structure",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="06_local_assembly/class_i directory")
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Markdown report path")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    samples = load_analytical_samples(args.registry)
    sample_ids = [row["sample_id"] for row in samples]
    by_id = {row["sample_id"]: row for row in samples}

    summaries = []
    copy_rows = []
    hashes = {}
    for sample_id in sample_ids:
        sample_dir = args.root / sample_id
        status_path = sample_dir / f"{sample_id}.assembly_status.json"
        paf_path = sample_dir / "classI.paf"
        if not status_path.is_file():
            raise SystemExit(f"missing assembly status for {sample_id}: {status_path}")
        status = json.loads(status_path.read_text(encoding="utf-8"))
        paf_text = paf_path.read_text(encoding="utf-8") if paf_path.is_file() else ""
        copies = cluster_alignments(parse_paf(paf_text))
        registry_row = by_id[sample_id]
        summary = summarize_sample(
            sample_id=sample_id,
            copies=copies,
            n_reads=int(status.get("n_reads", 0)),
            n_bases=int(status.get("n_bases", 0)),
            contig_lengths=[status.get("longest_contig", 0)] if status.get("n_contigs") else [],
            status=status.get("status", "failed"),
            hap1=registry_row.get("pcr_sbt_drb3_hap1", ""),
            hap2=registry_row.get("pcr_sbt_drb3_hap2", ""),
            estimated_coverage=status.get("estimated_coverage"),
            n50_value=status.get("n50"),
        )
        if status.get("n_contigs"):
            summary["n_contigs"] = int(status["n_contigs"])
            summary["total_assembly_length"] = int(status.get("total_assembly_length", 0))
            summary["longest_contig"] = int(status.get("longest_contig", 0))
        summaries.append(summary)
        hashes[sample_id] = sha256_file(status_path)
        for copy in copies:
            copy_rows.append(
                {
                    "sample_id": sample_id,
                    "contig": copy.contig,
                    "ref_start": copy.tstart,
                    "ref_end": copy.tend,
                    "class": "classical" if is_classical(copy.best.locus) else "non-classical",
                    "locus": copy.best.locus,
                    "best_allele": copy.best.allele,
                    "identity": round(copy.best.identity, 4),
                    "coverage": round(copy.best.coverage, 4),
                    "mapq": copy.best.mapq,
                    "n_distinct_alleles": copy.n_distinct_alleles,
                    "structure": locus_structure(copies, copy.best.locus),
                }
            )

    validate_cohort_ids(sample_ids, [row["sample_id"] for row in summaries])
    args.root.mkdir(parents=True, exist_ok=True)
    write_recovery_tsv(args.root / "classI_recovery.tsv", summaries)
    write_copies_tsv(args.root / "classI_copies.tsv", copy_rows)
    payload = {
        "sample_count": len(summaries),
        "passed_count": sum(1 for row in summaries if row["status"] == "passed"),
        "interval": CLASS_I_REGION,
        "status": "passed" if all(row["status"] == "passed" for row in summaries) else "failed",
    }
    (args.root / "assembly_status.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(args.root, {"input_hashes": hashes})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(summaries, interval=CLASS_I_REGION), encoding="utf-8")
    print(f"assembly report written: {args.output} ({payload['passed_count']}/{payload['sample_count']} passed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
