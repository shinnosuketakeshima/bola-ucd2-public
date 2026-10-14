"""Cohort QC for DRB3–rs29026690 wide-interval local assemblies."""

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

from scripts.codex_reanalysis.assemble_interval import WIDE_REGION
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file
from scripts.codex_reanalysis.summarize_mapping_qc import validate_cohort_ids
from scripts.codex_reanalysis.wide_linkage import (
    DRB3_REGION_FLANK,
    RS29026690_FLANK,
    MapHit,
    extract_drb3_panel_sequences,
    summarize_linkage,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_analytical_samples(registry: Path) -> list[dict[str, str]]:
    rows = read_tsv(registry)
    return [row for row in rows if row.get("analysis_status", "analytical") == "analytical"]


def parse_sam_hits(text: str) -> list[MapHit]:
    hits: list[MapHit] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        hits.append(
            MapHit(
                query=fields[0],
                contig=fields[2],
                pos=int(fields[3]),
                mapq=int(fields[4]),
            )
        )
    return hits


def map_queries_to_assembly(
    assembly: Path,
    queries_fasta: Path,
    *,
    minimap2: str,
    samtools: str,
    threads: int,
) -> list[MapHit]:
    with tempfile.NamedTemporaryFile(suffix=".sam", delete=False) as handle:
        sam_path = Path(handle.name)
    try:
        map_proc = subprocess.run(
            [
                minimap2,
                "-a",
                "-x",
                "map-ont",
                "-t",
                str(threads),
                str(assembly),
                str(queries_fasta),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if map_proc.returncode != 0:
            return []
        sam_path.write_text(map_proc.stdout, encoding="utf-8")
        view_proc = subprocess.run(
            [samtools, "view", "-F", "0x900", str(sam_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if view_proc.returncode != 0:
            return []
        return parse_sam_hits(view_proc.stdout)
    finally:
        sam_path.unlink(missing_ok=True)


def extract_reference_flank(reference: Path, region: str, samtools: str) -> str:
    proc = subprocess.run(
        [samtools, "faidx", str(reference), region],
        capture_output=True,
        text=True,
        check=True,
    )
    return proc.stdout


def summarize_sample(
    *,
    sample_id: str,
    assembly: Path,
    panel_text: str,
    reference: Path,
    hap1: str,
    hap2: str,
    status: dict,
    minimap2: str,
    samtools: str,
    threads: int,
) -> dict[str, str]:
    allele_fasta = extract_drb3_panel_sequences(panel_text, hap1, hap2)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, encoding="utf-8") as handle:
        allele_path = Path(handle.name)
        handle.write(allele_fasta)
    rs_fasta = extract_reference_flank(reference, RS29026690_FLANK, samtools)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False, encoding="utf-8") as handle:
        rs_path = Path(handle.name)
        handle.write(rs_fasta)
    try:
        allele_hits = (
            map_queries_to_assembly(
                assembly,
                allele_path,
                minimap2=minimap2,
                samtools=samtools,
                threads=threads,
            )
            if allele_fasta.strip()
            else []
        )
        rs_hits = map_queries_to_assembly(
            assembly,
            rs_path,
            minimap2=minimap2,
            samtools=samtools,
            threads=threads,
        )
    finally:
        allele_path.unlink(missing_ok=True)
        rs_path.unlink(missing_ok=True)

    row = summarize_linkage(hap1, hap2, allele_hits, rs_hits)
    row["sample_id"] = sample_id
    row["status"] = status.get("status", "failed")
    row["n_reads"] = str(status.get("n_reads", 0))
    row["n_contigs"] = str(status.get("n_contigs", 0))
    row["n50"] = str(status.get("n50", 0))
    return row


def render_report(rows: Sequence[dict[str, str]], interval: str = WIDE_REGION) -> str:
    bridged = sum(1 for row in rows if row["same_contig"] == "yes")
    lines = [
        "# ARS-UCD2.0 Wide-Interval Local-Assembly Report (DRB3–rs29026690)",
        "",
        f"Interval: `{interval}` (DRB3 through class III on ARS-UCD2.0 BTA23). "
        "Reads were extracted from the UCD2.0 target BAM (primary alignments, `-F 0x900`). "
        "Assembler: Flye 2.9.6 `--nano-raw --keep-haplotypes --min-overlap 3000 --genome-size 2.3m`.",
        "",
        "Linkage is assessed by mapping PCR-SBT DRB3 allele sequences from the IPD-MHC panel "
        f"and an rs29026690 flanking window (`{RS29026690_FLANK}`) onto each assembly with "
        f"`minimap2 -x map-ont`. DRB3 region flank `{DRB3_REGION_FLANK}` is recorded for "
        "provenance but not used in the cohort table.",
        "",
        f"Cohort: {len(rows)} animals. DRB3 and rs29026690 on the same contig: {bridged}/{len(rows)}.",
        "",
        "## Per-sample linkage",
        "",
        "| sample_id | PCR-SBT DRB3 | status | reads | contigs | N50 | DRB3 contig | rs29026690 contig | same contig | DRB3 resolution | note |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample_id']} | {row['drb3_diplotype']} | {row['status']} | {row['n_reads']} | "
            f"{row['n_contigs']} | {row['n50']} | {row['drb3_contig']} | {row['rs29026690_contig']} | "
            f"{row['same_contig']} | {row['drb3_resolution']} | {row['note']} |"
        )
    lines.extend(
        [
            "",
            "No animal bridged DRB3 and rs29026690 on one contig, consistent with the legacy "
            "WhatsHap phase-block discontinuity and the earlier UCD1.2 wide assembly (Table 7).",
            "",
        ]
    )
    return "\n".join(lines)


def write_linkage_tsv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    fieldnames = (
        "sample_id",
        "drb3_diplotype",
        "drb3_contig",
        "rs29026690_contig",
        "same_contig",
        "drb3_resolution",
        "note",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path("results/codex_reanalysis_v1/01_reference/GCF_002263795.3_ARS-UCD2.0_genomic.fna"),
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("data/raw/references/ipd_mhc_bola_panel.fasta"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    samples = load_analytical_samples(args.registry)
    sample_ids = [row["sample_id"] for row in samples]
    by_id = {row["sample_id"]: row for row in samples}
    panel_text = args.panel.read_text(encoding="utf-8")

    rows = []
    hashes = {}
    for sample_id in sample_ids:
        sample_dir = args.root / sample_id
        status_path = sample_dir / f"{sample_id}.assembly_status.json"
        assembly_path = sample_dir / "assembly.fasta"
        if not status_path.is_file():
            raise SystemExit(f"missing assembly status for {sample_id}: {status_path}")
        status = json.loads(status_path.read_text(encoding="utf-8"))
        registry_row = by_id[sample_id]
        if status.get("status") != "passed" or not assembly_path.is_file():
            rows.append(
                {
                    "sample_id": sample_id,
                    "drb3_diplotype": (
                        f"{registry_row.get('pcr_sbt_drb3_hap1', '')}/"
                        f"{registry_row.get('pcr_sbt_drb3_hap2', '')}"
                    ),
                    "drb3_contig": "-",
                    "rs29026690_contig": "-",
                    "same_contig": "no",
                    "drb3_resolution": "failed",
                    "note": "assembly failed",
                    "status": status.get("status", "failed"),
                    "n_reads": str(status.get("n_reads", 0)),
                    "n_contigs": str(status.get("n_contigs", 0)),
                    "n50": str(status.get("n50", 0)),
                }
            )
            hashes[sample_id] = sha256_file(status_path)
            continue
        row = summarize_sample(
            sample_id=sample_id,
            assembly=assembly_path,
            panel_text=panel_text,
            reference=args.reference,
            hap1=registry_row.get("pcr_sbt_drb3_hap1", ""),
            hap2=registry_row.get("pcr_sbt_drb3_hap2", ""),
            status=status,
            minimap2=args.minimap2,
            samtools=args.samtools,
            threads=args.threads,
        )
        rows.append(row)
        hashes[sample_id] = sha256_file(status_path)

    validate_cohort_ids(sample_ids, [row["sample_id"] for row in rows])
    passed = sum(1 for row in rows if row["status"] == "passed")
    report_text = render_report(rows)
    args.output.write_text(report_text, encoding="utf-8")
    write_linkage_tsv(args.table, rows)

    cohort_status = {
        "stage": "wide_assembly_qc",
        "status": "passed" if passed == len(rows) else "failed",
        "n_samples": len(rows),
        "n_passed": passed,
        "same_contig_count": sum(1 for row in rows if row["same_contig"] == "yes"),
        "sample_hashes": hashes,
    }
    status_path = args.root / "assembly_status.json"
    status_path.write_text(
        json.dumps(cohort_status, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        args.root,
        {
            "stage": "wide_assembly_qc",
            "input_hashes": hashes,
            "reference_hash": sha256_file(args.reference),
            "panel_hash": sha256_file(args.panel),
        },
    )
    print(f"wide assembly report written: {args.output} ({passed}/{len(rows)} passed)")
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
