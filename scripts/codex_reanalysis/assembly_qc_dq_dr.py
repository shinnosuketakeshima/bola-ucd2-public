"""Cohort QC report for DQ/DR haplotype-kept local assemblies."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.assemble_interval import DQ_DR_REGION, paf_filename
from scripts.codex_reanalysis.dq_dr_recovery import (
    classify_dq_interpretation,
    classify_drb3_pcr_concordance,
    count_family_clusters,
    parse_ipd_classii_paf,
)
from scripts.codex_reanalysis.probe_assignment import cluster_alignments, locus_structure
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file
from scripts.codex_reanalysis.summarize_mapping_qc import validate_cohort_ids


COPY_NUMBER_CAVEAT = (
    "DQA and DQB cluster counts are locus-family recovery summaries, not validated "
    "diploid copy numbers: `--keep-haplotypes` can split one gene across contigs, and "
    "collapse can merge two alleles onto one contig."
)

def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_analytical_samples(registry_path: Path) -> list[dict[str, str]]:
    return [row for row in read_tsv(registry_path) if row.get("analysis_status", "").strip() == "analytical"]


def summarize_sample(
    *,
    sample_id: str,
    copies: Sequence,
    status: dict,
    hap1: str,
    hap2: str,
) -> dict:
    counts = count_family_clusters(copies)
    drb3_clusters = [copy for copy in copies if copy.best.locus == "DRB3"]
    best_drb3 = "; ".join(
        f"{copy.best.allele}@{copy.best.identity:.2f}"
        for copy in sorted(drb3_clusters, key=lambda item: (-item.best.identity, item.contig))
    )
    return {
        "sample_id": sample_id,
        "pcr_sbt": f"{hap1}/{hap2}",
        "status": status.get("status", "failed"),
        "n_reads": int(status.get("n_reads", 0)),
        "n_contigs": int(status.get("n_contigs", 0)),
        "n50": int(status.get("n50", 0)),
        "n_dqa_clusters": counts["DQA"],
        "n_dqb_clusters": counts["DQB"],
        "n_drb3_clusters": counts["DRB3"],
        "interpretation": classify_dq_interpretation(
            counts["DQA"], counts["DQB"], flye_status=status.get("status", "failed")
        ),
        "drb3_concordance": classify_drb3_pcr_concordance(hap1, hap2, drb3_clusters),
        "best_drb3_alleles": best_drb3,
        "carrier_016_01": "*016:01" in {hap1, hap2},
    }


def render_report(summaries: Sequence[dict], interval: str = DQ_DR_REGION) -> str:
    carriers = [row for row in summaries if row["carrier_016_01"]]
    noncarriers = [row for row in summaries if not row["carrier_016_01"]]
    carrier_dqa = [row["n_dqa_clusters"] for row in carriers]
    carrier_dqb = [row["n_dqb_clusters"] for row in carriers]
    noncarrier_dqa = [row["n_dqa_clusters"] for row in noncarriers]
    noncarrier_dqb = [row["n_dqb_clusters"] for row in noncarriers]

    lines = [
        "# ARS-UCD2.0 DQ/DR Local-Assembly Report",
        "",
        f"Interval: `{interval}` (corrected DQA2/DQB1/DRB3 window on ARS-UCD2.0 BTA23). "
        "Reads were extracted from the UCD2.0 target BAM (primary alignments, `-F 0x900`). "
        "Assembler: Flye 2.9.6 `--nano-raw --keep-haplotypes --min-overlap 3000 --genome-size 0.6m`. "
        "IPD-MHC panel (`ipd_mhc_bola_panel.fasta`) mapped with `minimap2 -x splice -p 0.2 -N 30 --secondary=yes`.",
        "",
        COPY_NUMBER_CAVEAT
        + " DRB3 assembly matches are compared to PCR-SBT but do not replace PCR-SBT typing.",
        "",
        f"Cohort: {len(summaries)} animals. *016:01 carriers: {len(carriers)}/{len(summaries)}.",
        "",
        "## Per-sample recovery",
        "",
        "| sample_id | PCR-SBT DRB3 | status | reads | DQA clusters | DQB clusters | DRB3 clusters | interpretation | DRB3 vs PCR-SBT | best DRB3 (identity) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['sample_id']} | {row['pcr_sbt']} | {row['status']} | {row['n_reads']} | "
            f"{row['n_dqa_clusters']} | {row['n_dqb_clusters']} | {row['n_drb3_clusters']} | "
            f"{row['interpretation']} | {row['drb3_concordance']} | {row['best_drb3_alleles'] or '-'} |"
        )

    lines.extend(
        [
            "",
            "## *016:01 carrier vs non-carrier (descriptive only)",
            "",
            f"- Carrier DQA clusters: {carrier_dqa}",
            f"- Carrier DQB clusters: {carrier_dqb}",
            f"- Non-carrier DQA clusters: {noncarrier_dqa}",
            f"- Non-carrier DQB clusters: {noncarrier_dqb}",
            "",
            "No consistent DQ cluster pattern distinguishes carriers from non-carriers in this "
            "ascertained 13-animal cohort. This is not a population association test.",
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
        "n_contigs",
        "n50",
        "n_dqa_clusters",
        "n_dqb_clusters",
        "n_drb3_clusters",
        "interpretation",
        "drb3_concordance",
        "best_drb3_alleles",
        "carrier_016_01",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)


def write_copies_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "sample_id",
        "family",
        "contig",
        "ref_start",
        "ref_end",
        "best_allele",
        "identity",
        "coverage",
        "mapq",
        "structure",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    samples = load_analytical_samples(args.registry)
    sample_ids = [row["sample_id"] for row in samples]
    by_id = {row["sample_id"]: row for row in samples}
    paf_name = paf_filename("dq_dr")

    summaries = []
    copy_rows = []
    hashes = {}
    for sample_id in sample_ids:
        sample_dir = args.root / sample_id
        status_path = sample_dir / f"{sample_id}.assembly_status.json"
        if not status_path.is_file():
            raise SystemExit(f"missing assembly status for {sample_id}: {status_path}")
        status = json.loads(status_path.read_text(encoding="utf-8"))
        paf_path = sample_dir / paf_name
        paf_text = paf_path.read_text(encoding="utf-8") if paf_path.is_file() else ""
        copies = cluster_alignments(parse_ipd_classii_paf(paf_text))
        registry_row = by_id[sample_id]
        summary = summarize_sample(
            sample_id=sample_id,
            copies=copies,
            status=status,
            hap1=registry_row.get("pcr_sbt_drb3_hap1", ""),
            hap2=registry_row.get("pcr_sbt_drb3_hap2", ""),
        )
        summaries.append(summary)
        hashes[sample_id] = sha256_file(status_path)
        for copy in copies:
            copy_rows.append(
                {
                    "sample_id": sample_id,
                    "family": copy.best.locus,
                    "contig": copy.contig,
                    "ref_start": copy.tstart,
                    "ref_end": copy.tend,
                    "best_allele": copy.best.allele,
                    "identity": round(copy.best.identity, 4),
                    "coverage": round(copy.best.coverage, 4),
                    "mapq": copy.best.mapq,
                    "structure": locus_structure(copies, copy.best.locus),
                }
            )

    validate_cohort_ids(sample_ids, [row["sample_id"] for row in summaries])
    args.root.mkdir(parents=True, exist_ok=True)
    write_recovery_tsv(args.root / "dq_dr_recovery.tsv", summaries)
    write_copies_tsv(args.root / "dq_dr_copies.tsv", copy_rows)
    payload = {
        "sample_count": len(summaries),
        "passed_count": sum(1 for row in summaries if row["status"] == "passed"),
        "interval": DQ_DR_REGION,
        "status": "passed" if all(row["status"] == "passed" for row in summaries) else "failed",
    }
    (args.root / "assembly_status.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(args.root, {"input_hashes": hashes})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(summaries, interval=DQ_DR_REGION), encoding="utf-8")
    print(f"DQ/DR assembly report written: {args.output} ({payload['passed_count']}/{payload['sample_count']} passed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
