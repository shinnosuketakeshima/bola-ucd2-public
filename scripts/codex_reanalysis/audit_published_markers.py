"""Audit published BLV-associated marker genotypes on ARS-UCD2.0.

Covers the three published GWAS SNPs (rs110616206, rs29026690, rs17872126)
and the published TNF -824A/G promoter polymorphism (Konnai et al. 2006),
using the coordinates ``resolve_marker_coordinates.py`` verified onto
ARS-UCD2.0. No Clair3 model verified against this cohort's Guppy5/R9.4.1
chemistry is available on this host, so every marker's caller status is
``exploratory`` per
``docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md`` §5.5:
a genotype is reported only alongside its two independent-pileup depths, and
"accepted" here means the two pileups agree and depth clears the configured
minimum -- it does not mean caller-verified.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.call_small_variants import genotype_marker
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

FOCAL_ALLELE = "*016:01"
DEFAULT_MINIMUM_DEPTH = 10
SHA256_HEX = re.compile(r"^[0-9a-f]{64}$")
UNRESOLVED = "coordinate_unresolved"
LEGACY_TABLE5_COLUMNS = {
    "rs110616206": "rs110616206_CNTN3",
    "rs29026690": "rs29026690_ABHD16A",
    "rs17872126": "rs17872126_PRR3",
}


def classify_marker(
    primary: str,
    pileup: str,
    depth: int,
    minimum_depth: int = DEFAULT_MINIMUM_DEPTH,
    pileup_depth: int | None = None,
) -> str:
    """Gate a genotype call on cross-pileup agreement and minimum depth.

    Never returns ``accepted`` on a caller/model result -- there is none
    verified for this chemistry here. It only reports that the two
    independent pileups agree at sufficient depth; the report text makes the
    ``exploratory`` caller status explicit regardless of this outcome.

    Both pileups must clear ``minimum_depth``. Passing only ``depth`` keeps
    the historical single-depth tests working by applying that value to both.
    """
    if pileup_depth is None:
        pileup_depth = depth
    if primary == "no_data" or pileup == "no_data":
        return "no_data"
    if primary != pileup:
        return "conflicting"
    if min(depth, pileup_depth) < minimum_depth:
        return "limited"
    return "accepted"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_analytical_samples(registry_path: Path) -> list[dict[str, str]]:
    rows = read_tsv(registry_path)
    return [row for row in rows if row.get("analysis_status", "").strip() == "analytical"]


def is_016_carrier(sample_row: Sequence[str]) -> bool:
    return FOCAL_ALLELE in sample_row


def as_carrier(value: object) -> bool:
    """Coerce a TSV/JSON carrier flag. The string ``False`` is not a carrier."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "yes", "1"}


def normalize_diplotype(genotype: str) -> str:
    """Sort the two alleles so ``T/C`` and ``C/T`` compare equal."""
    if not genotype or "/" not in genotype:
        return genotype
    left, right = genotype.split("/", 1)
    return "/".join(sorted([left, right]))


def legacy_table5_concordance(genotype_rows: Sequence[dict], legacy_rows: Sequence[dict]) -> dict:
    """Compare UCD2.0 pileup alleles with the ARS-UCD1.2 draft Table 5.

    Limited-depth rows are included in the allele check; they are not
    ``accepted`` evidence and must not be used to recompute r².
    """
    legacy_by_sample = {row["sample_id"]: row for row in legacy_rows}
    per_marker: dict[str, dict[str, int]] = {
        marker_id: {
            "n": 0,
            "n_match": 0,
            "n_mismatch": 0,
            "n_limited_matching": 0,
            "n_accepted_matching": 0,
        }
        for marker_id in LEGACY_TABLE5_COLUMNS
    }
    n_compared = 0
    n_mismatch = 0
    n_limited_matching = 0
    for row in genotype_rows:
        column = LEGACY_TABLE5_COLUMNS.get(row["marker_id"])
        if column is None:
            continue
        expected = normalize_diplotype(legacy_by_sample[row["sample_id"]][column])
        observed = normalize_diplotype(row["primary_genotype"])
        stats = per_marker[row["marker_id"]]
        stats["n"] += 1
        n_compared += 1
        if expected != observed:
            stats["n_mismatch"] += 1
            n_mismatch += 1
            continue
        stats["n_match"] += 1
        if row["evidence_status"] == "limited":
            stats["n_limited_matching"] += 1
            n_limited_matching += 1
        elif row["evidence_status"] == "accepted":
            stats["n_accepted_matching"] += 1
    return {
        "n_compared": n_compared,
        "n_mismatch": n_mismatch,
        "n_limited_matching": n_limited_matching,
        "per_marker": per_marker,
    }


def mapping_bam_hashes(mapping_root: Path, sample_ids: Sequence[str]) -> dict[str, str]:
    """Read BAM SHA-256 digests recorded by mapping; never treat a path as a hash."""
    hashes: dict[str, str] = {}
    for sample_id in sample_ids:
        status_path = Path(mapping_root) / f"{sample_id}.mapping_status.json"
        payload = json.loads(status_path.read_text(encoding="utf-8"))
        digest = payload.get("artifact_hashes", {}).get("target_bam", "")
        if not SHA256_HEX.fullmatch(str(digest)):
            raise ValueError(f"{status_path}: target_bam is not a SHA-256 digest")
        hashes[sample_id] = digest
    return hashes


def require_complete_marker_matrix(
    rows: Sequence[dict], sample_ids: Sequence[str], marker_ids: Sequence[str]
) -> None:
    """Fail unless every sample x marker pair appears exactly once."""
    expected = {(marker_id, sample_id) for marker_id in marker_ids for sample_id in sample_ids}
    observed = [(row["marker_id"], row["sample_id"]) for row in rows]
    if len(observed) != len(set(observed)):
        raise ValueError("duplicate sample x marker rows")
    missing = expected - set(observed)
    extra = set(observed) - expected
    if missing or extra:
        raise ValueError(f"incomplete marker matrix: missing {sorted(missing)} extra {sorted(extra)}")


def _empty_genotype() -> dict:
    return {
        "primary_genotype": "no_data",
        "primary_depth": 0,
        "primary_ref_count": 0,
        "primary_alt_count": 0,
        "pileup_genotype": "no_data",
        "pileup_depth": 0,
        "pileup_ref_count": 0,
        "pileup_alt_count": 0,
    }


def row_from_genotype(
    marker: dict[str, str],
    sample: dict[str, str],
    result: dict | None,
    minimum_depth: int,
) -> dict:
    """Build one audit row. ``result is None`` means the coordinate was not resolved."""
    carrier = is_016_carrier((sample.get("pcr_sbt_drb3_hap1", ""), sample.get("pcr_sbt_drb3_hap2", "")))
    if marker.get("resolution_status") != "verified_identical_chromosome" or result is None:
        payload = _empty_genotype()
        evidence_status = UNRESOLVED
        position: int | str = ""
    else:
        payload = {
            "primary_genotype": result["primary_genotype"],
            "primary_depth": result["primary_depth"],
            "primary_ref_count": result["primary_ref_count"],
            "primary_alt_count": result["primary_alt_count"],
            "pileup_genotype": result["pileup_genotype"],
            "pileup_depth": result["pileup_depth"],
            "pileup_ref_count": result["pileup_ref_count"],
            "pileup_alt_count": result["pileup_alt_count"],
        }
        evidence_status = classify_marker(
            payload["primary_genotype"],
            payload["pileup_genotype"],
            payload["primary_depth"],
            minimum_depth,
            pileup_depth=payload["pileup_depth"],
        )
        position = int(marker["ucd2_0_position"])
    return {
        "marker_id": marker["marker_id"],
        "sample_id": sample["sample_id"],
        "contig": marker.get("ucd2_0_contig", ""),
        "position": position,
        **payload,
        "caller_status": "exploratory",
        "evidence_status": evidence_status,
        "carrier_016_01": carrier,
    }


def build_marker_genotype_rows(
    markers: Sequence[dict[str, str]],
    samples: Sequence[dict[str, str]],
    mapping_root: Path,
    reference_fasta: Path,
    bcftools: str,
    samtools: str,
    mapq: int,
    baseq: int,
    minimum_depth: int,
) -> list[dict]:
    rows = []
    for marker in markers:
        if marker["resolution_status"] != "verified_identical_chromosome":
            for sample in samples:
                rows.append(row_from_genotype(marker, sample, None, minimum_depth))
            continue

        contig = marker["ucd2_0_contig"]
        position = int(marker["ucd2_0_position"])
        for sample in samples:
            bam = Path(mapping_root) / f"{sample['sample_id']}.targets.bam"
            result = genotype_marker(
                str(bam), str(reference_fasta), contig, position, bcftools, samtools, mapq, baseq
            )
            rows.append(row_from_genotype(marker, sample, result, minimum_depth))
    sample_ids = [sample["sample_id"] for sample in samples]
    marker_ids = [marker["marker_id"] for marker in markers]
    require_complete_marker_matrix(rows, sample_ids, marker_ids)
    return rows


def write_marker_genotypes_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "marker_id",
        "sample_id",
        "contig",
        "position",
        "primary_genotype",
        "primary_depth",
        "primary_ref_count",
        "primary_alt_count",
        "pileup_genotype",
        "pileup_depth",
        "pileup_ref_count",
        "pileup_alt_count",
        "caller_status",
        "evidence_status",
        "carrier_016_01",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_report(
    rows: Sequence[dict],
    markers: Sequence[dict],
    legacy_rows: Sequence[dict] | None = None,
) -> str:
    lines = [
        "# ARS-UCD2.0 Published-Marker Genotype Audit",
        "",
        "**Caller status: exploratory for every marker.** No Clair3 model verified against this "
        "cohort's Guppy5/R9.4.1 chemistry is available on this host "
        "(`docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md` §5.5), so no "
        "genotype below is a caller/model result. Each is a BCFtools allele-depth pileup (`primary`) "
        "cross-checked against an independently coded raw-pileup base counter (`pileup`); "
        "`evidence_status=accepted` means those two pileups agree at sufficient depth, not that a "
        "verified caller confirmed it. `carrier_016_01` cross-tabs are exploratory hypothesis-generating "
        "counts in a small, *016:01-ascertained 13-animal cohort, not a population association test.",
        "",
        "AD values are **ref/alt counts against the linear-reference base**, not against the first "
        "allele of the genotype label. A `C/C (0/6)` row at a T-reference marker is homozygous "
        "alternate C, not homozygous reference.",
        "",
    ]
    by_marker: dict[str, list[dict]] = {}
    for row in rows:
        by_marker.setdefault(row["marker_id"], []).append(row)

    marker_meta = {marker["marker_id"]: marker for marker in markers}
    for marker_id, marker_rows in by_marker.items():
        meta = marker_meta.get(marker_id, {})
        lines.append(f"## {marker_id}")
        lines.append("")
        lines.append(f"- Coordinate: `{meta.get('ucd2_0_contig', '?')}:{meta.get('ucd2_0_position', '?')}` "
                      f"(resolution: {meta.get('resolution_status', '?')})")
        ref_allele = meta.get("ucd2_0_reference_allele", "")
        if ref_allele:
            lines.append(f"- Reference allele {ref_allele} (ARS-UCD2.0 + strand)")
        lines.append(f"- Citation: {meta.get('citation', '')}")
        if meta.get("notes"):
            lines.append(f"- Notes: {meta['notes']}")
        accepted = sum(1 for row in marker_rows if row["evidence_status"] == "accepted")
        conflicting = sum(1 for row in marker_rows if row["evidence_status"] == "conflicting")
        lines.append(
            f"- {accepted}/{len(marker_rows)} samples: two-pileup agreement at sufficient depth "
            f"(`accepted`); {conflicting} `conflicting`."
        )
        lines.append("")
        lines.append("| sample_id | primary_genotype (AD) | primary_depth | pileup_genotype (manual) | "
                      "pileup_depth | evidence_status | *016:01 carrier |")
        lines.append("|---|---|---|---|---|---|---|")
        for row in sorted(marker_rows, key=lambda item: item["sample_id"]):
            primary_ad = f"{row.get('primary_ref_count', '')}/{row.get('primary_alt_count', '')}"
            lines.append(
                f"| {row['sample_id']} | {row['primary_genotype']} ({primary_ad}) | "
                f"{row.get('primary_depth', '')} | "
                f"{row.get('pileup_genotype', '')} | {row.get('pileup_depth', '')} | "
                f"{row['evidence_status']} | "
                f"{'yes' if as_carrier(row.get('carrier_016_01', False)) else 'no'} |"
            )
        lines.append("")

        genotype_by_carrier: dict[bool, dict[str, int]] = {True: {}, False: {}}
        for row in marker_rows:
            if row["evidence_status"] != "accepted":
                continue
            bucket = genotype_by_carrier[as_carrier(row.get("carrier_016_01", False))]
            bucket[row["primary_genotype"]] = bucket.get(row["primary_genotype"], 0) + 1
        lines.append(
            f"*016:01 carriers (accepted calls only): {genotype_by_carrier[True]}. "
            f"Non-carriers: {genotype_by_carrier[False]}."
        )
        lines.append("")
    if legacy_rows is not None:
        summary = legacy_table5_concordance(rows, legacy_rows)
        lines.extend(
            [
                "## Concordance with legacy Table 5",
                "",
                "Alleles are compared after sorting (`T/C` = `C/T`) against "
                "`docs/tables/table5_gwas_snp_genotypes.tsv`. Limited-depth rows are included "
                "in this allele check; they are not `accepted` evidence. The manuscript "
                "r² = 0.74 between *016:01 and rs29026690-T is a legacy Table 5 calculation, "
                "not a claim re-established from accepted UCD2.0 pileups.",
                "",
                f"Allele labels match in {summary['n_compared'] - summary['n_mismatch']}/{summary['n_compared']} "
                f"sample × SNP cells"
                + ("" if summary["n_mismatch"] == 0 else f" ({summary['n_mismatch']} mismatches)")
                + ".",
                "",
                "| marker | n | allele_match | mismatch | limited_and_matching | accepted_and_matching |",
                "|---|---|---|---|---|---|",
            ]
        )
        for marker_id, stats in summary["per_marker"].items():
            if stats["n"] == 0:
                continue
            lines.append(
                f"| {marker_id} | {stats['n']} | {stats['n_match']} | {stats['n_mismatch']} | "
                f"{stats['n_limited_matching']} | {stats['n_accepted_matching']} |"
            )
        lines.append("")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", required=True, type=Path, help="marker_coordinate_registry.tsv")
    parser.add_argument("--registry", required=True, type=Path, help="00_input_audit/sample_registry.tsv")
    parser.add_argument("--mapping-root", required=True, type=Path, help="02_mapping directory")
    parser.add_argument("--reference", required=True, type=Path, help="ARS-UCD2.0 FASTA")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--baseq", type=int, default=20)
    parser.add_argument("--minimum-depth", type=int, default=DEFAULT_MINIMUM_DEPTH)
    parser.add_argument("--reference-seal", type=Path, default=None, help="01_reference/mapping_reference_seal.json")
    parser.add_argument("--output", required=True, type=Path, help="04_small_variants stage directory")
    parser.add_argument("--report", required=True, type=Path, help="Markdown report destination")
    parser.add_argument(
        "--legacy-table5",
        type=Path,
        default=Path("docs/tables/table5_gwas_snp_genotypes.tsv"),
        help="legacy ARS-UCD1.2 Table 5 genotypes for allele concordance",
    )
    args = parser.parse_args(argv)

    markers = read_tsv(args.markers)
    samples = load_analytical_samples(args.registry)
    if len(samples) != 13:
        raise SystemExit(f"expected 13 analytical samples, got {len(samples)}")
    rows = build_marker_genotype_rows(
        markers,
        samples,
        args.mapping_root,
        args.reference,
        args.bcftools,
        args.samtools,
        args.mapq,
        args.baseq,
        args.minimum_depth,
    )

    args.output.mkdir(parents=True, exist_ok=True)
    write_marker_genotypes_tsv(args.output / "marker_genotypes.tsv", rows)

    status_payload = {
        "marker_count": len(markers),
        "sample_count": len(samples),
        "row_count": len(rows),
        "accepted_count": sum(1 for row in rows if row["evidence_status"] == "accepted"),
        "conflicting_count": sum(1 for row in rows if row["evidence_status"] == "conflicting"),
        "limited_count": sum(1 for row in rows if row["evidence_status"] == "limited"),
        "unresolved_count": sum(1 for row in rows if row["evidence_status"] == UNRESOLVED),
        "caller_status": "exploratory",
        "status": "passed" if not any(row["evidence_status"] == UNRESOLVED for row in rows) else "decision_required",
    }
    (args.output / "marker_audit_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    sample_ids = [sample["sample_id"] for sample in samples]
    reference_hashes = {}
    if args.reference_seal is not None:
        seal = json.loads(args.reference_seal.read_text(encoding="utf-8"))
        fasta_name = Path(args.reference).name
        digest = seal["resources"][fasta_name]["sha256"]
        if not SHA256_HEX.fullmatch(digest):
            raise SystemExit(f"{args.reference_seal}: fasta hash is not SHA-256")
        reference_hashes["fasta"] = digest
    write_provenance(
        args.output,
        {
            "input_hashes": {
                "markers": sha256_file(args.markers),
                "sample_registry": sha256_file(args.registry),
                "bams": mapping_bam_hashes(args.mapping_root, sample_ids),
            },
            "reference_hashes": reference_hashes,
        },
    )

    args.report.parent.mkdir(parents=True, exist_ok=True)
    legacy_rows = read_tsv(args.legacy_table5) if args.legacy_table5.is_file() else None
    args.report.write_text(render_report(rows, markers, legacy_rows=legacy_rows), encoding="utf-8")
    print(f"marker audit written: {args.report} ({status_payload['accepted_count']}/{len(rows)} accepted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
