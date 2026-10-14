"""Resolve published marker coordinates onto ARS-UCD2.0 without silent reuse.

The three published GWAS SNPs (rs110616206, rs29026690, rs17872126) and the
published TNF -824A/G promoter polymorphism were all originally located on
ARS-UCD1.2 (see ``docs/CNTN3_BLV_LOCUS.md`` and
``docs/CLASSIII_016_TAG_SNP_EXPLORATION.md`` #6). Per the reanalysis design
(``docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md`` §5.2),
a UCD1.2/UCD2.0 contig difference is a scientific decision point, not a
liftover to perform silently. This script makes that decision explicit: it
verifies -- by whole-chromosome SHA-256, not by assuming -- whether the
BTA22/BTA23 sequence is byte-identical between the two assemblies for the
chromosome each marker sits on, and only then reuses the legacy position.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file


class MarkerCoordinateError(ValueError):
    """Raised when a marker's chromosome cannot be verified identical across assemblies."""


def read_marker_table(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_faidx_command(fasta_path: str, contig: str, samtools: str = "samtools") -> list[str]:
    return [samtools, "faidx", str(fasta_path), contig]


def parse_fasta_sequence(text: str) -> str:
    """Concatenate a single-record FASTA body into one uppercase sequence."""
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise MarkerCoordinateError("faidx output did not contain a FASTA record")
    return "".join(lines[1:]).upper()


def parse_faidx_base(text: str) -> str:
    """Read the single-base FASTA that ``samtools faidx contig:pos-pos`` emits."""
    sequence = parse_fasta_sequence(text)
    if len(sequence) != 1:
        raise MarkerCoordinateError(f"expected one base from faidx, got {len(sequence)}")
    return sequence


def annotate_reference_alleles(
    resolved: Sequence[dict], fetch_base
) -> list[dict]:
    """Attach the UCD2.0 reference base. ``fetch_base(contig, position) -> str``."""
    annotated = []
    for row in resolved:
        payload = dict(row)
        if row.get("resolution_status") == "verified_identical_chromosome" and row.get("ucd2_0_position"):
            payload["ucd2_0_reference_allele"] = fetch_base(row["ucd2_0_contig"], int(row["ucd2_0_position"]))
        else:
            payload["ucd2_0_reference_allele"] = ""
        annotated.append(payload)
    return annotated


def hash_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def decide_resolution(ucd1_2_sha256: str, ucd2_0_sha256: str) -> str:
    """A UCD1.2/UCD2.0 contig mismatch is a decision point, never a silent reuse."""
    return "verified_identical_chromosome" if ucd1_2_sha256 == ucd2_0_sha256 else "decision_required"


def identity_from_comparison_rows(rows: Sequence[dict[str, str]]) -> dict[str, dict]:
    """Build the contig-pair identity map from ``sequence_comparison.tsv``.

    Hash equality is the decision, not the ``status`` column: a row labelled
    identical with mismatched hashes is still ``decision_required``.
    """
    checked: dict[str, dict] = {}
    for row in rows:
        identity = {
            "ucd1_2_contig": row["old_fasta_name"],
            "ucd2_0_contig": row["new_accession"],
            "ucd1_2_sha256": row["old_sha256"],
            "ucd2_0_sha256": row["new_sha256"],
            "identical": row["old_sha256"] == row["new_sha256"],
        }
        checked[f"{row['old_fasta_name']}->{row['new_accession']}"] = identity
        checked[f"{row['old_accession']}->{row['new_accession']}"] = {
            **identity,
            "ucd1_2_contig": row["old_accession"],
        }
    return checked


def sha256_of_contig(fasta_path: Path, contig: str, samtools: str = "samtools") -> str:
    """Hash one contig's uppercase sequence, independent of line wrapping or case."""
    result = subprocess.run(
        build_faidx_command(str(fasta_path), contig, samtools),
        capture_output=True,
        check=True,
        text=True,
    )
    return hash_sequence(parse_fasta_sequence(result.stdout))


def fetch_reference_base(fasta_path: Path, contig: str, position: int, samtools: str = "samtools") -> str:
    """Return the single UCD2.0 base at a 1-based closed position."""
    region = f"{contig}:{position}-{position}"
    result = subprocess.run(
        build_faidx_command(str(fasta_path), region, samtools),
        capture_output=True,
        check=True,
        text=True,
    )
    return parse_faidx_base(result.stdout)


def verify_chromosome_identity(
    ucd1_2_fasta: Path,
    ucd2_0_fasta: Path,
    ucd1_2_contig: str,
    ucd2_0_contig: str,
    samtools: str = "samtools",
) -> dict:
    ucd1_2_sha256 = sha256_of_contig(ucd1_2_fasta, ucd1_2_contig, samtools)
    ucd2_0_sha256 = sha256_of_contig(ucd2_0_fasta, ucd2_0_contig, samtools)
    return {
        "ucd1_2_contig": ucd1_2_contig,
        "ucd2_0_contig": ucd2_0_contig,
        "ucd1_2_sha256": ucd1_2_sha256,
        "ucd2_0_sha256": ucd2_0_sha256,
        "identical": ucd1_2_sha256 == ucd2_0_sha256,
    }


def resolve_markers_from_identity(
    markers: Sequence[dict[str, str]], checked: dict[str, dict]
) -> list[dict]:
    """Build resolved marker rows from precomputed per-chromosome identity checks.

    Pure and independent of subprocess/filesystem access, so the resolution
    logic (one row per marker, decision-required on any mismatch) can be
    tested without genomic data.
    """
    resolved: list[dict] = []
    for marker in markers:
        ucd1_2_contig = marker["ucd1_2_contig"]
        ucd2_0_contig = marker["ucd2_0_contig"]
        key = f"{ucd1_2_contig}->{ucd2_0_contig}"
        if key not in checked:
            raise MarkerCoordinateError(f"no chromosome comparison for {key}")
        identity = checked[key]
        status = decide_resolution(identity["ucd1_2_sha256"], identity["ucd2_0_sha256"])
        resolved.append(
            {
                "marker_id": marker["marker_id"],
                "ucd1_2_contig": ucd1_2_contig,
                "ucd1_2_position": int(marker["ucd1_2_position"]),
                "ucd2_0_contig": ucd2_0_contig,
                "ucd2_0_position": int(marker["ucd1_2_position"]) if status == "verified_identical_chromosome" else None,
                "resolution_status": status,
                "citation": marker.get("citation", ""),
                "notes": marker.get("notes", ""),
            }
        )
    return resolved


def resolve_markers(
    markers: Sequence[dict[str, str]],
    ucd1_2_fasta: Path,
    ucd2_0_fasta: Path,
    samtools: str = "samtools",
) -> tuple[list[dict], dict[str, dict]]:
    """Resolve every marker's UCD2.0 coordinate, verifying each contig once."""
    checked: dict[str, dict] = {}
    for marker in markers:
        key = f"{marker['ucd1_2_contig']}->{marker['ucd2_0_contig']}"
        if key not in checked:
            checked[key] = verify_chromosome_identity(
                ucd1_2_fasta, ucd2_0_fasta, marker["ucd1_2_contig"], marker["ucd2_0_contig"], samtools
            )
    return resolve_markers_from_identity(markers, checked), checked


def write_resolved_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "marker_id",
        "ucd1_2_contig",
        "ucd1_2_position",
        "ucd2_0_contig",
        "ucd2_0_position",
        "resolution_status",
        "ucd2_0_reference_allele",
        "citation",
        "notes",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", required=True, type=Path)
    parser.add_argument("--assembly-comparison", type=Path, help="01_reference/assembly_comparison/sequence_comparison.tsv")
    parser.add_argument("--ucd1-2-fasta", type=Path)
    parser.add_argument("--ucd2-0-fasta", type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--output", required=True, type=Path, help="04_small_variants/coordinates directory")
    args = parser.parse_args(argv)

    markers = read_marker_table(args.markers)
    if args.assembly_comparison is not None:
        comparison_rows = read_marker_table(args.assembly_comparison)
        checked = identity_from_comparison_rows(comparison_rows)
        resolved = resolve_markers_from_identity(markers, checked)
        input_hashes = {"assembly_comparison": sha256_file(args.assembly_comparison)}
    elif args.ucd1_2_fasta is not None and args.ucd2_0_fasta is not None:
        resolved, checked = resolve_markers(markers, args.ucd1_2_fasta, args.ucd2_0_fasta, args.samtools)
        input_hashes = {
            key: {"ucd1_2_sha256": value["ucd1_2_sha256"], "ucd2_0_sha256": value["ucd2_0_sha256"]}
            for key, value in checked.items()
        }
    else:
        raise SystemExit("provide --assembly-comparison or both --ucd1-2-fasta and --ucd2-0-fasta")

    if args.ucd2_0_fasta is not None:
        fasta = args.ucd2_0_fasta
        samtools = args.samtools
        resolved = annotate_reference_alleles(
            resolved,
            lambda contig, position, fasta=fasta, samtools=samtools: fetch_reference_base(
                fasta, contig, position, samtools
            ),
        )

    args.output.mkdir(parents=True, exist_ok=True)
    write_resolved_tsv(args.output / "marker_coordinate_registry.tsv", resolved)

    unresolved = [row["marker_id"] for row in resolved if row["resolution_status"] == "decision_required"]
    status_payload = {
        "marker_count": len(resolved),
        "chromosome_checks": checked,
        "unresolved_markers": unresolved,
        "status": "passed" if not unresolved else "decision_required",
    }
    (args.output / "marker_coordinate_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        args.output,
        {
            "input_hashes": input_hashes,
            "reference_hashes": {
                key: {"ucd1_2_sha256": value["ucd1_2_sha256"], "ucd2_0_sha256": value["ucd2_0_sha256"]}
                for key, value in checked.items()
            },
        },
    )
    print(f"marker coordinates resolved: {len(resolved)} markers, {len(unresolved)} decision_required")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
