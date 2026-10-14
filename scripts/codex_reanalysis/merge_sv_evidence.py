"""Classify Sniffles2 structural-variant calls and merge them across the cohort.

Only Sniffles2 is available as a caller in this workflow today. Per the
reanalysis design (`docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md`,
SV section): "A second caller or independent split-read/assembly evidence
will be used for manuscript-relevant SVs." A lone Sniffles2 call is therefore
``exploratory`` on its own; it is promoted to ``supported`` only by a second
caller, independent split-read support, or assembly support (Task 11) --
never by Sniffles2's own read-support count restated under a different name.

An event most of the cohort shares is judged a reference difference (this
individual diverging from the linear reference the same way most of the herd
does), not a private structural variant -- the same majority-sharing
heuristic already used for BTN/DQ region SVs in the legacy analysis
(`docs/HANDOFF_PLAN_PROGRESS.md` §2d: "多くは共通変異(≥8頭=参照差)").
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

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

SUPPORTED = "supported"
EXPLORATORY = "exploratory"
RECURRENT_REFERENCE_DIFFERENCE = "recurrent_reference_difference"

_BND_MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


class SvParseError(ValueError):
    """Raised when a Sniffles2 VCF record cannot be parsed."""


def classify_sv(
    callers: int,
    split_read_support: int,
    assembly_support: bool,
    minimum_split_read_support: int = 1,
) -> str:
    """A single caller's own read-support count never self-promotes a call.

    ``supported`` requires a second, independent source of evidence: a
    second caller (``callers >= 2``), independent split-read support at or
    above ``minimum_split_read_support``, or assembly support. Any one of
    those is sufficient; a lone Sniffles2 call with none of them stays
    ``exploratory`` regardless of its own SUPPORT/QUAL values.
    """
    independent_support = split_read_support >= minimum_split_read_support or assembly_support
    if callers >= 2 or independent_support:
        return SUPPORTED
    return EXPLORATORY


def is_recurrent_reference_difference(sample_count: int, cohort_size: int) -> bool:
    """True when more than half the cohort shares the same breakpoint event."""
    if cohort_size <= 0:
        raise ValueError("cohort_size must be positive")
    if sample_count < 0:
        raise ValueError("sample_count cannot be negative")
    return sample_count > cohort_size / 2


def _parse_info(info_field: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    for entry in info_field.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def _parse_bnd_mate(alt: str) -> tuple[str | None, int | None]:
    match = _BND_MATE_RE.search(alt)
    if match is None:
        return None, None
    return match.group(1), int(match.group(2))


def parse_sniffles_vcf(text: str) -> list[dict]:
    """Parse a Sniffles2 VCF body into one record per structural-variant call.

    Reads only the fields this project's classification and reporting use;
    verified against a real Sniffles2 2.8.0 run's header and PRECISE/
    IMPRECISE DEL/INS/INV/DUP/BND records (see the fixture in
    ``tests/codex_reanalysis/test_sv_evidence.py``), not invented against
    the spec alone.
    """
    records = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            raise SvParseError(f"VCF data line has fewer than 10 columns: {line!r}")
        chrom, pos, variant_id, _ref, alt, _qual, filter_field, info_field, format_keys, sample_value = fields[:10]
        info = _parse_info(info_field)
        sv_type = info.get("SVTYPE")
        if not isinstance(sv_type, str):
            raise SvParseError(f"VCF record is missing SVTYPE: {line!r}")
        key_index = {key: index for index, key in enumerate(format_keys.split(":"))}
        sample_fields = sample_value.split(":")

        def sample_value_for(key: str) -> str | None:
            index = key_index.get(key)
            if index is None or index >= len(sample_fields):
                return None
            return sample_fields[index]

        mate_chrom, mate_pos = (None, None)
        if sv_type == "BND":
            mate_chrom, mate_pos = _parse_bnd_mate(alt)

        records.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "id": variant_id,
                "sv_type": sv_type,
                "sv_len": int(info["SVLEN"]) if "SVLEN" in info else None,
                "end": int(info["END"]) if "END" in info else None,
                "support": int(info["SUPPORT"]) if "SUPPORT" in info else 0,
                "support_sa": int(info["SUPPORT_SA"]) if "SUPPORT_SA" in info else 0,
                "precise": "PRECISE" in info,
                "filter": filter_field,
                "chr2": info.get("CHR2", mate_chrom),
                "mate_pos": mate_pos,
                "genotype": sample_value_for("GT"),
                "ref_reads": int(sample_value_for("DR")) if sample_value_for("DR") not in (None, ".") else None,
                "alt_reads": int(sample_value_for("DV")) if sample_value_for("DV") not in (None, ".") else None,
            }
        )
    return records


def build_call_status(sample_id: str, vcf_path: Path, sniffles_version: str) -> dict:
    """Per-sample status for one Sniffles2 call, gated like every other stage output."""
    records = parse_sniffles_vcf(vcf_path.read_text(encoding="utf-8")) if vcf_path.is_file() else []
    return {
        "sample_id": sample_id,
        "sv_count": len(records),
        "sv_type_counts": dict(sorted(_count_by_key(records, "sv_type").items())),
        "sniffles_version": sniffles_version,
        "status": "passed" if vcf_path.is_file() else "failed",
    }


def _count_by_key(records: Sequence[dict], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for record in records:
        value = record[key]
        counts[value] = counts.get(value, 0) + 1
    return counts


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def group_by_breakpoint(
    calls: Sequence[dict], position_tolerance: int = 50
) -> list[dict]:
    """Cluster per-sample calls into one row per (chrom, sv_type, breakpoint bucket).

    Sniffles2 reports each sample's own estimate of a breakpoint's position;
    exact positions rarely agree across samples for the same real event, so
    calls are bucketed by rounding position down to the nearest
    ``position_tolerance``. This is a clustering heuristic, not an
    identity claim -- two real, distinct nearby events can collide into one
    bucket, which is why every row also carries its own min/max position
    span for review.
    """
    grouped: dict[tuple[str, str, int], list[dict]] = {}
    for call in calls:
        bucket_key = (call["chrom"], call["sv_type"], call["pos"] // position_tolerance)
        grouped.setdefault(bucket_key, []).append(call)

    rows = []
    for (chrom, sv_type, _bucket), members in grouped.items():
        sample_ids = sorted({member["sample_id"] for member in members})
        positions = [member["pos"] for member in members]
        rows.append(
            {
                "chrom": chrom,
                "sv_type": sv_type,
                "position_min": min(positions),
                "position_max": max(positions),
                "sample_count": len(sample_ids),
                "sample_ids": sample_ids,
                "total_support": sum(member["support"] for member in members),
                "total_support_sa": sum(member["support_sa"] for member in members),
                "any_imprecise": any(not member["precise"] for member in members),
                "calls": members,
            }
        )
    return sorted(rows, key=lambda row: (row["chrom"], row["position_min"], row["sv_type"]))


def build_sv_table(
    breakpoint_groups: Sequence[dict],
    cohort_size: int,
    callers_by_group: Sequence[int] | None = None,
    assembly_support_by_group: Sequence[bool] | None = None,
    independent_split_read_support_by_group: Sequence[int] | None = None,
    minimum_split_read_support: int = 1,
) -> list[dict]:
    """Attach a claim status to each breakpoint group.

    ``callers_by_group``, ``assembly_support_by_group`` and
    ``independent_split_read_support_by_group`` let a later stage supply
    evidence from a second SV caller, Task 11 assembly recovery, or a
    split-read count that is not Sniffles2's own ``SUPPORT_SA``. All three
    default to "not available" (one caller, no assembly, no independent
    split reads). Sniffles2 ``SUPPORT_SA`` is recorded but never used to
    self-promote a call.
    """
    callers_by_group = callers_by_group or [1] * len(breakpoint_groups)
    assembly_support_by_group = assembly_support_by_group or [False] * len(breakpoint_groups)
    independent_split_read_support_by_group = independent_split_read_support_by_group or [0] * len(
        breakpoint_groups
    )
    if not (
        len(callers_by_group)
        == len(assembly_support_by_group)
        == len(independent_split_read_support_by_group)
        == len(breakpoint_groups)
    ):
        raise ValueError("per-group evidence lists must match the number of breakpoint groups")

    rows = []
    for group, callers, assembly_support, independent_split in zip(
        breakpoint_groups,
        callers_by_group,
        assembly_support_by_group,
        independent_split_read_support_by_group,
    ):
        recurrent = is_recurrent_reference_difference(group["sample_count"], cohort_size)
        status = (
            RECURRENT_REFERENCE_DIFFERENCE
            if recurrent
            else classify_sv(
                callers=callers,
                split_read_support=independent_split,
                assembly_support=assembly_support,
                minimum_split_read_support=minimum_split_read_support,
            )
        )
        rows.append(
            {
                "chrom": group["chrom"],
                "sv_type": group["sv_type"],
                "position_min": group["position_min"],
                "position_max": group["position_max"],
                "sample_count": group["sample_count"],
                "cohort_size": cohort_size,
                "callers": callers,
                "split_read_support": independent_split,
                "sniffles_support_sa": group["total_support_sa"],
                "assembly_support": assembly_support,
                "status": status,
            }
        )
    return rows


def write_sv_table_tsv(path: Path, rows: Sequence[dict]) -> None:
    fieldnames = (
        "chrom",
        "sv_type",
        "position_min",
        "position_max",
        "sample_count",
        "cohort_size",
        "callers",
        "split_read_support",
        "sniffles_support_sa",
        "assembly_support",
        "status",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def render_report(rows: Sequence[dict], sample_ids: Sequence[str]) -> str:
    supported = [row for row in rows if row["status"] == SUPPORTED]
    exploratory = [row for row in rows if row["status"] == EXPLORATORY]
    recurrent = [row for row in rows if row["status"] == RECURRENT_REFERENCE_DIFFERENCE]
    lines = [
        "# ARS-UCD2.0 Structural Variant Evidence Report",
        "",
        f"Cohort: {len(sample_ids)} animals. {len(rows)} breakpoint groups from Sniffles2.",
        "",
        "Callset scope: Sniffles2 2.8.0 on the target BAM, restricted to `NC_037349.1` "
        "(CNTN3) and `NC_037350.1` (BoLA/MHC). Control-gene islands were left out of "
        "scope after an unrestricted run hit a Sniffles2 coverage bug "
        "(`UnboundLocalError` on `end`) on ALB/PRNP/RPP30 contigs in 8/13 samples. "
        "There is no second SV caller yet, so `supported` is expected to be 0; "
        "recurrent events are reference-difference labels, not private SVs.",
        "",
        "**Only `supported` events are eligible for locus-specific manuscript statements.** "
        "`supported` requires a second caller, independent split-read support, or assembly "
        "support beyond Sniffles2's own call; `exploratory` is a single-caller call with none "
        "of those. `recurrent_reference_difference` events are shared by more than half the "
        "cohort and are far more likely to be this reference's divergence from the herd than a "
        "real private structural variant, so they are reported separately rather than as "
        "individual findings.",
        "",
        f"- supported: {len(supported)}",
        f"- exploratory: {len(exploratory)}",
        f"- recurrent_reference_difference: {len(recurrent)}",
        "",
        "## Supported events",
        "",
        "| chrom | sv_type | position_min | position_max | sample_count | callers | split_read_support | assembly_support |",
        "|---|---|---|---|---|---|---|---|",
    ]
    if supported:
        for row in supported:
            lines.append(
                f"| {row['chrom']} | {row['sv_type']} | {row['position_min']} | {row['position_max']} | "
                f"{row['sample_count']}/{row['cohort_size']} | {row['callers']} | {row['split_read_support']} | "
                f"{row['assembly_support']} |"
            )
    else:
        lines.append("| - | - | - | - | - | - | - | none |")
    lines.append("")
    return "\n".join(lines)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("record-call", "merge"), default="merge")
    parser.add_argument("--root", type=Path, help="05_sv directory containing per-sample VCFs (merge mode)")
    parser.add_argument("--registry", type=Path, help="00_input_audit/sample_registry.tsv (merge mode)")
    parser.add_argument("--position-tolerance", type=int, default=50)
    parser.add_argument("--minimum-split-read-support", type=int, default=1)
    parser.add_argument("--output", type=Path, help="Markdown report destination (merge mode)")
    parser.add_argument("--sample-id", help="record-call mode")
    parser.add_argument("--vcf", type=Path, help="record-call mode")
    parser.add_argument("--sniffles-version", help="record-call mode")
    parser.add_argument("--output-dir", type=Path, help="record-call mode")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)

    if args.mode == "record-call":
        for name, value in (
            ("--sample-id", args.sample_id),
            ("--vcf", args.vcf),
            ("--sniffles-version", args.sniffles_version),
            ("--output-dir", args.output_dir),
        ):
            if value is None:
                raise SystemExit(f"record-call mode requires {name}")
        status = build_call_status(args.sample_id, args.vcf, args.sniffles_version)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / f"{args.sample_id}.sniffles_status.json").write_text(
            json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        write_provenance(
            args.output_dir,
            {"input_hashes": {"vcf": sha256_file(args.vcf)} if args.vcf.is_file() else {}},
        )
        print(f"Sniffles2 call recorded for {args.sample_id}: {status['sv_count']} SVs -> {args.output_dir}")
        return 0

    for name, value in (("--root", args.root), ("--registry", args.registry), ("--output", args.output)):
        if value is None:
            raise SystemExit(f"merge mode requires {name}")

    samples = read_tsv(args.registry)
    sample_ids = [row["sample_id"] for row in samples if row.get("analysis_status", "").strip() == "analytical"]

    calls = []
    for sample_id in sample_ids:
        vcf_path = args.root / sample_id / f"{sample_id}.sniffles2.vcf"
        if not vcf_path.is_file():
            continue
        for record in parse_sniffles_vcf(vcf_path.read_text(encoding="utf-8")):
            record["sample_id"] = sample_id
            calls.append(record)

    groups = group_by_breakpoint(calls, position_tolerance=args.position_tolerance)
    rows = build_sv_table(
        groups, cohort_size=len(sample_ids), minimum_split_read_support=args.minimum_split_read_support
    )

    args.root.mkdir(parents=True, exist_ok=True)
    write_sv_table_tsv(args.root / "sv_evidence_table.tsv", rows)

    status_payload = {
        "sample_count": len(sample_ids),
        "breakpoint_group_count": len(rows),
        "supported_count": sum(1 for row in rows if row["status"] == SUPPORTED),
        "exploratory_count": sum(1 for row in rows if row["status"] == EXPLORATORY),
        "recurrent_reference_difference_count": sum(
            1 for row in rows if row["status"] == RECURRENT_REFERENCE_DIFFERENCE
        ),
        "status": "passed",
    }
    (args.root / "sv_evidence_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        args.root,
        {
            "input_hashes": {
                sample_id: sha256_file(args.root / sample_id / f"{sample_id}.sniffles2.vcf")
                for sample_id in sample_ids
                if (args.root / sample_id / f"{sample_id}.sniffles2.vcf").is_file()
            }
        },
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(rows, sample_ids), encoding="utf-8")
    print(f"SV evidence report written: {args.output} ({status_payload['supported_count']}/{len(rows)} supported)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
