"""WGS-based DRB3 typing vs PCR-SBT on ARS-UCD2.0 target BAMs.

Ports the read-based half of ``docs/DRB3_WGS_CONCORDANCE.md``
(``results/annotation/drb3_wgs_readbased_concordance.tsv``) onto UCD2.0
target BAMs, and summarises the already-computed assembly-based
concordance from ``06_local_assembly/dq_dr/dq_dr_recovery.tsv``.

Method (read-based):
1. Extract primary alignments from the real DRB3 window on each target BAM.
2. Map those reads to the IPD-MHC BoLA-DRB3 allele panel.
3. Assign each read its best-identity allele; tally allele counts.
4. Compare top alleles and known-allele presence against PCR-SBT.

The manuscript claim is deliberately self-limiting: WGS alone does not
reliably reproduce PCR-SBT DRB3 types (central-allele mapping bias toward
*016:01 / *027:03). This script reproduces that negative result.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.dq_dr_recovery import normalize_drb3_allele
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

# Legacy concordance used NC_037350.1:25,720,000-25,738,000 (real coords).
# Gene model (0-based half-open) is 25,723,716-25,734,819; keep the same window.
DRB3_REGION = "NC_037350.1:25720000-25738000"
DEFAULT_MIN_ALN = 200
DEFAULT_MIN_IDENTITY = 0.70
EXTRACT_FLAGS = "0x900"
CENTRAL_BIAS_ALLELES = frozenset({"*016:01", "*027:03"})


class Drb3WgsTypingError(ValueError):
    """Raised when inputs cannot be analysed safely."""


@dataclass(frozen=True)
class AlleleHit:
    read: str
    allele: str
    identity: float
    nmatch: int
    aln_len: int


def parse_drb3_header(header: str) -> str | None:
    """Return normalized ``*016:01``-style allele from an IPD FASTA header, else None."""
    token = header[1:].split()[0] if header.startswith(">") else header.split()[0]
    if not token.startswith("BoLA-DRB3*"):
        return None
    try:
        return normalize_drb3_allele(token)
    except ValueError:
        return None


def write_drb3_panel(full_panel: Path, out_fa: Path) -> int:
    """Filter the full IPD panel down to BoLA-DRB3; return allele count."""
    n = 0
    keep = False
    with full_panel.open(encoding="utf-8") as src, out_fa.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                allele = parse_drb3_header(line)
                keep = allele is not None
                if keep:
                    n += 1
                    dst.write(line)
            elif keep:
                dst.write(line)
    if n == 0:
        raise Drb3WgsTypingError(f"no BoLA-DRB3 sequences in {full_panel}")
    return n


def bam_to_fasta(
    bam: Path,
    region: str,
    out_fa: Path,
    *,
    samtools: str,
) -> int:
    """Export primary DRB3-region reads as FASTA; return unique-read count."""
    view = subprocess.run(
        [samtools, "view", "-F", EXTRACT_FLAGS, str(bam), region],
        check=True,
        capture_output=True,
        text=True,
    )
    seen: set[str] = set()
    n = 0
    with out_fa.open("w", encoding="utf-8") as handle:
        for line in view.stdout.splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            qname = fields[0]
            seq = fields[9]
            if qname in seen:
                continue
            seen.add(qname)
            handle.write(f">{qname}\n{seq}\n")
            n += 1
    return n


def run_minimap2(
    panel_fa: Path,
    reads_fa: Path,
    out_paf: Path,
    *,
    minimap2: str,
    threads: int,
) -> str:
    with out_paf.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                minimap2,
                "-x",
                "map-ont",
                "-N",
                "50",
                "-p",
                "0.2",
                "--secondary=yes",
                "-t",
                str(threads),
                str(panel_fa),
                str(reads_fa),
            ],
            check=True,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    return out_paf.read_text(encoding="utf-8")


def parse_paf_best_per_read(
    paf_text: str,
    *,
    min_aln: int,
    min_identity: float,
) -> dict[str, AlleleHit]:
    """Return the best-identity DRB3 hit per read (ties broken by nmatch)."""
    best: dict[str, AlleleHit] = {}
    for line in paf_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        read = fields[0]
        allele_raw = fields[5].split()[0]
        allele = parse_drb3_header(allele_raw)
        if allele is None:
            continue
        nmatch, aln_len = int(fields[9]), int(fields[10])
        if aln_len < min_aln:
            continue
        divergence = None
        for tag in fields[12:]:
            if tag.startswith("dv:f:"):
                divergence = float(tag[5:])
                break
        identity = (1.0 - divergence) if divergence is not None else (nmatch / aln_len if aln_len else 0.0)
        if identity < min_identity:
            continue
        hit = AlleleHit(read=read, allele=allele, identity=identity, nmatch=nmatch, aln_len=aln_len)
        previous = best.get(read)
        if previous is None or (hit.identity, hit.nmatch) > (previous.identity, previous.nmatch):
            best[read] = hit
    return best


def format_top_alleles(counts: Counter[str], *, top_n: int = 3) -> str:
    if not counts:
        return "-"
    return ",".join(f"{allele}({n})" for allele, n in counts.most_common(top_n))


def known_presence_label(expected: set[str], observed: set[str], *, top_alleles: list[str]) -> str:
    """Describe whether PCR-SBT alleles appear among read-best hits (any rank)."""
    present = expected & observed
    if not expected:
        return "no_expected"
    if present == expected:
        if len(expected) == 1:
            return f"{next(iter(expected))} present"
        return "both"
    if present:
        allele = next(iter(present))
        in_top3 = allele in top_alleles[:3]
        suffix = "(top3)" if in_top3 else "(low)"
        return f"{allele} only{suffix}"
    return "none"


def summarise_read_typing(
    best_by_read: dict[str, AlleleHit],
    *,
    hap1: str,
    hap2: str,
    n_reads: int,
) -> dict[str, object]:
    expected = {normalize_drb3_allele(hap1), normalize_drb3_allele(hap2)}
    counts: Counter[str] = Counter(hit.allele for hit in best_by_read.values())
    top = [allele for allele, _n in counts.most_common()]
    top1 = top[0] if top else "-"
    observed = set(counts)
    presence = known_presence_label(expected, observed, top_alleles=top)
    top1_is_central_bias = top1 in CENTRAL_BIAS_ALLELES and top1 not in expected
    top1_matches_expected = top1 in expected
    exact_top_set = set(top[: len(expected)]) == expected and len(expected) <= len(top)
    return {
        "n_reads": n_reads,
        "n_reads_with_panel_hit": len(best_by_read),
        "top_alleles": format_top_alleles(counts),
        "top1_allele": top1,
        "top1_count": counts[top1] if top1 != "-" else 0,
        "known_alleles_present": presence,
        "top1_matches_expected": top1_matches_expected,
        "top1_is_central_bias_artifact": top1_is_central_bias,
        "exact_diplotype_from_top": exact_top_set,
        "allele_counts_json": json.dumps(dict(counts.most_common()), sort_keys=False),
    }


def load_sample_registry(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    analytical = [r for r in rows if r.get("analysis_status") == "analytical" and r.get("sample_id")]
    if not analytical:
        raise Drb3WgsTypingError(f"no analytical samples in {path}")
    return analytical


def load_assembly_concordance(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise Drb3WgsTypingError(f"missing assembly recovery TSV: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def render_report(
    read_rows: Sequence[dict[str, object]],
    assembly_rows: Sequence[dict[str, str]],
    *,
    region: str,
    min_identity: float,
    min_aln: int,
    n_panel_alleles: int,
) -> str:
    assembly_by_id = {r["sample_id"]: r for r in assembly_rows}
    n_bias = sum(1 for r in read_rows if r["top1_is_central_bias_artifact"])
    n_top1_ok = sum(1 for r in read_rows if r["top1_matches_expected"])
    n_exact_top = sum(1 for r in read_rows if r["exact_diplotype_from_top"])
    concordance_counts: Counter[str] = Counter(
        r.get("drb3_concordance", "missing") for r in assembly_rows
    )
    lines = [
        "# WGS-based DRB3 typing vs PCR-SBT (ARS-UCD2.0)",
        "",
        "Ports the read-based half of `docs/DRB3_WGS_CONCORDANCE.md` onto UCD2.0 "
        "target BAMs, and summarises assembly-based PCR-SBT concordance already "
        "computed by `dq_dr_recovery` / `ASSEMBLY_DQ_DR_REPORT.md`. The claim is "
        "self-limiting: this per-read scoring does not reproduce PCR-SBT DRB3 types. "
        "Superseded for interpretation by `03d_drb3_exon2_typing` (see the Conclusion).",
        "",
        "## Method (read-based)",
        "",
        f"- Region: `{region}` (primary alignments, `-F {EXTRACT_FLAGS}`)",
        f"- Panel: {n_panel_alleles} BoLA-DRB3 alleles from `ipd_mhc_bola_panel.fasta`",
        f"- Align: `minimap2 -x map-ont` (reads → panel), identity ≥ {min_identity:.2f}, "
        f"aln ≥ {min_aln} bp; per-read best allele by identity then nmatch",
        "",
        "## Read-based cohort summary",
        "",
        "| sample | PCR-SBT DRB3 | reads | panel hits | top alleles (count) | "
        "known present | top1 matches expected | central-bias artifact |",
        "| --- | --- | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in read_rows:
        lines.append(
            f"| {row['sample_id']} | {row['drb3']} | {row['n_reads']} | "
            f"{row['n_reads_with_panel_hit']} | {row['top_alleles']} | "
            f"{row['known_alleles_present']} | {row['top1_matches_expected']} | "
            f"{row['top1_is_central_bias_artifact']} |"
        )
    lines += [
        "",
        "## Read-based interpretation",
        "",
        f"- Animals whose top1 allele is a central-bias artifact (*016:01/*027:03 "
        f"not in PCR-SBT diplotype): **{n_bias}/{len(read_rows)}**",
        f"- Animals whose top1 allele is among PCR-SBT alleles: **{n_top1_ok}/{len(read_rows)}**",
        f"- Animals whose top-N allele set recovers the full PCR-SBT diplotype: "
        f"**{n_exact_top}/{len(read_rows)}**",
        "- True PCR-SBT alleles are often present at lower rank among reads. The "
        "*016:01 and *027:03 panel entries are full-length (801/780 bp) while most entries "
        "are exon 2 only (234-270 bp), so ranking by matching bases favours length; "
        "this is not sample mix-up or allele absence.",
        "",
        "## Assembly-based concordance (from DQ/DR local assembly)",
        "",
        "Source: `results/codex_reanalysis_v1/06_local_assembly/dq_dr/dq_dr_recovery.tsv` "
        "(IPD panel gate ≥95% identity, ≥80% coverage, MAPQ ≥20).",
        "",
        "| sample | PCR-SBT DRB3 | concordance | best DRB3 allele(s) |",
        "| --- | --- | --- | --- |",
    ]
    for row in read_rows:
        asm = assembly_by_id.get(str(row["sample_id"]), {})
        lines.append(
            f"| {row['sample_id']} | {row['drb3']} | "
            f"{asm.get('drb3_concordance', 'missing')} | "
            f"{asm.get('best_drb3_alleles') or '-'} |"
        )
    lines += [
        "",
        "### Assembly concordance counts",
        "",
    ]
    for label, n in sorted(concordance_counts.items()):
        lines.append(f"- `{label}`: {n}/{len(assembly_rows)}")
    homozygous = next((r for r in read_rows if r["drb3"] == "*010:01/*010:01"), None)
    lines += [
        "",
        "## Homozygous internal control (ID08 / ID08, *010:01/*010:01)",
        "",
    ]
    if homozygous is not None:
        asm = assembly_by_id.get(str(homozygous["sample_id"]), {})
        lines += [
            f"- Read-based top1: `{homozygous['top1_allele']}` "
            f"(n={homozygous['top1_count']}); "
            f"central-bias artifact={homozygous['top1_is_central_bias_artifact']}; "
            f"known presence=`{homozygous['known_alleles_present']}`",
            f"- DQ/DR assembly concordance: `{asm.get('drb3_concordance', 'missing')}` "
            f"(best `{asm.get('best_drb3_alleles') or '-'}`)",
            "- On ARS-UCD2.0, the homozygous animal does **not** recover *010:01 as the "
            "top call by either read-based mapping or DQ/DR assembly best-allele "
            "assignment — the central-allele bias applies even without heterozygous "
            "collapse. (The wide-interval assembly's `homozygous_control` label refers "
            "to single-probe placement mechanics, not IPD-panel allele identity.)",
            "",
        ]
    lines += [
        "## Conclusion",
        "",
        "Per-read best-allele scoring against the mixed-length IPD panel (and the DQ/DR "
        "assembly best-allele calls, which share the length preference) does not reproduce "
        "PCR-SBT DRB3 types. The failure is a scoring artifact: the exon-2 stage "
        "`03d_drb3_exon2_typing` scores every allele and read on one shared exon-2 window "
        "and recovers the PCR-SBT diplotype in most animals "
        "(`docs/codex_reanalysis/DRB3_EXON2_TYPING_REPORT.md`). The wide-interval probe "
        "offset (Table 7 / `ASSEMBLY_WIDE_REPORT.md`) concerns probe placement, not "
        "allele identity.",
        "",
    ]
    return "\n".join(lines)


def analyse_sample(
    sample_id: str,
    bam: Path,
    panel_fa: Path,
    work_dir: Path,
    *,
    samtools: str,
    minimap2: str,
    region: str,
    min_aln: int,
    min_identity: float,
    threads: int,
    hap1: str,
    hap2: str,
) -> dict[str, object]:
    work_dir.mkdir(parents=True, exist_ok=True)
    reads_fa = work_dir / f"{sample_id}.drb3.fasta"
    paf_path = work_dir / f"{sample_id}.reads_to_drb3.paf"
    n_reads = bam_to_fasta(bam, region, reads_fa, samtools=samtools)
    if n_reads == 0:
        return summarise_read_typing({}, hap1=hap1, hap2=hap2, n_reads=0)
    paf_text = run_minimap2(
        panel_fa, reads_fa, paf_path, minimap2=minimap2, threads=threads
    )
    best = parse_paf_best_per_read(paf_text, min_aln=min_aln, min_identity=min_identity)
    return summarise_read_typing(best, hap1=hap1, hap2=hap2, n_reads=n_reads)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument(
        "--assembly-recovery",
        type=Path,
        required=True,
        help="dq_dr_recovery.tsv from 06_local_assembly/dq_dr/",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--region", default=DRB3_REGION)
    parser.add_argument("--min-aln", type=int, default=DEFAULT_MIN_ALN)
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--samples", default="")
    args = parser.parse_args(list(argv) if argv is not None else None)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp"
    work.mkdir(exist_ok=True)

    panel_fa = work / "ipd_drb3_panel.fasta"
    n_panel = write_drb3_panel(args.panel, panel_fa)
    registry = load_sample_registry(args.registry)
    if args.samples.strip():
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        registry = [r for r in registry if r["sample_id"] in wanted]
        missing = wanted - {r["sample_id"] for r in registry}
        if missing:
            raise Drb3WgsTypingError(f"unknown sample_id(s): {sorted(missing)}")

    assembly_rows = load_assembly_concordance(args.assembly_recovery)
    rows: list[dict[str, object]] = []
    for entry in registry:
        sample_id = entry["sample_id"]
        bam = args.mapping_dir / f"{sample_id}.targets.bam"
        if not bam.exists():
            raise Drb3WgsTypingError(f"missing target BAM: {bam}")
        print(f"=== {sample_id} ===", flush=True)
        summary = analyse_sample(
            sample_id,
            bam,
            panel_fa,
            work / sample_id,
            samtools=args.samtools,
            minimap2=args.minimap2,
            region=args.region,
            min_aln=args.min_aln,
            min_identity=args.min_identity,
            threads=args.threads,
            hap1=entry["pcr_sbt_drb3_hap1"],
            hap2=entry["pcr_sbt_drb3_hap2"],
        )
        row = {
            "sample_id": sample_id,
            "drb3": f"{entry['pcr_sbt_drb3_hap1']}/{entry['pcr_sbt_drb3_hap2']}",
            **summary,
        }
        rows.append(row)
        print(
            f"  reads={summary['n_reads']} hits={summary['n_reads_with_panel_hit']} "
            f"top1={summary['top1_allele']}({summary['top1_count']}) "
            f"bias={summary['top1_is_central_bias_artifact']} "
            f"present={summary['known_alleles_present']}",
            flush=True,
        )

    tsv_path = out / "drb3_wgs_readbased_concordance.tsv"
    fieldnames = [
        "sample_id",
        "drb3",
        "n_reads",
        "n_reads_with_panel_hit",
        "top_alleles",
        "top1_allele",
        "top1_count",
        "known_alleles_present",
        "top1_matches_expected",
        "top1_is_central_bias_artifact",
        "exact_diplotype_from_top",
        "allele_counts_json",
    ]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    report = render_report(
        rows,
        assembly_rows,
        region=args.region,
        min_identity=args.min_identity,
        min_aln=args.min_aln,
        n_panel_alleles=n_panel,
    )
    report_path = Path("docs/codex_reanalysis/DRB3_WGS_TYPING_REPORT.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    (out / "DRB3_WGS_TYPING_REPORT.md").write_text(report, encoding="utf-8")

    status = {
        "status": "passed",
        "n_samples": len(rows),
        "region": args.region,
        "min_identity": args.min_identity,
        "min_aln": args.min_aln,
        "n_panel_alleles": n_panel,
        "tsv": str(tsv_path),
        "report": str(report_path),
        "n_central_bias_artifact_top1": sum(
            1 for r in rows if r["top1_is_central_bias_artifact"]
        ),
        "n_top1_matches_expected": sum(1 for r in rows if r["top1_matches_expected"]),
        "n_exact_diplotype_from_top": sum(1 for r in rows if r["exact_diplotype_from_top"]),
        "assembly_concordance_counts": dict(
            Counter(r.get("drb3_concordance", "missing") for r in assembly_rows)
        ),
    }
    (out / "drb3_wgs_typing_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    write_provenance(
        out,
        {
            "command": " ".join(sys.argv),
            "versions": {
                "python": sys.version.split()[0],
                "samtools": args.samtools,
                "minimap2": args.minimap2,
            },
            "input_hashes": {
                "registry": sha256_file(args.registry),
                "assembly_recovery": sha256_file(args.assembly_recovery),
            },
            "reference_hashes": {
                "ipd_panel": sha256_file(args.panel),
                "drb3_panel_subset": sha256_file(panel_fa),
            },
            "config_hash": "drb3_wgs_typing",
            "parameters": {
                "region": args.region,
                "min_identity": args.min_identity,
                "min_aln": args.min_aln,
                "bam_bytes": {
                    r["sample_id"]: (args.mapping_dir / f"{r['sample_id']}.targets.bam").stat().st_size
                    for r in rows
                },
            },
        },
    )
    print(f"Wrote {tsv_path}", flush=True)
    print(f"Wrote {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Drb3WgsTypingError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
