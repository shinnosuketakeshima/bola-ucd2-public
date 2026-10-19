"""Read-based in-cis DQ/DQB copy counting on ARS-UCD2.0 target BAMs.

Ports the legacy scratch analysis (`dq_copycount.sh` / `dq_cohort.py` →
`results/cnv_corrected/dq_readbased_copynumber.tsv`) described in
`docs/HANDOFF_PLAN_PROGRESS.md` §2e and `METHODS_EN.md`:

1. Extract primary alignments from the real DQ window on the UCD2.0 target BAM.
2. Keep long reads (≥ ``--min-read-len``).
3. Map those reads to the IPD-MHC BoLA DQA/DQB/DRB3 panel allele sequences.
4. On each read, cluster same-family hits whose query intervals fall within a
   gene-span merge window (default 15 kb). Exons of one gene merge; tandem
   paralogs farther apart stay separate copies.
5. Summarise the fraction of long reads with ≥2 in-cis DQA or DQB copies.

The manuscript claim is deliberately self-limiting: nearly all reads carry
DQA=1 / DQB=1 because DQ paralogs sit 49–141 kb apart while median DQ-region
read length is 22–30 kb. This script reproduces that negative result on the
UCD2.0 BAMs rather than inventing a positive copy-number call.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

# HANDOFF §2e real DQ window (DQA2–DRB3), 0-based half-open → samtools 1-based closed.
DQ_REGION = "NC_037350.1:25550000-25720000"
DEFAULT_GENE_SPAN = 15_000
DEFAULT_MIN_READ_LEN = 15_000
DEFAULT_MIN_ALN = 200
DEFAULT_MIN_IDENTITY = 0.70
FAMILIES = ("DQA", "DQB")
PANEL_FAMILIES = frozenset({"DQA", "DQB", "DRB3"})
EXTRACT_FLAGS = "0x900"


class DqReadCopycountError(ValueError):
    """Raised when inputs cannot be analysed safely."""


@dataclass(frozen=True)
class Hit:
    read: str
    family: str
    allele: str
    qstart: int
    qend: int
    identity: float
    aln_len: int


def parse_panel_family(header: str) -> str | None:
    """Return DQA/DQB/DRB3 from an IPD FASTA header, else None.

    Accepts ``>BoLA-DQA*001:01 ...`` and ignores non-BoLA / class-I loci.
    """
    token = header[1:].split()[0] if header.startswith(">") else header.split()[0]
    if not token.startswith("BoLA-"):
        return None
    body = token[5:]
    if "*" not in body:
        return None
    family = body.split("*", 1)[0]
    if family in PANEL_FAMILIES:
        return family
    return None


def write_classii_panel(full_panel: Path, out_fa: Path) -> dict[str, int]:
    """Filter the full IPD panel down to BoLA-DQA / BoLA-DQB / BoLA-DRB3."""
    counts = {family: 0 for family in sorted(PANEL_FAMILIES)}
    keep = False
    family: str | None = None
    with full_panel.open(encoding="utf-8") as src, out_fa.open("w", encoding="utf-8") as dst:
        for line in src:
            if line.startswith(">"):
                family = parse_panel_family(line)
                keep = family is not None
                if keep:
                    counts[family] += 1
                    dst.write(line)
            elif keep:
                dst.write(line)
    if sum(counts.values()) == 0:
        raise DqReadCopycountError(f"no BoLA-DQA/DQB/DRB3 sequences in {full_panel}")
    return counts


def merge_intervals(intervals: Sequence[tuple[int, int]], gap: int) -> list[tuple[int, int]]:
    """Merge intervals that overlap or sit within ``gap`` bp of each other."""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def count_family_copies(
    hits: Sequence[Hit],
    *,
    gene_span: int,
) -> dict[str, int]:
    """Return {family: n_copies} after gene-span clustering on one read."""
    by_family: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for hit in hits:
        if hit.family in FAMILIES:
            by_family[hit.family].append((hit.qstart, hit.qend))
    return {family: len(merge_intervals(by_family.get(family, []), gene_span)) for family in FAMILIES}


def parse_paf_hits(
    paf_text: str,
    *,
    min_aln: int,
    min_identity: float,
) -> dict[str, list[Hit]]:
    """Parse ``minimap2 panel.fa reads.fa`` PAF into per-read Hit lists.

    Query = read, target = panel allele. Identity prefers ``dv:f:`` when present.
    """
    by_read: dict[str, list[Hit]] = defaultdict(list)
    for line in paf_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        read = fields[0]
        allele = fields[5].split()[0]
        family = parse_panel_family(allele)
        if family is None:
            continue
        qstart, qend = int(fields[2]), int(fields[3])
        if qend < qstart:
            qstart, qend = qend, qstart
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
        by_read[read].append(
            Hit(
                read=read,
                family=family,
                allele=allele,
                qstart=qstart,
                qend=qend,
                identity=identity,
                aln_len=aln_len,
            )
        )
    return by_read


def summarise_sample(
    by_read: dict[str, list[Hit]],
    *,
    n_long_reads: int,
    gene_span: int,
) -> dict[str, float | int]:
    """Build the legacy-compatible summary row fields for one sample.

    ``pct_*_ge2`` is percent of *long reads* (legacy denominator), not of
    reads that hit the panel — matching ``dq_readbased_copynumber.tsv``.
    """
    max_dqa = 0
    max_dqb = 0
    n_dqa_ge2 = 0
    n_dqb_ge2 = 0
    n_with_family_hit = 0
    for hits in by_read.values():
        copies = count_family_copies(hits, gene_span=gene_span)
        if copies["DQA"] == 0 and copies["DQB"] == 0:
            continue
        n_with_family_hit += 1
        max_dqa = max(max_dqa, copies["DQA"])
        max_dqb = max(max_dqb, copies["DQB"])
        if copies["DQA"] >= 2:
            n_dqa_ge2 += 1
        if copies["DQB"] >= 2:
            n_dqb_ge2 += 1
    denom = n_long_reads
    return {
        "long_reads": n_long_reads,
        "reads_with_dqa_or_dqb_hit": n_with_family_hit,
        "pct_DQA_ge2": round(100.0 * n_dqa_ge2 / denom, 1) if denom else 0.0,
        "max_DQA_cis": max_dqa,
        "pct_DQB_ge2": round(100.0 * n_dqb_ge2 / denom, 1) if denom else 0.0,
        "max_DQB_cis": max_dqb,
        "n_DQA_ge2": n_dqa_ge2,
        "n_DQB_ge2": n_dqb_ge2,
    }


def load_sample_registry(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    analytical = [r for r in rows if r.get("analysis_status") == "analytical" and r.get("sample_id")]
    if not analytical:
        raise DqReadCopycountError(f"no analytical samples in {path}")
    return analytical


def format_drb3(hap1: str, hap2: str) -> str:
    return f"{hap1}/{hap2}"


def bam_to_long_fasta(
    bam: Path,
    region: str,
    out_fa: Path,
    *,
    samtools: str,
    min_read_len: int,
) -> int:
    """Export primary DQ-region reads ≥ min_read_len as FASTA; return count."""
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
            if qname in seen or len(seq) < min_read_len:
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
    """Map reads (query) to the class-II panel (target); return PAF text."""
    with out_paf.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
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
    _ = result  # stderr kept for debugging if we re-raise
    return out_paf.read_text(encoding="utf-8")


def render_report(
    rows: Sequence[dict[str, object]],
    *,
    region: str,
    gene_span: int,
    min_read_len: int,
    min_identity: float,
    min_aln: int,
) -> str:
    lines = [
        "# Read-based in-cis DQ/DQB copy counting (ARS-UCD2.0)",
        "",
        "Ports the legacy analysis behind `results/cnv_corrected/dq_readbased_copynumber.tsv` "
        "(`HANDOFF_PLAN_PROGRESS.md` §2e) onto UCD2.0 target BAMs. The claim is "
        "self-limiting: this method cannot resolve DQ copy number because paralogs "
        "sit farther apart than typical read lengths.",
        "",
        "## Method",
        "",
        f"- Region: `{region}` (primary alignments, `-F {EXTRACT_FLAGS}`)",
        f"- Reads: length ≥ {min_read_len:,} bp",
        f"- Panel: IPD-MHC BoLA-DQA / BoLA-DQB / BoLA-DRB3 from `ipd_mhc_bola_panel.fasta`",
        f"- Align: `minimap2 -x map-ont` (reads → panel), identity ≥ {min_identity:.2f}, "
        f"aln ≥ {min_aln} bp",
        f"- Copy call: merge same-family query intervals within {gene_span:,} bp "
        "(gene-span correction; naive ~200 bp clustering double-counts exons)",
        "",
        "## Cohort summary",
        "",
        "| sample | PCR-SBT DRB3 | long reads | reads w/ DQA|DQB hit | "
        "% DQA≥2 | max DQA cis | % DQB≥2 | max DQB cis |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample_id']} | {row['drb3']} | {row['long_reads']} | "
            f"{row['reads_with_dqa_or_dqb_hit']} | {row['pct_DQA_ge2']} | "
            f"{row['max_DQA_cis']} | {row['pct_DQB_ge2']} | {row['max_DQB_cis']} |"
        )
    n_any_dqa = sum(1 for r in rows if int(r["max_DQA_cis"]) >= 2)
    n_any_dqb = sum(1 for r in rows if int(r["max_DQB_cis"]) >= 2)
    lines += [
        "",
        "## Interpretation",
        "",
        f"- Animals with any in-cis DQA≥2 read: **{n_any_dqa}/{len(rows)}**",
        f"- Animals with any in-cis DQB≥2 read: **{n_any_dqb}/{len(rows)}**",
        "- Percentages are near zero for nearly every animal; rare multicopy reads "
        "do not resolve array copy number. Haplotype-keeping local assembly remains "
        "the method that recovers DQ structure (`ASSEMBLY_DQ_DR_REPORT.md`).",
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
    min_read_len: int,
    gene_span: int,
    min_aln: int,
    min_identity: float,
    threads: int,
) -> dict[str, float | int]:
    work_dir.mkdir(parents=True, exist_ok=True)
    reads_fa = work_dir / f"{sample_id}.dq_long.fasta"
    paf_path = work_dir / f"{sample_id}.reads_to_panel.paf"
    n_long = bam_to_long_fasta(
        bam, region, reads_fa, samtools=samtools, min_read_len=min_read_len
    )
    if n_long == 0:
        return {
            "long_reads": 0,
            "reads_with_dqa_or_dqb_hit": 0,
            "pct_DQA_ge2": 0.0,
            "max_DQA_cis": 0,
            "pct_DQB_ge2": 0.0,
            "max_DQB_cis": 0,
            "n_DQA_ge2": 0,
            "n_DQB_ge2": 0,
        }
    paf_text = run_minimap2(
        panel_fa, reads_fa, paf_path, minimap2=minimap2, threads=threads
    )
    by_read = parse_paf_hits(paf_text, min_aln=min_aln, min_identity=min_identity)
    return summarise_sample(by_read, n_long_reads=n_long, gene_span=gene_span)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, required=True, help="02_mapping directory")
    parser.add_argument("--panel", type=Path, required=True, help="ipd_mhc_bola_panel.fasta")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--region", default=DQ_REGION)
    parser.add_argument("--min-read-len", type=int, default=DEFAULT_MIN_READ_LEN)
    parser.add_argument("--gene-span", type=int, default=DEFAULT_GENE_SPAN)
    parser.add_argument("--min-aln", type=int, default=DEFAULT_MIN_ALN)
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_IDENTITY)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument(
        "--samples",
        default="",
        help="Optional comma-separated sample_id subset (default: all analytical)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp"
    work.mkdir(exist_ok=True)

    panel_fa = work / "ipd_classii_dqa_dqb_drb3.fasta"
    panel_counts = write_classii_panel(args.panel, panel_fa)

    registry = load_sample_registry(args.registry)
    if args.samples.strip():
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        registry = [r for r in registry if r["sample_id"] in wanted]
        missing = wanted - {r["sample_id"] for r in registry}
        if missing:
            raise DqReadCopycountError(f"unknown sample_id(s): {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for entry in registry:
        sample_id = entry["sample_id"]
        bam = args.mapping_dir / f"{sample_id}.targets.bam"
        if not bam.exists():
            raise DqReadCopycountError(f"missing target BAM: {bam}")
        print(f"=== {sample_id} ===", flush=True)
        summary = analyse_sample(
            sample_id,
            bam,
            panel_fa,
            work / sample_id,
            samtools=args.samtools,
            minimap2=args.minimap2,
            region=args.region,
            min_read_len=args.min_read_len,
            gene_span=args.gene_span,
            min_aln=args.min_aln,
            min_identity=args.min_identity,
            threads=args.threads,
        )
        row = {
            "sample_id": sample_id,
            "drb3": format_drb3(entry["pcr_sbt_drb3_hap1"], entry["pcr_sbt_drb3_hap2"]),
            **summary,
        }
        rows.append(row)
        print(
            f"  long={summary['long_reads']} hit={summary['reads_with_dqa_or_dqb_hit']} "
            f"DQA%≥2={summary['pct_DQA_ge2']} maxDQA={summary['max_DQA_cis']} "
            f"DQB%≥2={summary['pct_DQB_ge2']} maxDQB={summary['max_DQB_cis']}",
            flush=True,
        )

    tsv_path = out / "dq_readbased_copynumber.tsv"
    fieldnames = [
        "sample_id",
        "drb3",
        "long_reads",
        "reads_with_dqa_or_dqb_hit",
        "pct_DQA_ge2",
        "max_DQA_cis",
        "pct_DQB_ge2",
        "max_DQB_cis",
        "n_DQA_ge2",
        "n_DQB_ge2",
    ]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    report = render_report(
        rows,
        region=args.region,
        gene_span=args.gene_span,
        min_read_len=args.min_read_len,
        min_identity=args.min_identity,
        min_aln=args.min_aln,
    )
    report_path = Path("docs/codex_reanalysis/DQ_READ_COPYCOUNT_REPORT.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    (out / "DQ_READ_COPYCOUNT_REPORT.md").write_text(report, encoding="utf-8")

    status = {
        "status": "passed",
        "n_samples": len(rows),
        "region": args.region,
        "gene_span": args.gene_span,
        "min_read_len": args.min_read_len,
        "min_identity": args.min_identity,
        "min_aln": args.min_aln,
        "panel_allele_counts": panel_counts,
        "tsv": str(tsv_path),
        "report": str(report_path),
        "animals_with_any_dqa_ge2": sum(1 for r in rows if int(r["max_DQA_cis"]) >= 2),
        "animals_with_any_dqb_ge2": sum(1 for r in rows if int(r["max_DQB_cis"]) >= 2),
    }
    (out / "dq_read_copycount_status.json").write_text(
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
            },
            "reference_hashes": {
                "ipd_panel": sha256_file(args.panel),
                "classii_panel_subset": sha256_file(panel_fa),
            },
            "config_hash": "dq_read_copycount",
            "parameters": {
                "region": args.region,
                "gene_span": args.gene_span,
                "min_read_len": args.min_read_len,
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
    except DqReadCopycountError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
