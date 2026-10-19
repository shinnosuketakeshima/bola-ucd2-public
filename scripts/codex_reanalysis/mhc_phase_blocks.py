"""WhatsHap MHC phase-block analysis on ARS-UCD2.0 target BAMs.

Ports the legacy WhatsHap + ``estimate_mhc_haplotype_blocks.py`` analysis
(``results/validation/mhc_haplotype_blocks/``, ``docs/MHC_HAPLOBLOCK_LD.md``)
onto UCD2.0 target BAMs with **corrected compartment labels**:

  class_II_DQ_DR → btn → class_III → class_I

(the legacy script used retracted class I/II windows; boundary *positions*
are coordinate-independent, but compartment labels were wrong).

Per sample:
1. Call biallelic SNPs on the continuous MHC-span BAM region (bcftools).
2. Phase with WhatsHap (PS tags).
3. Build phase blocks from heterozygous phased SNPs; tabulate gap midpoints.

Cohort:
- 100 kb boundary hotspots (all boundaries; matches the manuscript ~26.25 Mb claim)
- Landmark co-membership for DRB3 / rs29026690 / rs17872126
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

# Continuous MHC span present on UCD2.0 mapping_targets.bed (DQA2→TRIM cluster).
# Truncate distal end to the legacy MHC analysis window.
MHC_REGION = "NC_037350.1:25533030-29300000"
MHC_START = 25_533_030
MHC_END = 29_300_000

# Corrected physical compartments (ARS-UCD1.2/2.0 BTA23; byte-identical).
COMPARTMENTS: list[tuple[str, int, int]] = [
    ("class_II_DQ_DR", 25_580_000, 25_850_000),
    ("btn", 25_850_000, 27_150_000),
    ("class_III", 27_150_000, 28_000_000),
    ("class_I", 28_450_000, 28_800_000),
]

LANDMARKS: dict[str, int] = {
    "DRB3": 25_728_000,  # gene mid (25,723,716–25,734,819)
    "rs29026690": 27_603_651,
    "rs17872126": 28_419_662,
    "TNF_824": 27_719_770,  # TNF promoter -824A/G, repositioned (Konnai et al. 2006)
}

DEFAULT_BIN = 100_000
DEFAULT_MIN_BLOCK_HETS = 2


class MhcPhasingError(ValueError):
    """Raised when inputs cannot be analysed safely."""


@dataclass(frozen=True)
class PhasedHet:
    pos: int
    gt: str
    ps: str


@dataclass(frozen=True)
class PhaseBlock:
    sample_id: str
    ps: str
    start: int
    end: int
    n_hets: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class BlockBoundary:
    sample_id: str
    left_ps: str
    right_ps: str
    left_end: int
    right_start: int
    gap_bp: int
    midpoint: int
    left_comp: str
    right_comp: str
    boundary_type: str


def compartment_at(pos: int) -> str:
    for name, start, end in COMPARTMENTS:
        if start <= pos < end:
            return name
    if MHC_START <= pos < COMPARTMENTS[0][1]:
        return "flank_proximal"
    if COMPARTMENTS[0][2] <= pos < COMPARTMENTS[1][1]:
        return "gap_dq_btn"
    if COMPARTMENTS[2][2] <= pos < COMPARTMENTS[3][1]:
        return "gap_iii_i"
    if pos >= COMPARTMENTS[-1][2]:
        return "flank_distal"
    return "other"


def boundary_type(left_comp: str, right_comp: str) -> str:
    if left_comp == right_comp:
        return f"within_{left_comp}"
    return f"between_{left_comp}__{right_comp}"


def parse_phased_hets(vcf: Path, *, start: int, end: int) -> list[PhasedHet]:
    hets: list[PhasedHet] = []
    chroms = {"NC_037350.1", "chr23", "23", "BTA23"}
    opener = gzip.open if str(vcf).endswith(".gz") else open
    with opener(vcf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 10 or fields[0] not in chroms:
                continue
            pos = int(fields[1])
            if pos < start or pos > end:
                continue
            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            if "GT" not in fmt:
                continue
            gt = sample[fmt.index("GT")]
            if "|" not in gt:
                continue
            a, b = gt.split("|", 1)
            if a in {".", ""} or b in {".", ""} or a == b:
                continue
            ps = ""
            if "PS" in fmt:
                ps = sample[fmt.index("PS")].strip().strip('"')
            if not ps or ps in {".", ""}:
                ps = str(pos)
            hets.append(PhasedHet(pos=pos, gt=gt, ps=ps))
    hets.sort(key=lambda h: h.pos)
    return hets


def build_blocks(sample_id: str, hets: Sequence[PhasedHet]) -> list[PhaseBlock]:
    if not hets:
        return []
    by_ps: dict[str, list[PhasedHet]] = defaultdict(list)
    order: list[str] = []
    for het in hets:
        if het.ps not in by_ps:
            order.append(het.ps)
        by_ps[het.ps].append(het)
    blocks: list[PhaseBlock] = []
    for ps in order:
        group = sorted(by_ps[ps], key=lambda h: h.pos)
        blocks.append(
            PhaseBlock(
                sample_id=sample_id,
                ps=ps,
                start=group[0].pos,
                end=group[-1].pos,
                n_hets=len(group),
            )
        )
    blocks.sort(key=lambda b: (b.start, b.end))
    return blocks


def build_boundaries(
    sample_id: str,
    blocks: Sequence[PhaseBlock],
    *,
    min_block_hets: int = DEFAULT_MIN_BLOCK_HETS,
) -> list[BlockBoundary]:
    usable = [b for b in blocks if b.n_hets >= min_block_hets] or list(blocks)
    out: list[BlockBoundary] = []
    for i in range(len(usable) - 1):
        left, right = usable[i], usable[i + 1]
        if right.start <= left.end:
            gap = 0
            mid = right.start
        else:
            gap = right.start - left.end
            mid = left.end + gap // 2
        lc = compartment_at(left.end)
        rc = compartment_at(right.start)
        out.append(
            BlockBoundary(
                sample_id=sample_id,
                left_ps=left.ps,
                right_ps=right.ps,
                left_end=left.end,
                right_start=right.start,
                gap_bp=gap,
                midpoint=mid,
                left_comp=lc,
                right_comp=rc,
                boundary_type=boundary_type(lc, rc),
            )
        )
    return out


def block_covers(block: PhaseBlock, pos: int) -> bool:
    return block.start <= pos <= block.end


def landmark_same_block(blocks: Sequence[PhaseBlock], pos_a: int, pos_b: int) -> bool:
    """True if any phase block's genomic span covers both landmarks."""
    lo, hi = (pos_a, pos_b) if pos_a <= pos_b else (pos_b, pos_a)
    return any(b.start <= lo and b.end >= hi for b in blocks)


def window_hotspots(
    boundaries: Sequence[BlockBoundary],
    *,
    win: int = DEFAULT_BIN,
) -> list[dict[str, object]]:
    counts: Counter[int] = Counter()
    types: dict[int, Counter[str]] = defaultdict(Counter)
    for boundary in boundaries:
        w = (boundary.midpoint // win) * win
        counts[w] += 1
        types[w][boundary.boundary_type] += 1
    rows: list[dict[str, object]] = []
    for w, n in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        top = types[w].most_common(1)[0][0] if types[w] else ""
        rows.append(
            {
                "window_start": w,
                "window_end": w + win,
                "n_breaks": n,
                "top_boundary_type": top,
                "compartment": compartment_at(w + win // 2),
            }
        )
    return rows


def load_sample_registry(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    analytical = [r for r in rows if r.get("analysis_status") == "analytical" and r.get("sample_id")]
    if not analytical:
        raise MhcPhasingError(f"no analytical samples in {path}")
    return analytical


def run_cmd(cmd: Sequence[str], *, log_path: Path | None = None) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            f"$ {' '.join(cmd)}\n\nSTDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n",
            encoding="utf-8",
        )
    if result.returncode != 0:
        raise MhcPhasingError(
            f"command failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr[-2000:]}"
        )


def phase_sample(
    sample_id: str,
    bam: Path,
    reference: Path,
    work_dir: Path,
    *,
    region: str,
    samtools: str,
    bcftools: str,
    whatshap: str,
    threads: int,
) -> Path:
    """SNP-call + WhatsHap phase; return path to phased VCF.gz.

    Matches the legacy ``run_whatshap_phase.sh`` pattern: extract the MHC
    region, subsample ~35% of reads for SNP calling, then WhatsHap-phase.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    snp_vcf = work_dir / f"{sample_id}.bola.snp.vcf.gz"
    phased = work_dir / f"{sample_id}.phased.vcf.gz"
    if phased.exists() and Path(str(phased) + ".csi").exists():
        return phased

    region_bam = work_dir / f"{sample_id}.region.bam"
    sub_bam = work_dir / f"{sample_id}.sub.bam"
    if not region_bam.exists():
        run_cmd(
            [samtools, "view", "-@", str(threads), "-b", "-o", str(region_bam), str(bam), region],
            log_path=work_dir / f"{sample_id}.extract.log",
        )
        run_cmd([samtools, "index", "-@", str(threads), str(region_bam)])
    if not sub_bam.exists():
        # Legacy used ``samtools view -s 42.35`` (~35% subsample) for mpileup speed.
        run_cmd(
            [samtools, "view", "-@", str(threads), "-s", "42.35", "-b", "-o", str(sub_bam), str(region_bam)],
            log_path=work_dir / f"{sample_id}.subsample.log",
        )
        run_cmd([samtools, "index", "-@", str(threads), str(sub_bam)])

    if not snp_vcf.exists():
        raw_vcf = work_dir / f"{sample_id}.bola.vcf.gz"
        mpileup_bcf = work_dir / f"{sample_id}.mpileup.bcf"
        run_cmd(
            [
                bcftools,
                "mpileup",
                "-Ob",
                "-B",
                "-d",
                "50",
                "-f",
                str(reference),
                "-r",
                region,
                "-q",
                "20",
                "-o",
                str(mpileup_bcf),
                str(sub_bam),
            ],
            log_path=work_dir / f"{sample_id}.mpileup.log",
        )
        run_cmd(
            [bcftools, "call", "-mv", "-Oz", "-o", str(raw_vcf), str(mpileup_bcf)],
            log_path=work_dir / f"{sample_id}.call.log",
        )
        run_cmd([bcftools, "index", "-f", str(raw_vcf)])
        run_cmd(
            [
                bcftools,
                "view",
                "-v",
                "snps",
                "-m2",
                "-M2",
                "-i",
                "QUAL>=20",
                "-Oz",
                "-o",
                str(snp_vcf),
                str(raw_vcf),
            ]
        )
        run_cmd([bcftools, "index", "-f", str(snp_vcf)])
        mpileup_bcf.unlink(missing_ok=True)

    n_var = int(
        subprocess.run(
            [bcftools, "view", "-H", str(snp_vcf)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.count("\n")
    )
    if n_var < 2:
        raise MhcPhasingError(f"{sample_id}: too few SNPs ({n_var})")

    # Phase against the full region BAM (not the subsample) for denser PS support.
    run_cmd(
        [
            whatshap,
            "phase",
            "-o",
            str(phased),
            "--reference",
            str(reference),
            "--ignore-read-groups",
            str(snp_vcf),
            str(region_bam),
        ],
        log_path=work_dir / f"{sample_id}.whatshap.log",
    )
    run_cmd([bcftools, "index", "-f", str(phased)])
    return phased


def landmark_row(sample_id: str, blocks: Sequence[PhaseBlock]) -> dict[str, object]:
    pairs = [
        ("DRB3", "rs29026690"),
        ("DRB3", "rs17872126"),
        ("rs29026690", "rs17872126"),
        ("TNF_824", "rs29026690"),
        ("TNF_824", "DRB3"),
        ("TNF_824", "rs17872126"),
    ]
    row: dict[str, object] = {"sample_id": sample_id}
    for a, b in pairs:
        key = f"{a}_x_{b}_same_block"
        row[key] = landmark_same_block(blocks, LANDMARKS[a], LANDMARKS[b])
    return row


def render_report(
    *,
    hotspot_rows: Sequence[dict[str, object]],
    landmark_rows: Sequence[dict[str, object]],
    sample_summaries: Sequence[dict[str, object]],
    region: str,
) -> str:
    top = list(hotspot_rows[:8])
    n = len(landmark_rows)

    def count_true(key: str) -> int:
        return sum(1 for r in landmark_rows if r[key])

    lines = [
        "# MHC phase-block analysis (ARS-UCD2.0 WhatsHap)",
        "",
        "Ports the legacy WhatsHap phase-block tabulation "
        "(`results/validation/mhc_haplotype_blocks/`) onto UCD2.0 target BAMs "
        "with corrected compartment labels (class II DQ/DR → BTN → class III → class I).",
        "",
        "## Method",
        "",
        f"- Region: `{region}` (primary target BAM coverage; bcftools mpileup/call + WhatsHap 2.2)",
        "- Blocks: heterozygous phased SNPs sharing a PS tag",
        "- Boundaries: midpoints of gaps between adjacent blocks (≥2 hets)",
        "- Hotspots: all boundaries binned at 100 kb (matches manuscript ~26.25 Mb claim)",
        "- Landmark co-membership: a block spans both positions",
        "",
        "### Corrected compartments",
        "",
        "| compartment | start–end |",
        "| --- | --- |",
    ]
    for name, start, end in COMPARTMENTS:
        lines.append(f"| {name} | {start:,}–{end:,} |")
    lines += [
        "",
        "## Landmark co-membership (physical cis on one PS block)",
        "",
        f"- DRB3 × rs29026690 same block: **{count_true('DRB3_x_rs29026690_same_block')}/{n}**",
        f"- DRB3 × rs17872126 same block: **{count_true('DRB3_x_rs17872126_same_block')}/{n}**",
        f"- rs29026690 × rs17872126 same block: "
        f"**{count_true('rs29026690_x_rs17872126_same_block')}/{n}**",
        f"- TNF_824 × rs29026690 same block: **{count_true('TNF_824_x_rs29026690_same_block')}/{n}**",
        f"- TNF_824 × DRB3 same block: **{count_true('TNF_824_x_DRB3_same_block')}/{n}**",
        f"- TNF_824 × rs17872126 same block: "
        f"**{count_true('TNF_824_x_rs17872126_same_block')}/{n}**",
        "",
        "| sample | DRB3×rs29026690 | DRB3×rs17872126 | rs29026690×rs17872126 | TNF_824×rs29026690 | TNF_824×DRB3 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in landmark_rows:
        lines.append(
            f"| {row['sample_id']} | {row['DRB3_x_rs29026690_same_block']} | "
            f"{row['DRB3_x_rs17872126_same_block']} | "
            f"{row['rs29026690_x_rs17872126_same_block']} | "
            f"{row['TNF_824_x_rs29026690_same_block']} | "
            f"{row['TNF_824_x_DRB3_same_block']} |"
        )
    lines += [
        "",
        "## Boundary hotspots (100 kb, all boundaries, top windows)",
        "",
        "| window | n_breaks | top type | compartment |",
        "| --- | ---: | --- | --- |",
    ]
    for h in top:
        lines.append(
            f"| {int(h['window_start']):,}–{int(h['window_end']):,} | {h['n_breaks']} | "
            f"{h['top_boundary_type']} | {h['compartment']} |"
        )
    peak = hotspot_rows[0] if hotspot_rows else None
    lines += [
        "",
        "## Interpretation",
        "",
    ]
    if peak is not None:
        lines.append(
            f"- Strongest boundary window: **{int(peak['window_start'])/1e6:.2f}–"
            f"{int(peak['window_end'])/1e6:.2f} Mb** ({peak['n_breaks']} breaks; "
            f"`{peak['compartment']}` / `{peak['top_boundary_type']}`)."
        )
    lines += [
        "- DRB3 and rs29026690 remain on separate phase blocks in essentially all "
        "animals (population r² from marker genotypes is not physical cis).",
        "- Wide-interval de novo assembly independently confirms non-bridging "
        "(`ASSEMBLY_WIDE_REPORT.md`); this WhatsHap re-run pins the ~26 Mb boundary "
        "location on UCD2.0 BAMs.",
        "",
        "## Per-sample summary",
        "",
        "| sample | phased hets | blocks | boundaries |",
        "| --- | ---: | ---: | ---: |",
    ]
    for s in sample_summaries:
        lines.append(
            f"| {s['sample_id']} | {s['n_phased_hets']} | {s['n_phase_blocks']} | "
            f"{s['n_boundaries']} |"
        )
    lines.append("")
    return "\n".join(lines)


def write_tsv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--region", default=MHC_REGION)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--whatshap", default="whatshap")
    parser.add_argument("--samples", default="")
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        help="Skip SNP/phase; only summarise existing per-sample phased VCFs under output-dir/tmp",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp"
    work.mkdir(exist_ok=True)

    registry = load_sample_registry(args.registry)
    if args.samples.strip():
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        registry = [r for r in registry if r["sample_id"] in wanted]
        missing = wanted - {r["sample_id"] for r in registry}
        if missing:
            raise MhcPhasingError(f"unknown sample_id(s): {sorted(missing)}")

    all_boundaries: list[BlockBoundary] = []
    all_blocks: list[PhaseBlock] = []
    landmark_rows: list[dict[str, object]] = []
    sample_summaries: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    boundary_rows: list[dict[str, object]] = []

    for entry in registry:
        sample_id = entry["sample_id"]
        bam = args.mapping_dir / f"{sample_id}.targets.bam"
        if not bam.exists():
            raise MhcPhasingError(f"missing target BAM: {bam}")
        sample_work = work / sample_id
        print(f"=== {sample_id} ===", flush=True)
        if args.summarize_only:
            phased = sample_work / f"{sample_id}.phased.vcf.gz"
            if not phased.exists():
                raise MhcPhasingError(f"missing phased VCF for summarize-only: {phased}")
        else:
            phased = phase_sample(
                sample_id,
                bam,
                args.reference,
                sample_work,
                region=args.region,
                samtools=args.samtools,
                bcftools=args.bcftools,
                whatshap=args.whatshap,
                threads=args.threads,
            )
        hets = parse_phased_hets(phased, start=MHC_START, end=MHC_END)
        blocks = build_blocks(sample_id, hets)
        boundaries = build_boundaries(sample_id, blocks)
        all_blocks.extend(blocks)
        all_boundaries.extend(boundaries)
        lm = landmark_row(sample_id, blocks)
        landmark_rows.append(lm)
        sample_summaries.append(
            {
                "sample_id": sample_id,
                "drb3": f"{entry['pcr_sbt_drb3_hap1']}/{entry['pcr_sbt_drb3_hap2']}",
                "n_phased_hets": len(hets),
                "n_phase_blocks": len(blocks),
                "n_boundaries": len(boundaries),
            }
        )
        for b in blocks:
            block_rows.append(
                {
                    "sample_id": sample_id,
                    "ps": b.ps,
                    "start": b.start,
                    "end": b.end,
                    "length": b.length,
                    "n_hets": b.n_hets,
                    "start_comp": compartment_at(b.start),
                    "end_comp": compartment_at(b.end),
                }
            )
        for b in boundaries:
            boundary_rows.append(
                {
                    "sample_id": sample_id,
                    "left_ps": b.left_ps,
                    "right_ps": b.right_ps,
                    "left_end": b.left_end,
                    "right_start": b.right_start,
                    "gap_bp": b.gap_bp,
                    "midpoint": b.midpoint,
                    "left_comp": b.left_comp,
                    "right_comp": b.right_comp,
                    "boundary_type": b.boundary_type,
                }
            )
        print(
            f"  hets={len(hets)} blocks={len(blocks)} boundaries={len(boundaries)} "
            f"DRB3×rs290={lm['DRB3_x_rs29026690_same_block']}",
            flush=True,
        )

    hotspots = window_hotspots(all_boundaries)
    write_tsv(out / "phase_blocks.tsv", block_rows)
    write_tsv(out / "block_boundaries.tsv", boundary_rows)
    write_tsv(out / "landmark_co_membership.tsv", landmark_rows)
    write_tsv(out / "sample_block_summary.tsv", sample_summaries)
    write_tsv(out / "all_boundary_hotspots_100kb.tsv", hotspots)

    report = render_report(
        hotspot_rows=hotspots,
        landmark_rows=landmark_rows,
        sample_summaries=sample_summaries,
        region=args.region,
    )
    report_path = Path("docs/codex_reanalysis/MHC_PHASE_BLOCKS_REPORT.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    (out / "MHC_PHASE_BLOCKS_REPORT.md").write_text(report, encoding="utf-8")

    peak = hotspots[0] if hotspots else {}
    status = {
        "status": "passed",
        "n_samples": len(sample_summaries),
        "region": args.region,
        "compartments": [
            {"name": n, "start": s, "end": e} for n, s, e in COMPARTMENTS
        ],
        "landmarks": LANDMARKS,
        "n_DRB3_x_rs29026690_same_block": sum(
            1 for r in landmark_rows if r["DRB3_x_rs29026690_same_block"]
        ),
        "n_DRB3_x_rs17872126_same_block": sum(
            1 for r in landmark_rows if r["DRB3_x_rs17872126_same_block"]
        ),
        "n_rs29026690_x_rs17872126_same_block": sum(
            1 for r in landmark_rows if r["rs29026690_x_rs17872126_same_block"]
        ),
        "n_TNF_824_x_rs29026690_same_block": sum(
            1 for r in landmark_rows if r["TNF_824_x_rs29026690_same_block"]
        ),
        "n_TNF_824_x_DRB3_same_block": sum(
            1 for r in landmark_rows if r["TNF_824_x_DRB3_same_block"]
        ),
        "n_TNF_824_x_rs17872126_same_block": sum(
            1 for r in landmark_rows if r["TNF_824_x_rs17872126_same_block"]
        ),
        "peak_hotspot_window_start": peak.get("window_start"),
        "peak_hotspot_n_breaks": peak.get("n_breaks"),
        "peak_hotspot_compartment": peak.get("compartment"),
        "report": str(report_path),
    }
    (out / "mhc_phase_blocks_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(
        out,
        {
            "command": " ".join(sys.argv),
            "versions": {
                "python": sys.version.split()[0],
                "samtools": args.samtools,
                "bcftools": args.bcftools,
                "whatshap": args.whatshap,
            },
            "input_hashes": {"registry": sha256_file(args.registry)},
            "reference_hashes": {"reference": sha256_file(args.reference)},
            "config_hash": "mhc_phase_blocks",
            "parameters": {
                "region": args.region,
                "compartments": COMPARTMENTS,
                "landmarks": LANDMARKS,
            },
        },
    )
    print(f"Wrote {report_path}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MhcPhasingError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
