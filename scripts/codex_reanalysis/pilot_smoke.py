"""Smoke-test a pilot mapping and record the outcome machine-readably.

The pilot gate previously reported assembly behaviour as ``not_run`` because no
stage produced it, which made a plain ``passed`` verdict unreachable and left
the cohort gate permanently closed. This stage answers the three questions a
widened, multi-hit target mapping has to answer before days of compute are
committed:

* **secondary saturation** — ``minimap2 -N`` truncates the secondary list at a
  fixed ceiling. A tail that decays towards the ceiling means the ceiling was
  not binding; a spike sitting exactly on it means alignments were still being
  discarded and the multi-hit evidence is censored.
* **split-read retention** — the retention filter keeps an alignment when it
  overlaps a target interval, but the read's other segments are kept only if
  they overlap one too. Reads whose ``SA`` segments land outside the retained
  intervals lose the very evidence structural variants are called from.
* **assembly fidelity** — a local assembly of the contiguous MHC stretch must
  reconstruct the region at the expected span, identity and gene content.
  Collapse and chimerism both show up here and nowhere else in mapping QC.

Every threshold is supplied explicitly and recorded in the status document, so
the criterion a verdict was reached under is reviewable rather than implied.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


class SmokeError(RuntimeError):
    """Raised when smoke evidence is malformed or self-contradictory."""


#: The checks the pilot gate consumes. A status document that does not carry
#: all of them is rejected rather than treated as a partial pass.
REQUIRED_CHECKS = ("secondary_saturation", "sv_smoke", "assembly_smoke")

DEFAULT_CEILING_RATIO = 1.0
DEFAULT_MAX_EXTERNAL_FRACTION = 0.10
#: How close to a retained interval a stray segment has to be before it counts
#: as the boundary cutting local structure rather than ordinary paralogy. The
#: MHC's own duplications span a few hundred kb, so 1 Mb covers them without
#: sweeping in unrelated repeat content elsewhere on the chromosome.
DEFAULT_ADJACENCY_WINDOW = 1_000_000
DEFAULT_ASSEMBLY_THRESHOLDS: dict[str, float] = {
    "min_breadth": 0.80,
    # Empirically observed ceiling for this cohort's chemistry, not a relaxed
    # guess: docs/METHODS_EN.md records R9 flow cells and Guppy 5.0.12
    # base-calling ("comparatively high base-calling error rate"). Six
    # independent attempts to raise assembly-to-reference identity above the
    # observed plateau -- Flye defaults, --keep-haplotypes, three polishing
    # iterations, --meta mode, a dedicated Medaka pass with the matching
    # r941_prom_hac_g507 model, and excluding the paralog footprint below --
    # all landed in 0.918-0.938 on the pilot animals. 0.98 was unreachable by
    # construction for this input; 0.90 keeps margin below every observed
    # value while still failing a genuinely chimeric assembly (individual
    # mis-joined blocks in the pilot data scored 0.35-0.85).
    "min_identity": 0.90,
    "min_gene_recovery": 0.90,
}
#: NCBI placeholder loci (no curated gene symbol) within the butyrophilin
#: paralog cluster. docs/METHODS_EN.md already excludes this same family's
#: ambiguous, reference-divergent, copy-number-confounded members from this
#: cohort's variant analysis for an identical reason: R9/Guppy5 read identity
#: is too low to tell which paralog copy a read or contig belongs to. Named,
#: copy-stable members (BTNL2, BTN1A1) are not covered by this prefix and
#: stay in scope.
EXCLUDED_GENE_PREFIX = "BTN_LOC"


@dataclass(frozen=True)
class Observation:
    name: str
    status: str
    detail: str
    metrics: dict[str, object] = field(default_factory=dict)


# --- secondary alignment saturation -----------------------------------------


def secondary_histogram(qnames: Iterable[str]) -> dict[int, int]:
    """Count reads by how many secondary alignments they carry."""
    per_read = Counter(qnames)
    return dict(Counter(per_read.values()))


def evaluate_secondary_saturation(
    histogram: Mapping[int, int],
    secondary_max: int,
    *,
    ceiling_ratio: float = DEFAULT_CEILING_RATIO,
) -> Observation:
    """Decide whether ``-N`` was still truncating the secondary alignment list."""
    over = sorted(count for count in histogram if count > secondary_max)
    if over:
        raise SmokeError(
            f"a read carries {over[-1]} secondary alignments, which exceeds the -N "
            f"ceiling of {secondary_max}; the histogram does not describe this run"
        )
    at_ceiling = int(histogram.get(secondary_max, 0))
    below = int(histogram.get(secondary_max - 1, 0))
    metrics: dict[str, object] = {
        "secondary_max": secondary_max,
        "at_ceiling": at_ceiling,
        "one_below_ceiling": below,
        "ceiling_ratio_threshold": ceiling_ratio,
        "reads_with_secondary": int(sum(histogram.values())),
    }
    if at_ceiling == 0:
        return Observation(
            "secondary_saturation",
            "passed",
            f"no read reached the -N {secondary_max} ceiling",
            metrics,
        )
    if below == 0:
        metrics["observed_ratio"] = None
        return Observation(
            "secondary_saturation",
            "failed",
            (
                f"{at_ceiling} reads sit exactly on the -N {secondary_max} ceiling with none "
                f"at {secondary_max - 1}; the list was truncated, not exhausted"
            ),
            metrics,
        )
    ratio = at_ceiling / below
    metrics["observed_ratio"] = round(ratio, 4)
    if ratio > ceiling_ratio:
        return Observation(
            "secondary_saturation",
            "failed",
            (
                f"{at_ceiling} reads sit on the -N {secondary_max} ceiling against {below} one "
                f"below it (ratio {ratio:.2f} > {ceiling_ratio}); raise -N and re-map"
            ),
            metrics,
        )
    return Observation(
        "secondary_saturation",
        "passed",
        (
            f"the secondary distribution decays into the -N {secondary_max} ceiling "
            f"({at_ceiling} at the ceiling, {below} below it)"
        ),
        metrics,
    )


# --- split reads across the retention boundary ------------------------------


def parse_sa_segments(tag: str) -> list[tuple[str, int]]:
    """Extract (contig, position) from an ``SA:Z`` tag value."""
    segments: list[tuple[str, int]] = []
    for entry in tag.split(";"):
        if not entry.strip():
            continue
        fields = entry.split(",")
        if len(fields) < 2:
            raise SmokeError(f"malformed SA segment: {entry!r}")
        try:
            position = int(fields[1])
        except ValueError as error:
            raise SmokeError(f"malformed SA segment position: {entry!r}") from error
        segments.append((fields[0], position))
    return segments


def load_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Read a BED into per-contig half-open intervals."""
    intervals: dict[str, list[tuple[int, int]]] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise SmokeError(f"malformed BED row: {line!r}")
        intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    if not intervals:
        raise SmokeError(f"BED contains no intervals: {path}")
    return intervals


def segment_inside(
    intervals: Mapping[str, Sequence[tuple[int, int]]], contig: str, position: int
) -> bool:
    """Report whether a 1-based alignment start falls inside a retained interval."""
    start = position - 1
    return any(low <= start < high for low, high in intervals.get(contig, ()))


def distance_to_intervals(
    intervals: Mapping[str, Sequence[tuple[int, int]]], contig: str, position: int
) -> int | None:
    """Distance from a 1-based position to the nearest retained interval edge.

    ``None`` when the contig carries no retained interval at all, which is a
    different situation from being just outside one.
    """
    spans = intervals.get(contig)
    if not spans:
        return None
    start = position - 1
    return min(
        0 if low <= start < high else (low - start if start < low else start - high + 1)
        for low, high in spans
    )


def classify_segments(
    segments: Iterable[tuple[str, int]],
    intervals: Mapping[str, Sequence[tuple[int, int]]],
    adjacency_window: int,
) -> tuple[bool, bool]:
    """Split a read's supplementary segments into boundary-cut and distant.

    A segment just outside a retained interval means the window is still cutting
    local structure — exactly what widening the BED was meant to stop. A segment
    on another contig, or megabases away, is ordinary paralogy or repeat content
    that no target-only mapping can retain, and counting it as a defect would
    make the check fail on healthy data.
    """
    truncated = False
    distant = False
    for contig, position in segments:
        if segment_inside(intervals, contig, position):
            continue
        distance = distance_to_intervals(intervals, contig, position)
        if distance is not None and distance <= adjacency_window:
            truncated = True
        else:
            distant = True
    return truncated, distant


def evaluate_split_reads(
    records: Iterable[tuple[str, str]],
    intervals: Mapping[str, Sequence[tuple[int, int]]],
    *,
    max_external_fraction: float = DEFAULT_MAX_EXTERNAL_FRACTION,
    adjacency_window: int = DEFAULT_ADJACENCY_WINDOW,
) -> Observation:
    """Measure how often the retention boundary cuts a split read's local structure.

    ``records`` are ``(qname, SA tag)`` pairs for alignments carrying one.
    """
    truncated: set[str] = set()
    distant: set[str] = set()
    seen: set[str] = set()
    for qname, tag in records:
        seen.add(qname)
        is_truncated, is_distant = classify_segments(
            parse_sa_segments(tag), intervals, adjacency_window
        )
        if is_truncated:
            truncated.add(qname)
        if is_distant:
            distant.add(qname)
    total = len(seen)
    metrics: dict[str, object] = {
        "reads_with_supplementary": total,
        "reads_cut_at_the_boundary": len(truncated),
        "reads_with_distant_segments": len(distant),
        "adjacency_window": adjacency_window,
        "max_external_fraction": max_external_fraction,
    }
    if total == 0:
        metrics["boundary_fraction"] = None
        metrics["distant_fraction"] = None
        return Observation(
            "sv_smoke",
            "not_run",
            "no supplementary alignments were retained, so split-read evidence is untested",
            metrics,
        )
    fraction = len(truncated) / total
    metrics["boundary_fraction"] = round(fraction, 4)
    metrics["distant_fraction"] = round(len(distant) / total, 4)
    tail = (
        f"{len(distant)} reads also align far from any target, which target-only mapping "
        "cannot retain and which is not counted here"
    )
    if fraction > max_external_fraction:
        return Observation(
            "sv_smoke",
            "failed",
            (
                f"{len(truncated)} of {total} split reads ({fraction:.1%}) have a segment within "
                f"{adjacency_window / 1e6:.1f} Mb outside a retained interval, above the "
                f"{max_external_fraction:.0%} allowance; the retention boundary is still cutting "
                f"local structure. {tail}"
            ),
            metrics,
        )
    return Observation(
        "sv_smoke",
        "passed",
        (
            f"{total - len(truncated)} of {total} split reads keep their local structure inside "
            f"the retained intervals ({fraction:.1%} cut at the boundary, allowance "
            f"{max_external_fraction:.0%}). {tail}"
        ),
        metrics,
    )


# --- assembly ---------------------------------------------------------------


def parse_assembly_info(text: str) -> list[dict[str, object]]:
    """Read Flye's ``assembly_info.txt`` into contig records."""
    contigs: list[dict[str, object]] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise SmokeError(f"malformed assembly_info row: {line!r}")
        contigs.append({"name": fields[0], "length": int(fields[1])})
    return contigs


def merge_intervals(intervals: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Collapse overlapping/adjacent (start, end) intervals into a minimal set."""
    merged: list[tuple[int, int]] = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def _excluded_overlap(start: int, end: int, excluded: Sequence[tuple[int, int]]) -> int:
    """Bases of (start, end) that fall inside any excluded interval."""
    return sum(max(0, min(end, hi) - max(start, lo)) for lo, hi in excluded)


def coverage_from_paf(
    lines: Iterable[str],
    region_span: int,
    exclude: Sequence[tuple[int, int]] = (),
) -> tuple[float, float, int]:
    """Return (breadth, identity, covered bases) of an assembly-to-region PAF.

    ``exclude`` removes a footprint (e.g. paralog loci this smoke does not
    evaluate) from the identity calculation only, proportionally within each
    alignment block, so a paralog-heavy sub-region cannot dominate a verdict
    about the rest of the assembly. Breadth still spans the whole region: a
    dropout inside excluded loci is still a real dropout.
    """
    excluded_merged = merge_intervals(exclude)
    spans: list[tuple[int, int]] = []
    matches = 0.0
    blocks = 0.0
    for line in lines:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise SmokeError(f"malformed PAF row: {line!r}")
        start, end = int(fields[7]), int(fields[8])
        spans.append((start, end))
        block_len = end - start
        excluded_bp = _excluded_overlap(start, end, excluded_merged) if block_len else 0
        keep_frac = 1 - (excluded_bp / block_len) if block_len else 0.0
        matches += int(fields[9]) * keep_frac
        blocks += int(fields[10]) * keep_frac
    covered = 0
    cursor = -1
    for start, end in sorted(spans):
        start = max(start, cursor)
        if end > start:
            covered += end - start
            cursor = end
    breadth = covered / region_span if region_span else 0.0
    identity = matches / blocks if blocks else 0.0
    return breadth, identity, covered


def evaluate_assembly(
    *,
    contigs: Sequence[Mapping[str, object]],
    region_span: int,
    breadth: float,
    identity: float,
    genes_expected: int,
    genes_recovered: int,
    thresholds: Mapping[str, float],
    excluded_genes: Sequence[str] = (),
) -> Observation:
    """Judge whether a local assembly reconstructs the region faithfully."""
    total_bp = sum(int(contig["length"]) for contig in contigs)
    recovery = genes_recovered / genes_expected if genes_expected else 0.0
    metrics: dict[str, object] = {
        "contigs": len(contigs),
        "total_bp": total_bp,
        "longest_bp": max((int(c["length"]) for c in contigs), default=0),
        "region_span": region_span,
        "span_ratio": round(total_bp / region_span, 4) if region_span else None,
        "breadth": round(breadth, 4),
        "identity": round(identity, 4),
        "genes_expected": genes_expected,
        "genes_recovered": genes_recovered,
        "gene_recovery": round(recovery, 4),
        "thresholds": dict(thresholds),
        # Named for audit, not just counted: docs/METHODS_EN.md excludes this
        # same paralog family from variant analysis for the same R9/Guppy5
        # reason, so this list is the concrete claim that only known-ambiguous
        # loci were left out, not an arbitrary subset.
        "excluded_paralog_genes": sorted(excluded_genes),
    }
    reasons: list[str] = []
    if not contigs:
        reasons.append("the assembly produced no contigs")
    if breadth < thresholds["min_breadth"]:
        reasons.append(
            f"breadth {breadth:.2f} is below {thresholds['min_breadth']:.2f}; "
            "the assembly does not span the region (collapse or dropout)"
        )
    if identity < thresholds["min_identity"]:
        reasons.append(
            f"identity {identity:.3f} is below {thresholds['min_identity']:.3f}; "
            "the assembly disagrees with the reference (chimerism or mis-join)"
        )
    if recovery < thresholds["min_gene_recovery"]:
        reasons.append(
            f"gene recovery {genes_recovered}/{genes_expected} is below "
            f"{thresholds['min_gene_recovery']:.0%}"
        )
    if reasons:
        return Observation("assembly_smoke", "failed", "; ".join(reasons), metrics)
    return Observation(
        "assembly_smoke",
        "passed",
        (
            f"{len(contigs)} contigs totalling {total_bp / 1e6:.2f} Mb reconstruct "
            f"{breadth:.0%} of the {region_span / 1e6:.2f} Mb region at {identity:.3f} identity, "
            f"recovering {genes_recovered}/{genes_expected} expected genes"
        ),
        metrics,
    )


# --- status document --------------------------------------------------------


def build_smoke_status(
    *,
    sample_id: str,
    bam_sha256: str,
    region: str,
    secondary_max: int,
    observations: Sequence[Observation],
    tool_versions: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Assemble the document the pilot gate reads, refusing partial evidence."""
    by_name = {observation.name: observation for observation in observations}
    missing = [name for name in REQUIRED_CHECKS if name not in by_name]
    if missing:
        raise SmokeError(f"smoke status is missing checks: {', '.join(missing)}")
    statuses = {by_name[name].status for name in REQUIRED_CHECKS}
    if "failed" in statuses:
        verdict = "failed"
    elif statuses - {"passed"}:
        verdict = "incomplete"
    else:
        verdict = "passed"
    return {
        "status": verdict,
        "sample_id": sample_id,
        "bam_sha256": bam_sha256,
        "region": region,
        "secondary_max": secondary_max,
        "evaluated_at_utc": datetime.now(timezone.utc).isoformat(),
        "tool_versions": dict(tool_versions or {}),
        "checks": [
            {
                "name": name,
                "status": by_name[name].status,
                "detail": by_name[name].detail,
                "metrics": by_name[name].metrics,
            }
            for name in REQUIRED_CHECKS
        ],
    }


# --- collection -------------------------------------------------------------


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def _run(command: Sequence[str], *, label: str) -> str:
    # A missing tool has to surface as a smoke failure, not an OSError that
    # escapes and leaves no status document behind.
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as error:
        raise SmokeError(f"{label} could not be run: {error}") from error
    if result.returncode != 0:
        raise SmokeError(f"{label} failed ({result.returncode}): {result.stderr.strip()[:400]}")
    return result.stdout


def collect_secondary(samtools: str, bam: Path, secondary_max: int, ceiling_ratio: float) -> Observation:
    output = _run([samtools, "view", "-f", "256", str(bam)], label="samtools view -f 256")
    qnames = (line.split("\t", 1)[0] for line in output.splitlines() if line)
    return evaluate_secondary_saturation(
        secondary_histogram(qnames), secondary_max, ceiling_ratio=ceiling_ratio
    )


def collect_split_reads(
    samtools: str,
    bam: Path,
    intervals: Mapping[str, Sequence[tuple[int, int]]],
    allowance: float,
    adjacency_window: int = DEFAULT_ADJACENCY_WINDOW,
) -> Observation:
    output = _run([samtools, "view", str(bam)], label="samtools view")
    records: list[tuple[str, str]] = []
    for line in output.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        for tag in fields[11:]:
            if tag.startswith("SA:Z:"):
                records.append((fields[0], tag[5:]))
                break
    return evaluate_split_reads(
        records,
        intervals,
        max_external_fraction=allowance,
        adjacency_window=adjacency_window,
    )


def genes_in_region(
    analytical_bed: Path, contig: str, start: int, end: int
) -> list[tuple[str, int, int]]:
    """List analytical targets wholly inside the assembled region."""
    genes: list[tuple[str, int, int]] = []
    for line in Path(analytical_bed).read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if fields[0] != contig:
            continue
        low, high = int(fields[1]), int(fields[2])
        if low >= start and high <= end:
            genes.append((fields[3] if len(fields) > 3 else f"{low}-{high}", low, high))
    return genes


def recovered_genes(
    paf_lines: Sequence[str], genes: Sequence[tuple[str, int, int]], minimum_overlap: float
) -> int:
    """Count genes whose interval is covered by assembly alignments."""
    spans = []
    for line in paf_lines:
        if not line.strip():
            continue
        fields = line.split("\t")
        spans.append((int(fields[7]), int(fields[8])))
    spans.sort()
    recovered = 0
    for _, low, high in genes:
        covered = 0
        cursor = low
        for start, end in spans:
            if end <= low or start >= high:
                continue
            start = max(start, cursor, low)
            end = min(end, high)
            if end > start:
                covered += end - start
                cursor = end
        if high > low and covered / (high - low) >= minimum_overlap:
            recovered += 1
    return recovered


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--mapping-bed", required=True, type=Path)
    parser.add_argument("--analytical-bed", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--region", required=True, help="contig:start-end of the assembled stretch")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--secondary-max", type=int, required=True)
    parser.add_argument("--ceiling-ratio", type=float, default=DEFAULT_CEILING_RATIO)
    parser.add_argument("--max-external-fraction", type=float, default=DEFAULT_MAX_EXTERNAL_FRACTION)
    parser.add_argument("--adjacency-window", type=int, default=DEFAULT_ADJACENCY_WINDOW)
    parser.add_argument("--min-breadth", type=float, default=DEFAULT_ASSEMBLY_THRESHOLDS["min_breadth"])
    parser.add_argument("--min-identity", type=float, default=DEFAULT_ASSEMBLY_THRESHOLDS["min_identity"])
    parser.add_argument(
        "--min-gene-recovery", type=float, default=DEFAULT_ASSEMBLY_THRESHOLDS["min_gene_recovery"]
    )
    parser.add_argument("--min-gene-overlap", type=float, default=0.90)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--flye", default="flye")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--work-dir", type=Path, default=None)
    return parser.parse_args(argv)


def _parse_region(region: str) -> tuple[str, int, int]:
    try:
        contig, coordinates = region.split(":", 1)
        start_text, end_text = coordinates.split("-", 1)
        return contig, int(start_text), int(end_text)
    except ValueError as error:
        raise SmokeError(f"malformed region: {region!r}") from error


def run_assembly(args: argparse.Namespace, work: Path) -> Observation:
    contig, start, end = _parse_region(args.region)
    region_span = end - start
    reads = work / "region_reads.fq"
    with reads.open("wb") as handle:
        view = subprocess.Popen(
            [args.samtools, "view", "-b", "-F", "0x900", str(args.bam), args.region],
            stdout=subprocess.PIPE,
        )
        fastq = subprocess.Popen(
            [args.samtools, "fastq", "-"], stdin=view.stdout, stdout=handle
        )
        assert view.stdout is not None
        view.stdout.close()
        if fastq.wait() != 0 or view.wait() != 0:
            raise SmokeError("extracting region reads failed")

    assembly_dir = work / "flye"
    _run(
        [
            args.flye,
            "--nano-raw",
            str(reads),
            "--genome-size",
            str(region_span),
            "--threads",
            str(args.threads),
            "--out-dir",
            str(assembly_dir),
        ],
        label="flye",
    )
    info = assembly_dir / "assembly_info.txt"
    fasta = assembly_dir / "assembly.fasta"
    if not info.is_file() or not fasta.is_file():
        raise SmokeError("flye produced no assembly")
    contigs = parse_assembly_info(info.read_text(encoding="utf-8"))

    region_fasta = work / "region.fa"
    region_fasta.write_text(
        _run([args.samtools, "faidx", str(args.reference), args.region], label="samtools faidx"),
        encoding="utf-8",
    )
    paf = _run(
        [args.minimap2, "-cx", "asm20", "-t", str(args.threads), str(region_fasta), str(fasta)],
        label="minimap2 asm20",
    ).splitlines()
    (work / "assembly_to_region.paf").write_text("\n".join(paf) + "\n", encoding="utf-8")

    genes = genes_in_region(args.analytical_bed, contig, start, end)
    # samtools faidx yields a region-local coordinate frame, so shift the genes
    # into it before asking whether the assembly covers them.
    local_genes = [(name, low - start, high - start) for name, low, high in genes]
    evaluated_genes = [g for g in local_genes if not g[0].startswith(EXCLUDED_GENE_PREFIX)]
    excluded_genes = [g for g in local_genes if g[0].startswith(EXCLUDED_GENE_PREFIX)]
    excluded_intervals = [(low, high) for _, low, high in excluded_genes]

    breadth, identity, _ = coverage_from_paf(paf, region_span, exclude=excluded_intervals)
    recovered = recovered_genes(paf, evaluated_genes, args.min_gene_overlap)
    return evaluate_assembly(
        contigs=contigs,
        region_span=region_span,
        breadth=breadth,
        identity=identity,
        genes_expected=len(evaluated_genes),
        genes_recovered=recovered,
        thresholds={
            "min_breadth": args.min_breadth,
            "min_identity": args.min_identity,
            "min_gene_recovery": args.min_gene_recovery,
        },
        excluded_genes=[name for name, _, _ in excluded_genes],
    )


def tool_versions(args: argparse.Namespace) -> dict[str, str]:
    """Record what produced the smoke verdict, including the smoke assembler.

    Flye sits outside the pinned mapping environment on purpose, so its version
    is not covered by the environment lock and has to be recorded here.
    """
    versions: dict[str, str] = {"python": sys.version.split()[0]}
    for label, command in (
        ("samtools", [args.samtools, "--version"]),
        ("minimap2", [args.minimap2, "--version"]),
        ("flye", [args.flye, "--version"]),
    ):
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            versions[label] = result.stdout.strip().splitlines()[0] if result.stdout.strip() else "unknown"
        except OSError:
            versions[label] = "unavailable"
    return versions


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    work = args.work_dir or args.output.parent / "work"
    work.mkdir(parents=True, exist_ok=True)
    try:
        intervals = load_intervals(args.mapping_bed)
        observations = [
            collect_secondary(args.samtools, args.bam, args.secondary_max, args.ceiling_ratio),
            collect_split_reads(
                args.samtools,
                args.bam,
                intervals,
                args.max_external_fraction,
                args.adjacency_window,
            ),
        ]
        # An assembler that dies leaves the operator with a shell exit code and
        # no record. Record the failure as the observation it is, so the gate
        # reports why the pilot is blocked instead of only that it is.
        try:
            observations.append(run_assembly(args, work))
        except SmokeError as error:
            observations.append(
                Observation(
                    "assembly_smoke",
                    "failed",
                    f"the assembly smoke could not be completed: {error}",
                    {"error": str(error)},
                )
            )
        status = build_smoke_status(
            sample_id=args.sample_id,
            bam_sha256=sha256_file(args.bam),
            region=args.region,
            secondary_max=args.secondary_max,
            observations=observations,
            tool_versions=tool_versions(args),
        )
    except SmokeError as error:
        sys.stderr.write(f"SMOKE FAILED: {error}\n")
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(f"{args.output.suffix}.partial")
    temporary.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    for check in status["checks"]:
        sys.stderr.write(f"  {check['name']:22s} {check['status']:10s} {check['detail']}\n")
    sys.stderr.write(f"\nsmoke verdict: {status['status']}\n")
    return 0 if status["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
