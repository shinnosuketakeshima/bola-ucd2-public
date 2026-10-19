"""Attribute assembled DQ sequences to reference positions and genotype DQ configurations.

The DQ/DR local assemblies (`06_local_assembly/dq_dr`) report which IPD-MHC
DQA/DQB/DRB3 sequence families each contig carries, but not *where* on the
haplotype they sit. A recovered DQA-family cluster therefore cannot, by
itself, show that a reference locus with no depth (e.g. DQA2) was recovered:
the cluster may be a different paralog. This stage places every family hit by
the reference sequence that flanks it on the same contig, and then genotypes
the configurations found on the contigs directly from reads in all animals.

1. **Anchors.** Each contig is compared with the ARS-UCD2.0 window. Tier 1
   uses `blastn` HSPs with identity >= `--min-identity` and length >=
   `--min-anchor`; tier 2 uses `minimap2 -x asm20` blocks >=
   `--min-divergent-block` whose matching bases over block length reach
   `--min-divergent-identity`, for contigs too divergent for tier 1. For each
   tier the heaviest colinear chain per contig (same orientation, increasing
   in both coordinates) is kept, so interspersed repeats, which give many
   short off-chain hits, do not enter it.
2. **Placement.** A family hit inside an anchor is projected to a reference
   coordinate (tier 1 first, then tier 2); a hit between anchors is
   bracketed by them. It is attributed to an annotated DQ feature only if the
   projection lies within `--feature-margin` of that feature.
3. **Events.** On the tier-1 chain, consecutive anchors whose reference gap
   exceeds their contig gap by more than `--min-event` mark a deletion
   relative to the reference; on the tier-2 chain, a contig gap exceeding the
   reference gap by that much marks an insertion. A deletion is
   read-supported when >= `--min-spanning` primary MAPQ >= 20 reads of the
   same animal, remapped to its own contigs, run continuously through the
   contig gap and `--flank` bp on both sides: no D/N/I longer than
   `--max-gap` in that window, none longer than `--boundary-max-gap` within
   `--boundary-margin` of the interval's ends, and aligned bases over at
   least `--min-aligned-fraction` of each flank and of the interval. End points
   enclosing the window are not enough, because a read can skip it with a
   large gap. Those reads are the assembly input, so this checks contig
   consistency; it is not independent validation.
4. **Configuration genotyping.** A three-sequence panel is built: the
   deletion junction from `--deletion-source`, the insertion and its DQA2
   hit from `--insertion-source`, and the reference over the same interval.
   Every animal's interval reads are mapped to the panel competitively, and
   reads running through each diagnostic point with `--flank` bp on both
   sides, by the same continuity rule, are counted. Counts show presence, not
   dosage. Every read whose end points enclose a window is written to
   `dq_support_reads.tsv` (read ID, strand, MAPQ, flank and interval aligned
   bases, longest gap, window-clipped CIGAR), with the endpoint-only count
   and a second `--max-gap-sensitivity` count kept beside the primary one.
5. **Origin of reference-DQA2 depth.** Primary MAPQ >= 20 reads aligned to
   reference DQA2 are remapped to the animal's own contigs and counted by
   whether they land within `--feature-margin` of a DQA hit placed at DQA2.
   This tells whether DQA2 depth comes from the DQA2-placed sequence or
   from reads of other contigs.
6. **DRB3 on contigs.** Each contig carrying a DRB3 hit is typed on the same
   shared exon-2 window as stage 03d (`drb3_exon2_typing`), by edit distance
   to the identical-window allele groups, so a configuration can be read
   together with the DRB3 group on the same contig.

Contigs are collapsed consensus sequences (Flye). Contig FASTA never leaves
`06_local_assembly`; this stage writes coordinates and counts only.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

CONTIG = "NC_037350.1"
DQ_FEATURES = ("BOLA-DQA2", "BOLA-DQB", "BOLA-DQA5", "BOLA-DQB1", "LOC100848815", "BOLA-DRB3")
FEATURE_LABELS = {
    "BOLA-DQA2": "DQA2",
    "BOLA-DQB": "DQB (reference BOLA-DQB)",
    "BOLA-DQA5": "DQA5",
    "BOLA-DQB1": "DQB1",
    "LOC100848815": "DQA1 (LOC100848815)",
    "BOLA-DRB3": "DRB3",
}
PANEL_PADDING = 10000
READ_FIELDS = ["sample_id", "analysis", "diagnostic", "target", "read_id", "strand", "mapq", "alignment_start",
               "alignment_end", "window_left", "window_right", "left_flank_aligned", "core_aligned",
               "right_flank_aligned", "core_length", "longest_deletion", "longest_insertion",
               "longest_gap_at_boundary", "window_identity", "window_cigar"]
REFERENCE_SPAN = re.compile(r"(\d+)([MDN=X])")


class AttributionError(ValueError):
    """Raised when inputs cannot support an attribution."""


@dataclass(frozen=True)
class Anchor:
    contig: str
    q_start: int  # 1-based, closed, on the contig
    q_end: int
    r_start: int  # 1-based, closed, chromosome coordinate (low < high)
    r_end: int
    strand: str  # "+" when contig and reference run in the same direction
    identity: float

    @property
    def length(self) -> int:
        return self.q_end - self.q_start + 1


def parse_blast_hsps(text: str, contig: str, offset: int, min_identity: float, min_length: int) -> list[Anchor]:
    """Parse `-outfmt "6 qstart qend sstart send pident length"` lines into anchors."""
    anchors = []
    for line in text.splitlines():
        if not line.strip():
            continue
        q_start, q_end, s_start, s_end, pident, length = line.split("\t")[:6]
        if float(pident) < min_identity or int(length) < min_length:
            continue
        s1, s2 = int(s_start) + offset, int(s_end) + offset
        anchors.append(
            Anchor(contig, int(q_start), int(q_end), min(s1, s2), max(s1, s2), "+" if s2 >= s1 else "-", float(pident))
        )
    return anchors


def parse_paf_blocks(text: str, offset: int, min_length: int, min_identity: float = 0.0) -> dict[str, list[Anchor]]:
    """Parse minimap2 PAF lines (0-based, half-open) into 1-based anchors per contig.

    The identity used here is matching bases over block length (gaps
    included), so a block stretched across a large deletion scores low and is
    dropped.
    """
    blocks: dict[str, list[Anchor]] = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        q_start, q_end = int(fields[2]) + 1, int(fields[3])
        if q_end - q_start + 1 < min_length:
            continue
        identity = 100.0 * int(fields[9]) / int(fields[10])
        if identity < min_identity:
            continue
        blocks.setdefault(fields[0], []).append(
            Anchor(fields[0], q_start, q_end, int(fields[7]) + 1 + offset, int(fields[8]) + offset, fields[4], identity)
        )
    return blocks


def chain_anchors(anchors: Sequence[Anchor], max_overlap: int = 500) -> list[Anchor]:
    """Return the heaviest colinear chain (by aligned contig length) of one contig's anchors."""
    best: list[Anchor] = []
    best_weight = 0
    for strand in ("+", "-"):
        items = sorted((a for a in anchors if a.strand == strand), key=lambda a: (a.q_start, a.q_end))
        weight = [a.length for a in items]
        previous = [-1] * len(items)
        for j, current in enumerate(items):
            for i in range(j):
                before = items[i]
                if before.q_end - max_overlap > current.q_start or before.q_start >= current.q_start:
                    continue
                if strand == "+":
                    colinear = before.r_end - max_overlap <= current.r_start and before.r_start < current.r_start
                else:
                    colinear = current.r_end - max_overlap <= before.r_start and current.r_start < before.r_start
                if colinear and weight[i] + current.length > weight[j]:
                    weight[j] = weight[i] + current.length
                    previous[j] = i
        if not items:
            continue
        end = max(range(len(items)), key=lambda k: weight[k])
        if weight[end] > best_weight:
            best_weight = weight[end]
            chain = []
            while end != -1:
                chain.append(items[end])
                end = previous[end]
            best = chain[::-1]
    return best


def project(chain: Sequence[Anchor], q_start: int, q_end: int) -> dict:
    """Place a contig interval on the reference using the chain that contains or brackets it."""
    middle = (q_start + q_end) // 2
    for anchor in chain:
        if anchor.q_start <= middle <= anchor.q_end:
            if anchor.strand == "+":
                position = anchor.r_start + (middle - anchor.q_start)
            else:
                position = anchor.r_end - (middle - anchor.q_start)
            return {"placement": "inside_anchor", "reference_position": position, "left": anchor, "right": anchor}
    left = [a for a in chain if a.q_end < middle]
    right = [a for a in chain if a.q_start > middle]
    if not left or not right:
        return {"placement": "beyond_chain", "reference_position": None,
                "left": left[-1] if left else None, "right": right[0] if right else None}
    return {"placement": "between_anchors", "reference_position": None, "left": left[-1], "right": right[0]}


def bracket(left: Anchor, right: Anchor) -> tuple[int, int]:
    """Reference interval between two consecutive chain anchors."""
    if left.strand == "+":
        return left.r_end, right.r_start
    return right.r_end, left.r_start


def attribute_feature(position: int | None, low: int | None, high: int | None,
                      features: dict[str, tuple[int, int]], margin: int) -> str:
    """Name the DQ feature a projection or bracket identifies, or say why none applies."""
    if position is not None:
        for name, (start, end) in features.items():
            if start - margin <= position <= end + margin:
                return FEATURE_LABELS[name]
        return "between_features"
    if low is None or high is None:
        return "unplaced"
    inside = [FEATURE_LABELS[n] for n, (s, e) in features.items() if s >= low and e <= high]
    return "bracketed:" + ("|".join(inside) if inside else "no_feature")


def find_events(chain: Sequence[Anchor], features: dict[str, tuple[int, int]], min_event: int,
                kind: str = "deletion") -> list[dict]:
    """Deletions (reference gap > contig gap) or insertions (contig gap > reference gap) between anchors."""
    events = []
    for left, right in zip(chain, chain[1:]):
        low, high = bracket(left, right)
        reference_gap = max(high - low - 1, 0)
        contig_gap = max(right.q_start - left.q_end - 1, 0)
        excess = reference_gap - contig_gap if kind == "deletion" else contig_gap - reference_gap
        if excess < min_event:
            continue
        removed = [FEATURE_LABELS[n] for n, (s, e) in features.items() if s > low and e < high]
        events.append(
            {
                "contig": left.contig,
                "event": kind,
                "contig_left": left.q_end,
                "contig_right": right.q_start,
                "contig_gap_bp": contig_gap,
                "reference_left": low,
                "reference_right": high,
                "reference_gap_bp": reference_gap,
                "net_bp": excess,
                "features_between": "|".join(removed) if removed else "none",
            }
        )
    return events


@dataclass(frozen=True)
class Alignment:
    read_id: str
    target: str
    start: int  # 1-based, closed, on the target
    end: int
    mapq: int
    flag: int
    cigar: str

    @property
    def strand(self) -> str:
        return "-" if self.flag & 0x10 else "+"

    def span(self) -> tuple[str, int, int, int, int]:
        return self.target, self.start, self.end, self.mapq, self.flag


def parse_sam_alignments(sam: str) -> list[Alignment]:
    """Every mapped record with its read ID and CIGAR, 1-based closed on the target."""
    alignments = []
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        fields = line.split("\t", 6)
        flag = int(fields[1])
        if flag & 0x4:
            continue
        start = int(fields[3])
        length = sum(int(n) for n, _ in REFERENCE_SPAN.findall(fields[5]))
        alignments.append(Alignment(fields[0], fields[2], start, start + length - 1, int(fields[4]), flag, fields[5]))
    return alignments


def parse_sam_spans(sam: str) -> list[tuple[str, int, int, int, int]]:
    """(target, start, end, mapq, flag) of every mapped record, 1-based closed on the target."""
    return [alignment.span() for alignment in parse_sam_alignments(sam)]


def window_profile(alignment: Alignment, left: int, right: int, flank: int, boundary_margin: int = 200) -> dict:
    """How the query actually covers [left - flank, right + flank] of the target.

    Aligned bases (M/=/X) are counted separately in the left flank, the
    diagnostic interval [left, right] and the right flank. The longest D/N
    overlapping the window and the longest I placed inside it are recorded,
    so an alignment whose end points enclose the window but which skips it
    with a large gap is visible as such; so is the longest D/N/I within
    `boundary_margin` bp of a diagnostic boundary (`left` and `right`), where
    a gap could carry the read across the boundary without query bases on
    it. The CIGAR is clipped to the window.
    """
    low, high = left - flank, right + flank
    boundaries = {left, right}
    longest_at_boundary = 0

    def near_boundary(first: int, last: int) -> bool:
        return any(first <= b + boundary_margin and last >= b - boundary_margin for b in boundaries)
    regions = ((low, left - 1), (left, right), (right + 1, high))
    aligned = [0, 0, 0]
    longest_deletion = longest_insertion = 0
    matches = columns = 0
    clipped: list[list] = []

    def keep(length: int, op: str) -> None:
        if clipped and clipped[-1][1] == op:
            clipped[-1][0] += length
        else:
            clipped.append([length, op])

    position = alignment.start
    for count, op in re.findall(r"(\d+)([MIDNSHP=X])", alignment.cigar):
        n = int(count)
        if op in "M=X":
            for k, (r_low, r_high) in enumerate(regions):
                aligned[k] += max(0, min(position + n - 1, r_high) - max(position, r_low) + 1)
            inside = max(0, min(position + n - 1, high) - max(position, low) + 1)
            if inside:
                keep(inside, "M")
                columns += inside
                matches += inside if op == "=" else 0
            position += n
        elif op in "DN":
            inside = max(0, min(position + n - 1, high) - max(position, low) + 1)
            if inside:
                longest_deletion = max(longest_deletion, n)
                keep(inside, op)
                columns += inside
            if near_boundary(position, position + n - 1):
                longest_at_boundary = max(longest_at_boundary, n)
            position += n
        elif op == "I" and low < position <= high:
            longest_insertion = max(longest_insertion, n)
            columns += n
            if near_boundary(position, position):
                longest_at_boundary = max(longest_at_boundary, n)
            keep(n, "I")
    return {
        "left_flank_aligned": aligned[0],
        "core_aligned": aligned[1],
        "right_flank_aligned": aligned[2],
        "core_length": right - left + 1,
        "longest_deletion": longest_deletion,
        "longest_insertion": longest_insertion,
        "longest_gap_at_boundary": longest_at_boundary,
        # matching bases over alignment columns; only defined for =/X (--eqx) CIGARs
        "window_identity": round(matches / columns, 4) if columns and re.search(r"\d[=X]", alignment.cigar) else "",
        "window_cigar": "".join(f"{n}{op}" for n, op in clipped),
    }


def encloses(alignment: Alignment, target: str, left: int, right: int, flank: int, min_mapq: int = 20) -> bool:
    """Primary, MAPQ >= min_mapq, and end points enclosing [left - flank, right + flank] (the pre-CIGAR rule)."""
    return (alignment.target == target and not alignment.flag & 0x904 and alignment.mapq >= min_mapq
            and alignment.start <= left - flank and alignment.end >= right + flank)


def continuous(profile: dict, flank: int, max_gap: int, min_aligned_fraction: float,
               boundary_max_gap: int = 50) -> bool:
    """True when the query itself runs through the window: no gap longer than max_gap in it,
    none longer than boundary_max_gap near a diagnostic boundary, and at least
    min_aligned_fraction of each flank and of the diagnostic interval (for a point, the
    point itself) covered by aligned bases. Small ONT indels are tolerated."""
    if profile["longest_deletion"] > max_gap or profile["longest_insertion"] > max_gap:
        return False
    if profile["longest_gap_at_boundary"] > min(boundary_max_gap, max_gap):
        return False
    if min(profile["left_flank_aligned"], profile["right_flank_aligned"]) < min_aligned_fraction * flank:
        return False
    return profile["core_aligned"] >= min_aligned_fraction * profile["core_length"]


def supporting_reads(alignments: Iterable[Alignment], target: str, left: int, right: int, flank: int,
                     max_gaps: Sequence[int], min_aligned_fraction: float, min_mapq: int = 20,
                     boundary_max_gap: int = 50, boundary_margin: int = 200) -> list[dict]:
    """One row per alignment that encloses the window, with whether it is continuous at each max_gap."""
    rows = []
    for alignment in alignments:
        if not encloses(alignment, target, left, right, flank, min_mapq):
            continue
        profile = window_profile(alignment, left, right, flank, boundary_margin)
        rows.append({
            "read_id": alignment.read_id, "target": target, "strand": alignment.strand, "mapq": alignment.mapq,
            "alignment_start": alignment.start, "alignment_end": alignment.end,
            "window_left": left - flank, "window_right": right + flank, **profile,
            **{f"continuous_gap{gap}": continuous(profile, flank, gap, min_aligned_fraction, boundary_max_gap)
               for gap in max_gaps},
        })
    return rows


def count_spanning(alignments: Iterable[Alignment], target: str, left: int, right: int, flank: int,
                   min_mapq: int = 20, max_gap: int = 200, min_aligned_fraction: float = 0.8,
                   boundary_max_gap: int = 50, boundary_margin: int = 200) -> int:
    """Count primary MAPQ >= min_mapq alignments whose query runs continuously through
    [left - flank, right + flank] of a target (see `continuous`)."""
    return sum(
        1 for row in supporting_reads(alignments, target, left, right, flank, (max_gap,), min_aligned_fraction,
                                      min_mapq, boundary_max_gap, boundary_margin)
        if row[f"continuous_gap{max_gap}"]
    )


def load_features(gff: Path, names: Sequence[str] = DQ_FEATURES) -> dict[str, tuple[int, int]]:
    features: dict[str, tuple[int, int]] = {}
    with gff.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != CONTIG or fields[2] != "gene":
                continue
            attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            if attributes.get("Name") in names:
                if attributes["Name"] in features:
                    raise AttributionError(f"duplicate gene record for {attributes['Name']}")
                features[attributes["Name"]] = (int(fields[3]), int(fields[4]))
    missing = [n for n in names if n not in features]
    if missing:
        raise AttributionError(f"DQ features missing from annotation: {', '.join(missing)}")
    return features


def load_family_hits(copies: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = {}
    with copies.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hits.setdefault(row["sample_id"], []).append(row)
    if not hits:
        raise AttributionError(f"no family hits in {copies}")
    return hits


def load_pcr_sbt(registry: Path) -> dict[str, tuple[str, str]]:
    with registry.open(encoding="utf-8", newline="") as handle:
        return {
            row["sample_id"]: (row["pcr_sbt_drb3_hap1"], row["pcr_sbt_drb3_hap2"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("analysis_status") == "analytical"
        }


def parse_source(value: str) -> tuple[str, str]:
    sample, _, contig = value.partition(":")
    if not sample or not contig:
        raise AttributionError(f"source must be SAMPLE:CONTIG, got {value!r}")
    return sample, contig


def _run(command: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(list(command), check=True, capture_output=True, text=True, **kwargs)


def fasta_sequence(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def contig_names(fasta: Path) -> list[str]:
    return [line[1:].split()[0] for line in fasta.read_text(encoding="utf-8").splitlines() if line.startswith(">")]


def tier1_chains(assembly: Path, window_fasta: Path, offset: int, workdir: Path, args) -> dict[str, list[Anchor]]:
    chains = {}
    for contig in contig_names(assembly):
        query = workdir / f"{contig}.fa"
        query.write_text(_run([args.samtools, "faidx", str(assembly), contig]).stdout, encoding="utf-8")
        output = _run(
            [args.blastn, "-query", str(query), "-subject", str(window_fasta), "-evalue", "1e-50",
             "-max_hsps", "1000", "-outfmt", "6 qstart qend sstart send pident length"]
        ).stdout
        chain = chain_anchors(parse_blast_hsps(output, contig, offset, args.min_identity, args.min_anchor))
        if chain:
            chains[contig] = chain
    return chains


def tier2_chains(assembly: Path, window_fasta: Path, offset: int, args) -> dict[str, list[Anchor]]:
    paf = _run([args.minimap2, "-x", "asm20", "-c", "--secondary=no", str(window_fasta), str(assembly)]).stdout
    return {
        contig: chain
        for contig, blocks in parse_paf_blocks(paf, offset, args.min_divergent_block, args.min_divergent_identity).items()
        if (chain := chain_anchors(blocks))
    }


def window_reads(bam: Path, region: str, destination: Path, samtools: str) -> Path:
    view = subprocess.run([samtools, "view", "-b", "-F", "0x904", str(bam), region], check=True, capture_output=True)
    fasta = subprocess.run([samtools, "fasta", "-"], input=view.stdout, check=True, capture_output=True)
    destination.write_bytes(fasta.stdout)
    return destination


def window_reads_mapq(bam: Path, region: str, destination: Path, samtools: str, min_mapq: int = 20) -> Path:
    view = subprocess.run([samtools, "view", "-b", "-F", "0x904", "-q", str(min_mapq), str(bam), region],
                          check=True, capture_output=True)
    fasta = subprocess.run([samtools, "fasta", "-"], input=view.stdout, check=True, capture_output=True)
    destination.write_bytes(fasta.stdout)
    return destination


def map_reads(target: Path, reads: Path, args) -> list[Alignment]:
    sam = _run([args.minimap2, "-a", "-x", "map-ont", "--secondary=no", "-t", str(args.threads), str(target), str(reads)]).stdout
    return parse_sam_alignments(sam)


def place_hit(hit: dict, tier1: list[Anchor], tier2: list[Anchor], features, margin: int) -> dict:
    q_start, q_end = int(hit["ref_start"]), int(hit["ref_end"])
    for tier, chain in ((1, tier1), (2, tier2)):
        placed = project(chain, q_start, q_end)
        if placed["placement"] == "inside_anchor":
            return {
                "placement": "inside_anchor", "tier": tier, "reference_position": placed["reference_position"],
                "bracket_left": "", "bracket_right": "",
                "attributed_locus": attribute_feature(placed["reference_position"], None, None, features, margin),
            }
    for tier, chain in ((1, tier1), (2, tier2)):
        placed = project(chain, q_start, q_end)
        if placed["placement"] == "between_anchors":
            low, high = bracket(placed["left"], placed["right"])
            return {
                "placement": "between_anchors", "tier": tier, "reference_position": "",
                "bracket_left": low, "bracket_right": high,
                "attributed_locus": attribute_feature(None, low, high, features, margin),
            }
    return {"placement": "unplaced", "tier": "", "reference_position": "", "bracket_left": "", "bracket_right": "",
            "attributed_locus": "unplaced"}


def dqa2_read_origin(alignments: Iterable[tuple[str, int, int, int, int]], dqa2_hits: Sequence[dict],
                     margin: int, min_mapq: int = 20) -> dict:
    """Split reads remapped to own contigs by whether they overlap a DQA2-placed DQA hit."""
    on_hit = elsewhere = 0
    for contig, start, end, mapq, flag in alignments:
        if flag & 0x904 or mapq < min_mapq:
            continue
        if any(contig == h["contig"] and start <= int(h["contig_end"]) + margin and end >= int(h["contig_start"]) - margin
               for h in dqa2_hits):
            on_hit += 1
        else:
            elsewhere += 1
    return {"reads_on_dqa2_placed_hit": on_hit, "reads_elsewhere": elsewhere}


def build_panel(args, features, deletion: dict, insertion: dict, dqa2_hit: dict, workdir: Path) -> tuple[Path, list[dict]]:
    """Write the three-sequence configuration panel and return its diagnostic points."""
    del_sample, del_contig = parse_source(args.deletion_source)
    ins_sample, ins_contig = parse_source(args.insertion_source)
    dqa2_start, dqa2_end = features["BOLA-DQA2"]

    del_lo = max(deletion["contig_left"] - PANEL_PADDING, 1)
    del_hi = deletion["contig_right"] + PANEL_PADDING
    ins_lo = max(min(insertion["contig_left"], int(dqa2_hit["ref_start"])) - PANEL_PADDING, 1)
    ins_hi = max(insertion["contig_right"], int(dqa2_hit["ref_end"])) + PANEL_PADDING
    ref_lo = min(insertion["reference_left"], deletion["reference_left"]) - PANEL_PADDING
    ref_hi = dqa2_end + PANEL_PADDING

    pieces = [
        ("deletion_config", _run([args.samtools, "faidx", str(args.assembly_dir / del_sample / "assembly.fasta"),
                                  f"{del_contig}:{del_lo}-{del_hi}"]).stdout),
        ("insertion_config", _run([args.samtools, "faidx", str(args.assembly_dir / ins_sample / "assembly.fasta"),
                                   f"{ins_contig}:{ins_lo}-{ins_hi}"]).stdout),
        ("reference_config", _run([args.samtools, "faidx", str(args.reference), f"{CONTIG}:{ref_lo}-{ref_hi}"]).stdout),
    ]
    panel = workdir / "configuration_panel.fa"
    panel.write_text(
        "".join(f">{name}\n{fasta_sequence(text)}\n" for name, text in pieces), encoding="utf-8"
    )
    diagnostics = [
        {"diagnostic": "deletion_junction", "target": "deletion_config",
         "left": deletion["contig_left"] - del_lo + 1, "right": deletion["contig_right"] - del_lo + 1,
         "meaning": "reads joining the reference sequence on both sides of the absent DQA2-DQB-DQA5 segment"},
        {"diagnostic": "insertion_left_edge", "target": "insertion_config",
         "left": insertion["contig_left"] - ins_lo + 1, "right": insertion["contig_left"] - ins_lo + 1,
         "meaning": "reads crossing from reference-colinear sequence into the non-reference segment"},
        {"diagnostic": "insertion_whole_segment", "target": "insertion_config",
         "left": insertion["contig_left"] - ins_lo + 1, "right": insertion["contig_right"] - ins_lo + 1,
         "meaning": "reads running through the whole non-reference segment and into sequence on both sides"},
        {"diagnostic": "insertion_right_edge", "target": "insertion_config",
         "left": insertion["contig_right"] - ins_lo + 1, "right": insertion["contig_right"] - ins_lo + 1,
         "meaning": "reads crossing from the non-reference segment into DQA2-colinear sequence"},
        {"diagnostic": "insertion_config_DQA2", "target": "insertion_config",
         "left": int(dqa2_hit["ref_start"]) - ins_lo + 1, "right": int(dqa2_hit["ref_end"]) - ins_lo + 1,
         "meaning": f"reads covering the DQA2-placed {dqa2_hit['best_allele']} hit on the insertion contig"},
        {"diagnostic": "reference_insertion_point", "target": "reference_config",
         "left": insertion["reference_left"] - ref_lo + 1, "right": insertion["reference_right"] - ref_lo + 1,
         "meaning": "reads continuous with the reference where the insertion contig carries extra sequence"},
        {"diagnostic": "reference_deletion_point", "target": "reference_config",
         "left": deletion["reference_left"] - ref_lo + 1, "right": deletion["reference_left"] - ref_lo + 1,
         "meaning": "reads continuous with the reference into DQA2 at the deletion's left edge"},
        {"diagnostic": "reference_DQA2", "target": "reference_config",
         "left": dqa2_start - ref_lo + 1, "right": dqa2_end - ref_lo + 1,
         "meaning": "reads covering reference DQA2 itself"},
    ]
    return panel, diagnostics


def drb3_typing_context(args):
    """Panel groups and the shared exon-2 window used by stage 03d."""
    from scripts.codex_reanalysis import drb3_exon2_typing as typing

    contig, exon_start, exon_end, _ = typing.exon_from_gff(args.gff)
    ref_exon = typing.reference_exon_sequence(args.reference, contig, exon_start, exon_end, args.samtools)
    panel = typing.cut_to_exon(typing.drb3_panel(args.drb3_panel), ref_exon)
    w_start, w_end, _ = typing.panel_window(panel, ref_exon)
    window = ref_exon[w_start:w_end]
    trimmed = {name: typing.trim_to(window, seq)[0] for name, seq in panel.items()}
    seqs, names = typing.group_identical(trimmed)
    return typing, window, seqs, names


def type_contig_drb3(typing, window: str, seqs, names, contig_sequence: str) -> dict:
    segment, anchor_edits = typing.trim_to(window, contig_sequence)
    if anchor_edits > typing.MAX_ANCHOR_ERROR * len(window):
        return {"drb3_group": "", "drb3_edits": "", "drb3_window_found": False}
    distances = typing.distance_matrix([segment], seqs)[0]
    best = min(distances)
    groups = [names[i] for i, d in enumerate(distances) if d == best]
    return {
        "drb3_group": " | ".join("/".join(g) for g in groups),
        "drb3_edits": best,
        "drb3_window_found": True,
    }


def write_tsv(path: Path, rows: list[dict], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assembly-dir", type=Path, required=True, help="06_local_assembly/dq_dr")
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--drb3-panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--deletion-source", required=True, help="SAMPLE:CONTIG carrying the deletion junction")
    parser.add_argument("--insertion-source", required=True, help="SAMPLE:CONTIG carrying the insertion and DQA2")
    parser.add_argument("--window", default=f"{CONTIG}:25400001-26000000")
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--min-anchor", type=int, default=2000)
    parser.add_argument("--min-divergent-block", type=int, default=5000)
    parser.add_argument("--min-divergent-identity", type=float, default=70.0,
                        help="tier 2: matching bases over block length, percent")
    parser.add_argument("--min-event", type=int, default=10000)
    parser.add_argument("--feature-margin", type=int, default=2000)
    parser.add_argument("--flank", type=int, default=2000)
    parser.add_argument("--min-spanning", type=int, default=3)
    parser.add_argument("--max-gap", type=int, default=200,
                        help="longest D/N/I tolerated inside a supporting read's window, bp")
    parser.add_argument("--boundary-max-gap", type=int, default=50,
                        help="longest D/N/I tolerated within --boundary-margin of a diagnostic boundary, bp")
    parser.add_argument("--boundary-margin", type=int, default=200)
    parser.add_argument("--max-gap-sensitivity", type=int, default=50,
                        help="a second (strict) max-gap value whose counts are reported alongside")
    parser.add_argument("--min-aligned-fraction", type=float, default=0.8,
                        help="aligned bases over length required in each flank and in the diagnostic interval")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--blastn", default="blastn")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    features = load_features(args.gff)
    copies = args.assembly_dir / "dq_dr_copies.tsv"
    hits = load_family_hits(copies)
    pcr_sbt = load_pcr_sbt(args.registry)
    samples = sorted(hits)
    if len(samples) != 13 or set(samples) != set(pcr_sbt):
        raise AttributionError(f"expected the 13 analytical animals, found {len(samples)} with family hits")
    offset = int(args.window.split(":")[1].split("-")[0]) - 1
    typing, drb3_window, drb3_seqs, drb3_names = drb3_typing_context(args)

    gaps = (args.max_gap, args.max_gap_sensitivity)
    primary, sensitivity = f"continuous_gap{args.max_gap}", f"continuous_gap{args.max_gap_sensitivity}"
    hit_rows: list[dict] = []
    event_rows: list[dict] = []
    read_rows: list[dict] = []
    comparison_rows: list[dict] = []
    reads_by_sample: dict[str, Path] = {}
    origin_rows: list[dict] = []
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        window_fasta = workdir / "window.fa"
        window_fasta.write_text(_run([args.samtools, "faidx", str(args.reference), args.window]).stdout, encoding="utf-8")
        for sample in samples:
            assembly = args.assembly_dir / sample / "assembly.fasta"
            bam = args.bam_dir / f"{sample}.targets.bam"
            if not assembly.is_file() or not bam.is_file():
                raise AttributionError(f"missing assembly or BAM for {sample}")
            sample_dir = workdir / sample
            sample_dir.mkdir()
            tier1 = tier1_chains(assembly, window_fasta, offset, sample_dir, args)
            tier2 = tier2_chains(assembly, window_fasta, offset, args)
            reads_by_sample[sample] = window_reads(bam, args.window, sample_dir / "reads.fa", args.samtools)
            dqa2_start, dqa2_end = features["BOLA-DQA2"]
            dqa2_reads = window_reads_mapq(bam, f"{CONTIG}:{dqa2_start}-{dqa2_end}", sample_dir / "dqa2_reads.fa",
                                           args.samtools)
            alignments = map_reads(assembly, reads_by_sample[sample], args)
            for hit in hits[sample]:
                row = {
                    "sample_id": sample,
                    "family": hit["family"],
                    "contig": hit["contig"],
                    "contig_start": hit["ref_start"],
                    "contig_end": hit["ref_end"],
                    "best_allele": hit["best_allele"],
                    **place_hit(hit, tier1.get(hit["contig"], []), tier2.get(hit["contig"], []), features,
                                args.feature_margin),
                    "drb3_group": "",
                    "drb3_edits": "",
                    "pcr_sbt": "/".join(pcr_sbt[sample]),
                }
                if hit["family"] == "DRB3":
                    sequence = fasta_sequence(_run([args.samtools, "faidx", str(assembly), hit["contig"]]).stdout)
                    typed = type_contig_drb3(typing, drb3_window, drb3_seqs, drb3_names, sequence)
                    row["drb3_group"], row["drb3_edits"] = typed["drb3_group"], typed["drb3_edits"]
                hit_rows.append(row)
            sample_hits = [r for r in hit_rows if r["sample_id"] == sample]
            dqa2_hits = [r for r in sample_hits if r["family"] == "DQA" and r["attributed_locus"] == "DQA2"]
            reference_dqa2_reads = sum(1 for line in dqa2_reads.read_text(encoding="utf-8").splitlines() if line.startswith(">"))
            origin = (dqa2_read_origin([a.span() for a in map_reads(assembly, dqa2_reads, args)], dqa2_hits,
                                       args.feature_margin)
                      if reference_dqa2_reads else {"reads_on_dqa2_placed_hit": 0, "reads_elsewhere": 0})
            origin_rows.append({"sample_id": sample, "reference_dqa2_reads": reference_dqa2_reads,
                                "dqa2_placed_hits": len(dqa2_hits), **origin})
            for tier, chains in ((1, tier1), (2, tier2)):
                for kind in ("deletion", "insertion"):
                    for contig, chain in chains.items():
                        for event in find_events(chain, features, args.min_event, kind):
                            reads = supporting_reads(alignments, contig, event["contig_left"], event["contig_right"],
                                                     args.flank, gaps, args.min_aligned_fraction, 20,
                                                     args.boundary_max_gap, args.boundary_margin)
                            spanning = sum(r[primary] for r in reads)
                            event_rows.append(
                                {"sample_id": sample, "tier": tier, **event, "spanning_reads": spanning,
                                 "read_supported": spanning >= args.min_spanning,
                                 "spanning_reads_endpoint_only": len(reads),
                                 f"spanning_reads_gap{args.max_gap_sensitivity}": sum(r[sensitivity] for r in reads)}
                            )
                            label = f"tier{tier}_{kind}:{event['contig_left']}-{event['contig_right']}"
                            read_rows.extend({"sample_id": sample, "analysis": "own_contig", "diagnostic": label, **r}
                                             for r in reads)

        deletion = _select_event(event_rows, args.deletion_source, "deletion", "DQA2", tier=1)
        insertion = _select_event(event_rows, args.insertion_source, "insertion", None, tier=2,
                                  before=features["BOLA-DQA2"][0])
        dqa2_hit = _select_dqa2_hit(hit_rows, args.insertion_source)
        panel, diagnostics = build_panel(args, features, deletion, insertion, dqa2_hit, workdir)
        genotype_rows = []
        for sample in samples:
            alignments = map_reads(panel, reads_by_sample[sample], args)
            row = {"sample_id": sample}
            for point in diagnostics:
                reads = supporting_reads(alignments, point["target"], point["left"], point["right"], args.flank,
                                         gaps, args.min_aligned_fraction, 20, args.boundary_max_gap,
                                         args.boundary_margin)
                row[point["diagnostic"]] = sum(r[primary] for r in reads)
                comparison_rows.append({
                    "sample_id": sample, "diagnostic": point["diagnostic"], "endpoint_only": len(reads),
                    f"continuous_gap{args.max_gap}": row[point["diagnostic"]],
                    f"continuous_gap{args.max_gap_sensitivity}": sum(r[sensitivity] for r in reads),
                })
                read_rows.extend({"sample_id": sample, "analysis": "panel", "diagnostic": point["diagnostic"], **r}
                                 for r in reads)
            genotype_rows.append(row)

    hit_fields = ["sample_id", "family", "contig", "contig_start", "contig_end", "best_allele", "placement", "tier",
                  "reference_position", "bracket_left", "bracket_right", "attributed_locus", "drb3_group",
                  "drb3_edits", "pcr_sbt"]
    event_fields = ["sample_id", "tier", "contig", "event", "contig_left", "contig_right", "contig_gap_bp", "reference_left",
                    "reference_right", "reference_gap_bp", "net_bp", "features_between", "spanning_reads",
                    "read_supported", "spanning_reads_endpoint_only", f"spanning_reads_gap{args.max_gap_sensitivity}"]
    write_tsv(args.output / "dq_family_attribution.tsv", hit_rows, hit_fields)
    write_tsv(args.output / "dq_contig_events.tsv", event_rows, event_fields)
    write_tsv(args.output / "dq_configuration_reads.tsv", genotype_rows,
              ["sample_id", *[p["diagnostic"] for p in diagnostics]])
    write_tsv(args.output / "dq_dqa2_read_origin.tsv", origin_rows,
              ["sample_id", "reference_dqa2_reads", "dqa2_placed_hits", "reads_on_dqa2_placed_hit", "reads_elsewhere"])
    write_tsv(args.output / "dq_configuration_support_comparison.tsv", comparison_rows,
              ["sample_id", "diagnostic", "endpoint_only", primary, sensitivity])
    write_tsv(args.output / "dq_support_reads.tsv", read_rows, READ_FIELDS + [primary, sensitivity])
    write_tsv(args.output / "dq_configuration_diagnostics.tsv", diagnostics,
              ["diagnostic", "target", "left", "right", "meaning"])
    parameters = {
        "window": args.window,
        "tier1": f"blastn HSP identity >= {args.min_identity}%, length >= {args.min_anchor} bp",
        "tier2": f"minimap2 -x asm20 -c block >= {args.min_divergent_block} bp, matches/block length "
                 f">= {args.min_divergent_identity}%",
        "min_event_bp": args.min_event,
        "feature_margin_bp": args.feature_margin,
        "spanning_flank_bp": args.flank,
        "min_spanning_reads": args.min_spanning,
        "support_rule": "primary MAPQ >= 20 alignment enclosing [left - flank, right + flank] with no D/N/I "
                        f"> {args.max_gap} bp in that window, none > {args.boundary_max_gap} bp within "
                        f"{args.boundary_margin} bp of a diagnostic boundary, and aligned (M/=/X) bases over >= "
                        f"{args.min_aligned_fraction:g} of each flank and of the diagnostic interval",
        "max_gap_bp": args.max_gap,
        "boundary_max_gap_bp": args.boundary_max_gap,
        "boundary_margin_bp": args.boundary_margin,
        "max_gap_sensitivity_bp": args.max_gap_sensitivity,
        "min_aligned_fraction": args.min_aligned_fraction,
        "read_filter": "target BAM primary reads in window (-F 0x904), remapped with minimap2 -x map-ont "
                       "--secondary=no; primary MAPQ >= 20",
        "deletion_source": args.deletion_source,
        "insertion_source": args.insertion_source,
        "drb3_window_nt": len(drb3_window),
    }
    status = {
        "status": "passed",
        "samples": len(samples),
        "family_hits": len(hit_rows),
        "events": len(event_rows),
        "read_supported_deletions": sum(
            1 for r in event_rows if r["event"] == "deletion" and r["tier"] == 1 and r["read_supported"]
        ),
        "selected_deletion": deletion,
        "selected_insertion": insertion,
        "features": {FEATURE_LABELS[k]: list(v) for k, v in features.items()},
        "parameters": parameters,
    }
    (args.output / "dq_attribution_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )
    write_provenance(
        args.output,
        {
            "versions": tool_versions(args),
            "inputs": {
                "family_hits": str(copies),
                "family_hits_sha256": sha256_file(copies),
                "reference_sha256": sha256_file(args.reference),
                "registry_sha256": sha256_file(args.registry),
                "drb3_panel_sha256": sha256_file(args.drb3_panel),
                "gff": str(args.gff),
            },
            "parameters": parameters,
        },
    )
    print(f"DQ attribution written: {args.output}")
    return 0


def tool_versions(args) -> dict:
    first_line = lambda command: _run(command).stdout.strip().splitlines()[0]
    return {
        "python": sys.version.split()[0],
        "minimap2": first_line([args.minimap2, "--version"]),
        "samtools": first_line([args.samtools, "--version"]),
        "blastn": first_line([args.blastn, "-version"]),
    }


def _select_event(rows: list[dict], source: str, kind: str, required_feature: str | None,
                  tier: int, before: int | None = None) -> dict:
    sample, contig = parse_source(source)
    for row in rows:
        if (row["sample_id"], row["contig"], row["event"], row["tier"]) != (sample, contig, kind, tier):
            continue
        if before is not None and row["reference_right"] > before:
            continue
        if required_feature and required_feature not in row["features_between"].split("|"):
            continue
        if kind == "deletion" and not row["read_supported"]:
            raise AttributionError(f"{source}: the {kind} is not read-supported")
        return row
    raise AttributionError(f"{source}: no {kind} event{' removing ' + required_feature if required_feature else ''}")


def _select_dqa2_hit(rows: list[dict], source: str) -> dict:
    sample, contig = parse_source(source)
    for row in rows:
        if (row["sample_id"], row["contig"], row["family"]) == (sample, contig, "DQA") and row["attributed_locus"] == "DQA2":
            return {"ref_start": row["contig_start"], "ref_end": row["contig_end"], "best_allele": row["best_allele"]}
    raise AttributionError(f"{source}: no DQA hit placed at DQA2")


if __name__ == "__main__":
    raise SystemExit(main())
