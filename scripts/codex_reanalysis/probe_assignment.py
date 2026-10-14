"""Assign IPD-MHC probes to local-assembly contigs.

IPD-panel hits count as locus evidence only at ≥95% identity, ≥80% query
coverage, and MAPQ ≥20 (`docs/superpowers/specs/2026-07-31-ucd2-paper-reanalysis-design.md`
and the project mapping-evidence rule). Hits below any of those stay
``unresolved`` and never become a locus call.

Cluster counts on haplotype-kept Flye contigs are recovery summaries, not
validated diploid copy numbers: ``--keep-haplotypes`` can split one gene
across contigs, and collapse can merge two alleles onto one contig.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Sequence


HIGH_CONFIDENCE = "high_confidence"
UNRESOLVED = "unresolved"

RECOVERED = "recovered"
COLLAPSED = "collapsed"
HAPLOTYPE_SPLIT = "haplotype_split"
DUPLICATED_ON_CONTIG = "duplicated_on_contig"
NOT_RECOVERED = "not_recovered"

MIN_IDENTITY = 0.95
MIN_COVERAGE = 0.80
MIN_MAPQ = 20
CLUSTER_GAP = 8000
MAX_SPAN = 12000

CLASSICAL_LOCI = frozenset({"BoLA-1", "BoLA-2", "BoLA-3", "BoLA-4", "BoLA-5", "BoLA-6"})


@dataclass(frozen=True)
class Alignment:
    allele: str
    locus: str
    contig: str
    tstart: int
    tend: int
    nmatch: int
    identity: float
    coverage: float
    mapq: int
    is_primary: bool


@dataclass(frozen=True)
class CopyCluster:
    contig: str
    tstart: int
    tend: int
    hits: tuple[Alignment, ...]
    best: Alignment
    n_distinct_alleles: int


def assign_probe(identity: float, coverage: float, mapq: int) -> str:
    """Classify one probe alignment. Low identity is never a locus call."""
    if identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE and mapq >= MIN_MAPQ:
        return HIGH_CONFIDENCE
    return UNRESOLVED


def parse_paf(text: str) -> list[Alignment]:
    """Parse a minimap2 PAF; identity comes from ``dv:f:``, not the CIGAR span.

    Splice mode inflates ``alnlen`` with introns, so ``nmatch/alnlen`` is not
    a usable identity. Coverage is aligned query bases over query length.
    """
    alignments: list[Alignment] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        query = fields[0]
        allele, locus = (query.split("|") + [""])[:2]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        contig = fields[5]
        tstart = int(fields[7])
        tend = int(fields[8])
        nmatch = int(fields[9])
        mapq = int(fields[11])
        coverage = (qend - qstart) / qlen if qlen else 0.0
        divergence = None
        is_primary = True
        for tag in fields[12:]:
            if tag.startswith("dv:f:"):
                divergence = float(tag[5:])
            elif tag.startswith("tp:A:"):
                is_primary = tag[5:] == "P"
        identity = (1.0 - divergence) if divergence is not None else 0.0
        alignments.append(
            Alignment(
                allele=allele,
                locus=locus,
                contig=contig,
                tstart=tstart,
                tend=tend,
                nmatch=nmatch,
                identity=identity,
                coverage=coverage,
                mapq=mapq,
                is_primary=is_primary,
            )
        )
    return alignments


def cluster_alignments(
    alignments: Sequence[Alignment],
    gap: int = CLUSTER_GAP,
    max_span: int = MAX_SPAN,
) -> list[CopyCluster]:
    """Group high-confidence hits that sit close together on one contig.

    Eligible hits are grouped by contig. Nearby hits (within ``gap``), including
    paralogous IPD alleles of a different locus name, become one copy whose
    locus is the single best allele. That stops a BoLA-1 panel sequence that
    also matches a BoLA-2 gene from inventing a second locus at the same
    coordinates. Two genes farther apart than ``gap`` stay separate copies.
    Over-long splice alignments (span > ``max_span``) are dropped first.
    """
    eligible = [
        hit
        for hit in alignments
        if assign_probe(hit.identity, hit.coverage, hit.mapq) == HIGH_CONFIDENCE
        and (hit.tend - hit.tstart) <= max_span
    ]
    by_contig: dict[str, list[Alignment]] = defaultdict(list)
    for hit in eligible:
        by_contig[hit.contig].append(hit)

    copies: list[CopyCluster] = []
    for contig, hits in by_contig.items():
        ordered = sorted(hits, key=lambda hit: (hit.tstart, hit.tend))
        current: list[Alignment] = []
        cluster_end: int | None = None
        for hit in ordered:
            if current and hit.tstart > (cluster_end or 0) + gap:
                copies.append(_finish_cluster(contig, current))
                current = []
                cluster_end = None
            current.append(hit)
            cluster_end = max(cluster_end or hit.tend, hit.tend)
        if current:
            copies.append(_finish_cluster(contig, current))
    return copies


def _finish_cluster(contig: str, hits: Sequence[Alignment]) -> CopyCluster:
    best = max(hits, key=lambda hit: (hit.nmatch, hit.identity, hit.is_primary, hit.mapq))
    return CopyCluster(
        contig=contig,
        tstart=min(hit.tstart for hit in hits),
        tend=max(hit.tend for hit in hits),
        hits=tuple(hits),
        best=best,
        n_distinct_alleles=len({hit.allele for hit in hits}),
    )


def classify_locus_structure(
    n_contigs: int,
    n_clusters: int,
    n_distinct_alleles_in_largest_cluster: int = 1,
) -> str:
    """Name the assembly geometry of one locus. Never a copy-number call.

    Several IPD alleles routinely hit one gene copy; that panel competition is
    not haplotype collapse, so ``n_distinct_alleles_in_largest_cluster`` does
    not change the geometry label.
    """
    if n_clusters <= 0:
        return NOT_RECOVERED
    if n_contigs > 1:
        return HAPLOTYPE_SPLIT
    if n_clusters > 1:
        return DUPLICATED_ON_CONTIG
    return RECOVERED


def locus_structure(copies: Sequence[CopyCluster], locus: str) -> str:
    members = [copy for copy in copies if copy.best.locus == locus]
    if not members:
        return NOT_RECOVERED
    largest = max(members, key=lambda copy: copy.n_distinct_alleles)
    return classify_locus_structure(
        n_contigs=len({copy.contig for copy in members}),
        n_clusters=len(members),
        n_distinct_alleles_in_largest_cluster=largest.n_distinct_alleles,
    )


def is_classical(locus: str) -> bool:
    return locus in CLASSICAL_LOCI
