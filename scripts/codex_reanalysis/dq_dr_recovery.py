"""Recover DQA-, DQB-, and DRB3-compatible loci from DQ/DR local assemblies.

IPD-MHC panel hits use the same evidence gate as class I (identity ≥95%,
coverage ≥80%, MAPQ ≥20). ``n_dqa_clusters`` and ``n_dqb_clusters`` are
locus-recovery summaries for the DQA and DQB gene families, not validated
diploid copy numbers. DRB3 is compared to PCR-SBT independently; assembly
does not replace PCR-SBT typing.
"""

from __future__ import annotations

import re
from typing import Sequence

from scripts.codex_reanalysis.probe_assignment import (
    HIGH_CONFIDENCE,
    Alignment,
    CopyCluster,
    assign_probe,
    cluster_alignments,
)

CLASS_II_FAMILIES = frozenset({"DQA", "DQB", "DRB3"})
DRB3_ALLELE_RE = re.compile(r"BoLA-DRB3(\*[^:]+(?::[^:]+)?)")

EXACT_TWO = "exact_two_allele_recovery"
ONE_ALLELE = "one_allele_recovery"
COLLAPSED = "collapsed"
AMBIGUOUS = "ambiguous_paralogous_mapping"
LOW_CONFIDENCE = "low_confidence_allele_match"
NO_RECOVERY = "no_reliable_recovery"

UNDER_ASSEMBLED = "under_assembled"
DIPLOID_BASELINE = "diploid_baseline"
EXPANDED = "expanded_dq"


def parse_classii_query(query: str) -> tuple[str, str, str]:
    """Return ``(allele, family, locus)`` from an IPD class II header token."""
    token = query.split()[0]
    if not token.startswith("BoLA-"):
        raise ValueError(f"unexpected IPD class II query: {query!r}")
    body = token[5:]
    if "*" not in body:
        raise ValueError(f"unexpected IPD class II query: {query!r}")
    family, _rest = body.split("*", 1)
    if family not in CLASS_II_FAMILIES:
        raise ValueError(f"unsupported class II family: {family}")
    return token, family, family


def parse_ipd_classii_paf(text: str) -> list[Alignment]:
    """Parse minimap2 PAF of an assembly against ``ipd_mhc_bola_panel.fasta``."""
    alignments: list[Alignment] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        allele, locus, _family = parse_classii_query(fields[0])
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


def normalize_drb3_allele(name: str) -> str:
    """Normalize ``BoLA-DRB3*016:01`` and ``*016:01`` to ``*016:01``."""
    if name.startswith("BoLA-DRB3"):
        return "*" + name.split("BoLA-DRB3*", 1)[1]
    match = DRB3_ALLELE_RE.search(name)
    if match:
        return match.group(1)
    if name.startswith("*"):
        return name
    raise ValueError(f"not a DRB3 allele label: {name!r}")


def count_family_clusters(copies: Sequence[CopyCluster]) -> dict[str, int]:
    counts = {family: 0 for family in CLASS_II_FAMILIES}
    for copy in copies:
        if copy.best.locus in counts:
            counts[copy.best.locus] += 1
    return counts


def best_alleles_by_family(copies: Sequence[CopyCluster]) -> dict[str, CopyCluster]:
    best: dict[str, CopyCluster] = {}
    for copy in copies:
        previous = best.get(copy.best.locus)
        if previous is None or (copy.best.nmatch, copy.best.identity) > (
            previous.best.nmatch,
            previous.best.identity,
        ):
            best[copy.best.locus] = copy
    return best


def classify_drb3_pcr_concordance(
    hap1: str,
    hap2: str,
    drb3_clusters: Sequence[CopyCluster],
    *,
    minimum_identity: float = 0.95,
) -> str:
    """Compare PCR-SBT DRB3 diplotype to assembly-recovered DRB3 alleles."""
    if not drb3_clusters:
        return NO_RECOVERY
    expected = {normalize_drb3_allele(hap1), normalize_drb3_allele(hap2)}
    matched = {
        normalize_drb3_allele(copy.best.allele)
        for copy in drb3_clusters
        if copy.best.identity >= minimum_identity
        and assign_probe(copy.best.identity, copy.best.coverage, copy.best.mapq) == HIGH_CONFIDENCE
    }
    if not matched:
        return NO_RECOVERY
    overlap = matched & expected
    if hap1 == hap2:
        if overlap:
            return ONE_ALLELE if len(drb3_clusters) == 1 else EXACT_TWO
        return LOW_CONFIDENCE
    if overlap == expected:
        if len(drb3_clusters) >= 2 and len({copy.contig for copy in drb3_clusters}) >= 2:
            return EXACT_TWO
        if len(drb3_clusters) == 1 and len({copy.best.allele for copy in drb3_clusters}) >= 2:
            return COLLAPSED
        return ONE_ALLELE
    if len(overlap) == 1:
        return ONE_ALLELE
    if len(matched) > len(expected):
        return AMBIGUOUS
    return LOW_CONFIDENCE


def classify_dq_interpretation(
    n_dqa: int,
    n_dqb: int,
    *,
    flye_status: str,
) -> str:
    """Name the DQ recovery pattern. Not a validated copy-number genotype."""
    if flye_status != "passed":
        return "failed"
    if min(n_dqa, n_dqb) < 2:
        return UNDER_ASSEMBLED
    if n_dqa >= 3 or n_dqb >= 3:
        return EXPANDED
    return DIPLOID_BASELINE
