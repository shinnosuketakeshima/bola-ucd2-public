"""Classify DRB3 and rs29026690 placement on wide-interval assemblies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from scripts.codex_reanalysis.dq_dr_recovery import normalize_drb3_allele

RS29026690_FLANK = "NC_037350.1:27601651-27605651"
DRB3_REGION_FLANK = "NC_037350.1:25720716-25737818"
OFFSET_MIN_BP = 2000
OFFSET_MAX_BP = 3000


@dataclass(frozen=True)
class MapHit:
    query: str
    contig: str
    pos: int
    mapq: int


def extract_drb3_panel_sequences(panel_text: str, hap1: str, hap2: str) -> str:
    """Return FASTA for PCR-SBT DRB3 alleles from the IPD-MHC panel."""
    targets = {f"BoLA-DRB3{allele}" for allele in {hap1, hap2}}
    keep = False
    lines: list[str] = []
    for line in panel_text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name in targets
        if keep:
            lines.append(line)
    return "\n".join(lines) + ("\n" if lines else "")


def best_hit_per_query(hits: Sequence[MapHit]) -> dict[str, MapHit]:
    best: dict[str, MapHit] = {}
    for hit in hits:
        previous = best.get(hit.query)
        if previous is None or (hit.mapq, hit.pos) > (previous.mapq, previous.pos):
            best[hit.query] = hit
    return best


def classify_drb3_resolution(
    hap1: str,
    hap2: str,
    allele_hits: Sequence[MapHit],
) -> tuple[str, str]:
    """Return ``(resolution_label, note)`` for Table 7 style reporting."""
    if hap1 == hap2:
        if not allele_hits:
            return "homozygous_control", "no DRB3 hits"
        best = max(allele_hits, key=lambda item: (item.mapq, item.pos))
        return "homozygous_control", f"single best hit mapq{best.mapq}"

    by_query = best_hit_per_query(allele_hits)
    matched = {
        normalize_drb3_allele(hit.query.replace("BoLA-DRB3", ""))
        for hit in by_query.values()
        if hit.query.startswith("BoLA-DRB3")
    }
    if len(matched) < 2:
        return "incomplete", "fewer than two allele hits"

    hits = list(by_query.values())
    contigs = {hit.contig for hit in hits}
    if len(contigs) > 1:
        return "split_contigs", "alleles on different contigs"

    positions = sorted(hit.pos for hit in hits)
    offset = positions[-1] - positions[0]
    mapqs = [hit.mapq for hit in hits]
    if offset == 0:
        low_mapq = [hit for hit in hits if hit.mapq < 10]
        if low_mapq:
            unreliable = low_mapq[0]
            partner = hits[0] if hits[0] is not unreliable else hits[1]
            return (
                "collapsed",
                f"{unreliable.query.split('*', 1)[-1]} mapq{unreliable.mapq} (unreliable) "
                f"vs {partner.query.split('*', 1)[-1]} mapq{partner.mapq}",
            )
        return "collapsed", f"identical position, both mapq{min(mapqs)}"

    if OFFSET_MIN_BP <= offset <= OFFSET_MAX_BP:
        degraded = [hit for hit in hits if hit.mapq < 58]
        if degraded:
            allele = degraded[0].query.split("*", 1)[-1]
            partner = hits[0] if hits[0] is not degraded[0] else hits[1]
            return (
                "collapsed",
                f"offset artifact, {allele} mapq{degraded[0].mapq} vs "
                f"{partner.query.split('*', 1)[-1]} mapq{partner.mapq}, offset~{offset}bp",
            )
        return "collapsed", f"offset artifact, both mapq{min(mapqs)}, offset~{offset}bp"

    return "collapsed", f"same contig offset {offset}bp"


def summarize_linkage(
    hap1: str,
    hap2: str,
    allele_hits: Sequence[MapHit],
    rs_hits: Sequence[MapHit],
) -> dict[str, str]:
    """Build one per-sample linkage row."""
    resolution, note = classify_drb3_resolution(hap1, hap2, allele_hits)
    drb3_contig = "-"
    if allele_hits:
        best = max(allele_hits, key=lambda item: (item.mapq, item.pos))
        drb3_contig = best.contig
    rs_contig = "-"
    if rs_hits:
        best_rs = max(rs_hits, key=lambda item: (item.mapq, item.pos))
        rs_contig = best_rs.contig
    same = (
        "yes"
        if drb3_contig != "-" and rs_contig != "-" and drb3_contig == rs_contig
        else "no"
    )
    return {
        "drb3_diplotype": f"{hap1}/{hap2}",
        "drb3_contig": drb3_contig,
        "rs29026690_contig": rs_contig,
        "same_contig": same,
        "drb3_resolution": resolution,
        "note": note,
    }
