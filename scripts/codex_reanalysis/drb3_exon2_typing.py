"""Exon-2 DRB3 typing from target BAMs, scored on a window every panel allele covers.

``drb3_wgs_typing.py`` (stage 03c) maps whole reads to the IPD DRB3 panel and
keeps each read's best allele by identity, then by matching bases. The panel
mixes lengths -- most entries are exon 2 only (234-270 bp), but *016:01 and
*027:03 are full CDS (780-801 bp) -- so the matching-base tie-break favours the
long entries. That is the "central-allele bias" 03c reports, and it is a
property of the scoring, not of the reads.

This stage removes the bias by construction:

1. The DRB3 exon-2 interval comes from the verified GFF3 (second exon of the
   BOLA-DRB3 transcript), never from a hard-coded coordinate.
2. Every panel allele is located inside the reference exon 2 and trimmed to
   the *window covered by every allele* (max of the starts, min of the ends).
   The window is a property of the panel alone; it does not look at PCR-SBT.
3. Reads that span exon 2 are cut to the same window.
4. Per sample, the diplotype is the allele pair minimising the summed
   per-read edit distance, each read charged to the closer allele of the pair.
   Alleles whose window sequence is identical form one group and cannot be
   told apart by exon 2; resolution is reported per group.

PCR-SBT is used only afterwards, as the comparison.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

DRB3_GENE = "BOLA-DRB3"
EXON_NUMBER = 2
EXTRACT_FLAGS = "0x900"
READ_FLANK = 15
# A trimmed read segment is kept only if its window anchors with at most this
# fraction of edits; it rejects reads that merely overlap the interval.
MAX_ANCHOR_ERROR = 0.35
# Robustness: shrink the derived window symmetrically by these margins.
ROBUSTNESS_MARGINS = (0, 8, 18, 38)
ALLELE_RE = re.compile(r"BoLA-DRB3\*(\d+):(\d+)")


class Drb3Exon2Error(ValueError):
    """Raised when inputs cannot be analysed safely."""


# --------------------------------------------------------------------------
# Pure helpers (no edlib needed)
# --------------------------------------------------------------------------
def normalize_allele(name: str) -> str:
    """``BoLA-DRB3*015:01:02`` and ``*015:01`` both become ``*015:01``."""
    match = ALLELE_RE.search(name)
    if match:
        return f"*{match.group(1)}:{match.group(2)}"
    match = re.match(r"\*(\d+):(\d+)", name)
    if match:
        return f"*{match.group(1)}:{match.group(2)}"
    raise Drb3Exon2Error(f"not a DRB3 allele label: {name!r}")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line)
    return {k: "".join(v).upper() for k, v in records.items()}


def drb3_panel(path: Path) -> dict[str, str]:
    panel = {k: v for k, v in read_fasta(path).items() if k.startswith("BoLA-DRB3*")}
    if not panel:
        raise Drb3Exon2Error(f"no BoLA-DRB3 sequences in {path}")
    return panel


def exon_from_gff(gff: Path, gene_name: str = DRB3_GENE, exon_number: int = EXON_NUMBER):
    """Return ``(contig, start_1based, end_1based, strand)`` of one exon of the
    first transcript of ``gene_name``, numbered in transcript order."""
    gene_id = None
    contig = strand = None
    mrna_id = None
    exons: list[tuple[int, int]] = []
    with gff.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = dict(
                item.split("=", 1) for item in cols[8].split(";") if "=" in item
            )
            if cols[2] == "gene" and attrs.get("Name") == gene_name:
                if gene_id is not None:
                    raise Drb3Exon2Error(f"more than one gene named {gene_name}")
                gene_id, contig, strand = attrs["ID"], cols[0], cols[6]
            elif gene_id and cols[2] == "mRNA" and attrs.get("Parent") == gene_id:
                if mrna_id is None:
                    mrna_id = attrs["ID"]
            elif mrna_id and cols[2] == "exon" and attrs.get("Parent") == mrna_id:
                exons.append((int(cols[3]), int(cols[4])))
    if gene_id is None:
        raise Drb3Exon2Error(f"gene {gene_name} not found in {gff}")
    if len(exons) < exon_number:
        raise Drb3Exon2Error(f"{gene_name} has {len(exons)} exons, need {exon_number}")
    exons.sort(reverse=(strand == "-"))
    start, end = exons[exon_number - 1]
    return contig, start, end, strand


def query_segment(pos1: int, cigar: str, seq: str, ref_start0: int, ref_end0: int):
    """Query bases aligned to ``[ref_start0, ref_end0)`` (0-based, half-open).

    ``pos1`` is the 1-based SAM POS. Returns None unless the alignment spans the
    whole interval."""
    ref = pos1 - 1
    qpos = 0
    q_first = q_last = None
    for length_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length_text)
        if op in "M=X":
            block_start, block_end = ref, ref + length
            lo, hi = max(block_start, ref_start0), min(block_end, ref_end0)
            if lo < hi:
                if q_first is None:
                    q_first = qpos + (lo - ref)
                q_last = qpos + (hi - ref)
            ref += length
            qpos += length
        elif op in "IS":
            qpos += length
        elif op in "DN":
            ref += length
        # H and P consume neither
    if q_first is None or pos1 - 1 > ref_start0 or ref < ref_end0:
        return None
    return seq[q_first:q_last]


def group_identical(sequences: dict[str, str]):
    """Group allele names by identical sequence. Returns ``(seqs, names)``
    where ``names[i]`` are the normalised allele names sharing ``seqs[i]``."""
    groups: dict[str, set[str]] = {}
    for name, seq in sequences.items():
        groups.setdefault(seq, set()).add(normalize_allele(name))
    seqs = list(groups)
    return seqs, [sorted(groups[s]) for s in seqs]


def diplotype_costs(distances: Sequence[Sequence[int]], n_groups: int):
    """Cost of every unordered group pair: sum over reads of the smaller
    distance to either member. Returns ``{(a, b): cost}`` with a <= b."""
    costs: dict[tuple[int, int], int] = {}
    for a, b in combinations_with_replacement(range(n_groups), 2):
        costs[(a, b)] = sum(min(row[a], row[b]) for row in distances)
    return costs


def rank_diplotypes(costs: dict[tuple[int, int], int]):
    """``(best_cost, tied_best_pairs, second_cost)``; second is the next
    distinct cost, or None when every pair ties."""
    ordered = sorted(set(costs.values()))
    best = ordered[0]
    tied = sorted(pair for pair, cost in costs.items() if cost == best)
    return best, tied, (ordered[1] if len(ordered) > 1 else None)


def group_indices(allele: str, names: Sequence[Sequence[str]]) -> list[int]:
    return [i for i, group in enumerate(names) if allele in group]


def expected_cost(costs, hap1_groups: Iterable[int], hap2_groups: Iterable[int]):
    options = [
        costs[tuple(sorted((a, b)))] for a in hap1_groups for b in hap2_groups
    ]
    return min(options) if options else None


def pair_label(pair: tuple[int, int], names: Sequence[Sequence[str]]) -> str:
    return " | ".join("/".join(names[i]) for i in pair)


# --------------------------------------------------------------------------
# edlib-dependent steps
# --------------------------------------------------------------------------
def _edlib():
    try:
        import edlib
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise Drb3Exon2Error(
            "edlib is required (pip install edlib / conda-forge edlib)"
        ) from exc
    return edlib


def infix_location(query: str, text: str):
    """``(start, end_exclusive, edit_distance)`` of the best infix hit."""
    edlib = _edlib()
    result = edlib.align(query, text, mode="HW", task="locations")
    start, end = result["locations"][0]
    return start, end + 1, result["editDistance"]


def trim_to(window: str, text: str):
    start, end, edits = infix_location(window, text)
    return text[start:end], edits


def cut_to_exon(panel: dict[str, str], reference_exon: str) -> dict[str, str]:
    """Cut full-CDS panel entries to the stretch matching the reference exon."""
    cut = {}
    for name, seq in panel.items():
        if len(seq) > len(reference_exon) + 30:
            start, end, _ = infix_location(reference_exon, seq)
            seq = seq[start:end]
        cut[name] = seq
    return cut


def panel_window(exon_panel: dict[str, str], reference_exon: str):
    """Intersection, in reference-exon coordinates, of every panel allele.

    Returns ``(start, end, rows)``; the window depends on the panel alone."""
    rows = []
    starts, ends = [], []
    for name, seq in sorted(exon_panel.items()):
        start, end, edits = infix_location(seq, reference_exon)
        rows.append(
            {"allele": name, "length": len(seq), "ref_start": start, "ref_end": end, "edits": edits}
        )
        starts.append(start)
        ends.append(end)
    return max(starts), min(ends), rows


def distance_matrix(segments: Sequence[str], seqs: Sequence[str]) -> list[list[int]]:
    edlib = _edlib()
    return [
        [edlib.align(q, t, mode="NW", task="distance")["editDistance"] for q in seqs]
        for t in segments
    ]


# --------------------------------------------------------------------------
# BAM access through samtools (no pysam)
# --------------------------------------------------------------------------
def spanning_segments(
    bam: Path, contig: str, start0: int, end0: int, *, samtools: str
) -> list[str]:
    """Read bases over ``[start0 - flank, end0 + flank)`` for primary
    alignments that span the exon with flanks, in reference orientation."""
    lo, hi = start0 - READ_FLANK, end0 + READ_FLANK
    region = f"{contig}:{lo + 1}-{hi}"
    proc = subprocess.run(
        [samtools, "view", "-F", EXTRACT_FLAGS, str(bam), region],
        check=True,
        capture_output=True,
        text=True,
    )
    segments = []
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 11 or cols[9] == "*":
            continue
        segment = query_segment(int(cols[3]), cols[5], cols[9].upper(), start0, end0)
        if segment is not None:
            flanked = query_segment(int(cols[3]), cols[5], cols[9].upper(), lo, hi)
            segments.append(flanked if flanked is not None else segment)
    return segments


def window_segments(
    raw_segments: Sequence[str], window: str
) -> tuple[list[str], int]:
    kept, dropped = [], 0
    for text in raw_segments:
        trimmed, edits = trim_to(window, text)
        if edits <= MAX_ANCHOR_ERROR * len(window):
            kept.append(trimmed)
        else:
            dropped += 1
    return kept, dropped


# --------------------------------------------------------------------------
# One sample
# --------------------------------------------------------------------------
def type_sample(
    sample_id: str,
    raw_segments: Sequence[str],
    window: str,
    seqs: Sequence[str],
    names: Sequence[Sequence[str]],
    hap1: str,
    hap2: str,
) -> dict[str, object]:
    reads, dropped = window_segments(raw_segments, window)
    if not reads:
        raise Drb3Exon2Error(f"{sample_id}: no reads span the DRB3 exon-2 window")
    distances = distance_matrix(reads, seqs)
    costs = diplotype_costs(distances, len(seqs))
    best, tied, second = rank_diplotypes(costs)
    g1 = group_indices(normalize_allele(hap1), names)
    g2 = group_indices(normalize_allele(hap2), names)
    if not g1 or not g2:
        raise Drb3Exon2Error(f"{sample_id}: PCR-SBT allele {hap1}/{hap2} absent from panel")
    exp_cost = expected_cost(costs, g1, g2)
    called = tied[0]
    read_error = sorted(min(row) for row in distances)
    return {
        "sample_id": sample_id,
        "pcr_sbt": "/".join(sorted((hap1, hap2))),
        "called_best_diplotype": pair_label(called, names),
        "n_spanning_reads": len(reads),
        "n_reads_dropped_unanchored": dropped,
        "median_read_edits_to_best_allele": read_error[len(read_error) // 2],
        "best_cost": best,
        "expected_cost": exp_cost,
        "gap_expected_minus_best": exp_cost - best,
        "concordant": exp_cost == best,
        "n_tied_best_pairs": len(tied),
        "margin_to_second_best": None if second is None else second - best,
        "pcr_sbt_group_sizes": f"{len(names[g1[0]])},{len(names[g2[0]])}",
    }


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------
def render_report(rows, robustness, window_info) -> str:
    n = len(rows)
    ok = sum(1 for r in rows if r["concordant"])
    lines = [
        "# DRB3 exon-2 typing from target BAMs (ARS-UCD2.0)",
        "",
        "Replaces the per-read, best-allele scoring of stage 03c, whose ranking by matching "
        "bases favours the full-length panel entries (*016:01, *027:03). Here every panel "
        "allele and every read is cut to one shared exon-2 window and the diplotype is fitted "
        "jointly. PCR-SBT is used only for comparison.",
        "",
        "## Window",
        "",
        f"- Exon 2 (GFF3): `{window_info['contig']}:{window_info['exon_start']}-{window_info['exon_end']}` "
        f"({window_info['strand']})",
        f"- Window covered by all {window_info['n_panel_alleles']} panel alleles: reference-exon "
        f"offsets {window_info['window_start']}-{window_info['window_end']} "
        f"({window_info['window_length']} nt)",
        f"- Distinct window sequences (groups): {window_info['n_groups']}",
        "",
        "## Concordance with PCR-SBT",
        "",
        f"- Best-cost diplotype equals the PCR-SBT diplotype (cost gap 0): **{ok}/{n}**",
        "",
        "| sample | PCR-SBT | best diplotype (sequence groups) | reads | cost gap | margin to 2nd best | concordant |",
        "| --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for r in rows:
        lines.append(
            f"| {r['sample_id']} | {r['pcr_sbt']} | {r['called_best_diplotype']} | "
            f"{r['n_spanning_reads']} | {r['gap_expected_minus_best']} | "
            f"{r['margin_to_second_best']} | {r['concordant']} |"
        )
    lines += ["", "## Window robustness", "", "| margin trimmed each side (nt) | window (nt) | concordant |", "| ---: | ---: | ---: |"]
    for margin, length, count, total in robustness:
        lines.append(f"| {margin} | {length} | {count}/{total} |")
    lines += [
        "",
        "## Limits",
        "",
        "- Resolution is by exon-2 window sequence group; alleles identical in the window "
        "(e.g. *016:01/*100:01) cannot be separated.",
        "- No heterozygosity penalty is applied; a homozygote can be fitted with a second "
        "allele that differs by a base or two when the data are noisy.",
        "- The window rule was fixed after the 03c scoring was found to be length-biased, and "
        "the method is evaluated on the same 13 animals; there is no held-out set.",
        "",
    ]
    return "\n".join(lines)


def write_tsv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_registry(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    analytical = [r for r in rows if r.get("analysis_status") == "analytical" and r.get("sample_id")]
    if not analytical:
        raise Drb3Exon2Error(f"no analytical samples in {path}")
    return analytical


def verified_bam(mapping_dir: Path, sample_id: str) -> tuple[Path, str]:
    """Return the target BAM and its SHA-256, refusing a BAM that does not match
    the hash its mapping status recorded."""
    status_path = mapping_dir / f"{sample_id}.mapping_status.json"
    bam = mapping_dir / f"{sample_id}.targets.bam"
    if not status_path.exists() or not bam.exists():
        raise Drb3Exon2Error(f"missing mapping status or BAM for {sample_id}")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "passed":
        raise Drb3Exon2Error(f"mapping status for {sample_id} is {status.get('status')!r}")
    digest = sha256_file(bam)
    recorded = status.get("artifact_hashes", {}).get("target_bam")
    if recorded != digest:
        raise Drb3Exon2Error(f"{sample_id}: BAM hash differs from its mapping status")
    return bam, digest


def reference_exon_sequence(fasta: Path, contig: str, start: int, end: int, samtools: str) -> str:
    proc = subprocess.run(
        [samtools, "faidx", str(fasta), f"{contig}:{start}-{end}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return "".join(proc.stdout.splitlines()[1:]).upper()


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=Path("docs/codex_reanalysis/DRB3_EXON2_TYPING_REPORT.md"))
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--samples", default="")
    args = parser.parse_args(list(argv) if argv is not None else None)

    edlib = _edlib()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    contig, exon_start, exon_end, strand = exon_from_gff(args.gff)
    if strand != "+":
        raise Drb3Exon2Error("minus-strand DRB3 would need reverse-complemented alleles")
    start0, end0 = exon_start - 1, exon_end
    ref_exon = reference_exon_sequence(args.reference_fasta, contig, exon_start, exon_end, args.samtools)
    panel = cut_to_exon(drb3_panel(args.panel), ref_exon)
    w_start, w_end, panel_rows = panel_window(panel, ref_exon)
    if w_end - w_start < 100:
        raise Drb3Exon2Error(f"panel-covered window is only {w_end - w_start} nt")

    registry = load_registry(args.registry)
    if args.samples.strip():
        wanted = {s.strip() for s in args.samples.split(",") if s.strip()}
        registry = [r for r in registry if r["sample_id"] in wanted]
        missing = wanted - {r["sample_id"] for r in registry}
        if missing:
            raise Drb3Exon2Error(f"unknown sample_id(s): {sorted(missing)}")

    per_margin: dict[int, list[dict[str, object]]] = {m: [] for m in ROBUSTNESS_MARGINS}
    bam_hashes: dict[str, str] = {}
    groups_by_margin: dict[int, int] = {}
    raw_by_sample: dict[str, list[str]] = {}
    for entry in registry:
        sample_id = entry["sample_id"]
        bam, digest = verified_bam(args.mapping_dir, sample_id)
        bam_hashes[sample_id] = digest
        raw_by_sample[sample_id] = spanning_segments(bam, contig, start0, end0, samtools=args.samtools)
        print(f"{sample_id}: {len(raw_by_sample[sample_id])} reads span exon 2", flush=True)

    for margin in ROBUSTNESS_MARGINS:
        lo, hi = w_start + margin, w_end - margin
        window = ref_exon[lo:hi]
        trimmed = {name: trim_to(window, seq)[0] for name, seq in panel.items()}
        seqs, names = group_identical(trimmed)
        groups_by_margin[margin] = len(seqs)
        for entry in registry:
            per_margin[margin].append(
                type_sample(
                    entry["sample_id"], raw_by_sample[entry["sample_id"]], window, seqs, names,
                    entry["pcr_sbt_drb3_hap1"], entry["pcr_sbt_drb3_hap2"],
                )
            )

    rows = per_margin[0]
    fields = list(rows[0])
    write_tsv(out / "drb3_exon2_typing.tsv", rows, fields)
    robustness = []
    robust_rows = []
    for margin in ROBUSTNESS_MARGINS:
        count = sum(1 for r in per_margin[margin] if r["concordant"])
        robustness.append((margin, w_end - w_start - 2 * margin, count, len(registry)))
        for r in per_margin[margin]:
            robust_rows.append({"margin": margin, "window_nt": w_end - w_start - 2 * margin, **r})
    write_tsv(out / "drb3_exon2_robustness.tsv", robust_rows, ["margin", "window_nt", *fields])
    write_tsv(out / "drb3_exon2_panel_window.tsv", panel_rows, ["allele", "length", "ref_start", "ref_end", "edits"])

    window_info = {
        "contig": contig, "exon_start": exon_start, "exon_end": exon_end, "strand": strand,
        "n_panel_alleles": len(panel), "window_start": w_start, "window_end": w_end,
        "window_length": w_end - w_start, "n_groups": groups_by_margin[0],
    }
    report = render_report(rows, robustness, window_info)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report, encoding="utf-8")
    (out / "DRB3_EXON2_TYPING_REPORT.md").write_text(report, encoding="utf-8")

    status = {
        "status": "passed",
        "n_samples": len(rows),
        "n_concordant": sum(1 for r in rows if r["concordant"]),
        "window": window_info,
        "robustness": [
            {"margin": m, "window_nt": n_nt, "n_concordant": c, "n_samples": t}
            for m, n_nt, c, t in robustness
        ],
        "tsv": str(out / "drb3_exon2_typing.tsv"),
        "report": str(args.report),
        "supersedes_for_interpretation": "03c_drb3_wgs_typing (kept unchanged)",
    }
    (out / "drb3_exon2_typing_status.json").write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_provenance(
        out,
        {
            "command": " ".join(sys.argv),
            "versions": {
                "python": sys.version.split()[0],
                "samtools": args.samtools,
                "edlib": getattr(edlib, "__version__", "unknown"),
            },
            "input_hashes": {
                "registry": sha256_file(args.registry),
                **{f"target_bam:{k}": v for k, v in sorted(bam_hashes.items())},
            },
            "reference_hashes": {
                "ipd_panel": sha256_file(args.panel),
                "reference_fasta": sha256_file(args.reference_fasta),
                "gff": sha256_file(args.gff),
            },
            "config_hash": "drb3_exon2_typing",
            "parameters": {
                "exon": f"{DRB3_GENE} exon {EXON_NUMBER}",
                "exon_interval_1based": f"{contig}:{exon_start}-{exon_end}",
                "window_rule": "intersection of all panel alleles in reference-exon coordinates",
                "window_offsets": [w_start, w_end],
                "robustness_margins": list(ROBUSTNESS_MARGINS),
                "max_anchor_error": MAX_ANCHOR_ERROR,
                "read_flank": READ_FLANK,
                "extract_flags": EXTRACT_FLAGS,
            },
        },
    )
    print(f"concordant {status['n_concordant']}/{status['n_samples']}; wrote {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
