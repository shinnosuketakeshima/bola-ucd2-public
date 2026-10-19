"""Compare the DQ deletion junction across the contigs that carry it.

Stage 06b (`attribute_dq_contigs`) reports a deletion junction per contig as
the pair of reference alignment anchors that bracket it. Contigs whose
anchors end and resume at the same reference coordinates share *flanking
alignment boundaries*; that alone does not show that they share one junction
sequence. This stage tests that for the junction selected by 06b:

1. Every tier-1, read-supported deletion with the selected event's reference
   boundaries is taken from `dq_contig_events.tsv`. The contig interval
   between its anchors, plus `--flank` bp on each side, is cut from the
   assembly and put in reference orientation.
2. Each interval is compared with the reference around the deletion by
   `blastn`. The longest internal hit (identity >= `--min-internal-identity`,
   length >= `--min-internal-length`), extended by hits adjacent to it on
   both sequences (within 200 bp; length >= a third of the minimum), gives
   the reference sequence the contig carries between the anchors. The query overlap between the left anchor
   and the first internal hit, and between the last internal hit and the
   right anchor, is the interval within which the breakpoint cannot be
   placed (microhomology or a shared repeat); the number of distinct window
   loci matching that overlap shows whether it is a repeat.
3. The between-anchor sequences are aligned globally (Edlib, NW) pair by
   pair, and the differences between the deletion-source contig and each
   other contig are binned along the source interval.

Outputs are coordinates, identities and counts; no sequence leaves
`06_local_assembly`. Contigs are collapsed Flye consensus from R9 reads, so
single-base differences are not interpreted as alleles.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.attribute_dq_contigs import CONTIG, AttributionError, fasta_sequence, parse_source
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def parse_hits(text: str, offset: int) -> list[dict]:
    """`-outfmt "6 qstart qend sstart send pident length"` lines with chromosome coordinates."""
    hits = []
    for line in text.splitlines():
        if not line.strip():
            continue
        q1, q2, s1, s2, pident, length = line.split("\t")[:6]
        r1, r2 = int(s1) + offset, int(s2) + offset
        hits.append({"q_start": int(q1), "q_end": int(q2), "r_start": min(r1, r2), "r_end": max(r1, r2),
                     "strand": "+" if r2 >= r1 else "-", "identity": float(pident), "length": int(length)})
    return hits


def junction_structure(hits: Sequence[dict], query_length: int, flank: int, reference_left: int,
                       reference_right: int, min_identity: float, min_length: int, adjacency: int = 200) -> dict:
    """Anchors, internal reference matches and boundary overlaps of one oriented junction interval.

    The query is `flank` bp of left anchor, the between-anchor sequence, and
    `flank` bp of right anchor, in reference orientation.
    """
    plus = [h for h in hits if h["strand"] == "+"]
    left = [h for h in plus if h["q_start"] <= 1 + 50 and h["r_end"] >= reference_left - 50
            and h["r_end"] <= reference_left + 2000]
    right = [h for h in plus if h["q_end"] >= query_length - 50 and h["r_start"] <= reference_right + 50
             and h["r_start"] >= reference_right - 2000]
    if not left or not right:
        raise AttributionError("junction interval lacks a left or right reference anchor")
    left_anchor = max(left, key=lambda h: h["length"])
    right_anchor = max(right, key=lambda h: h["length"])
    candidates = [h for h in plus if h["identity"] >= min_identity and h["length"] >= min_length // 3
                  and h["r_start"] > reference_left and h["r_end"] < reference_right
                  and h["q_end"] > flank and h["q_start"] < query_length - flank]
    internal = []
    seeds = [h for h in candidates if h["length"] >= min_length]
    if seeds:
        # the longest internal match, extended by matches adjacent to it on both sequences
        internal = [max(seeds, key=lambda h: h["length"])]
        extended = True
        while extended:
            extended = False
            for h in candidates:
                if h in internal:
                    continue
                first, last = internal[0], internal[-1]
                if abs(h["q_start"] - last["q_end"]) <= adjacency and abs(h["r_start"] - last["r_end"]) <= adjacency:
                    internal.append(h)
                    extended = True
                elif abs(first["q_start"] - h["q_end"]) <= adjacency and abs(first["r_start"] - h["r_end"]) <= adjacency:
                    internal.insert(0, h)
                    extended = True
    row = {
        "left_anchor_reference_end": left_anchor["r_end"],
        "right_anchor_reference_start": right_anchor["r_start"],
        "internal_hits": len(internal),
        "internal_reference_start": min((h["r_start"] for h in internal), default=""),
        "internal_reference_end": max((h["r_end"] for h in internal), default=""),
        "internal_min_identity": min((h["identity"] for h in internal), default=""),
        "left_boundary_overlap_bp": "",
        "right_boundary_overlap_bp": "",
    }
    if internal:
        row["left_boundary_overlap_bp"] = max(left_anchor["q_end"] - internal[0]["q_start"] + 1, 0)
        row["right_boundary_overlap_bp"] = max(internal[-1]["q_end"] - right_anchor["q_start"] + 1, 0)
        row["left_overlap_query"] = (internal[0]["q_start"], left_anchor["q_end"])
        row["right_overlap_query"] = (right_anchor["q_start"], internal[-1]["q_end"])
    return row


def overlap_loci(hits: Sequence[dict], interval: tuple[int, int], min_length: int = 300) -> int:
    """Distinct reference loci matched by >= min_length bp of a query interval (repeat indicator)."""
    low, high = interval
    if high < low:
        return 0
    loci = set()
    for h in hits:
        if min(h["q_end"], high) - max(h["q_start"], low) + 1 >= min(min_length, high - low + 1):
            loci.add(h["r_start"] // 1000)
    return len(loci)


def global_identity(a: str, b: str) -> tuple[int, int, float]:
    import edlib

    result = edlib.align(a, b, mode="NW", task="path")
    columns = len(edlib.getNiceAlignment(result, a, b)["matched_aligned"])
    return result["editDistance"], columns, 1 - result["editDistance"] / columns


def difference_bins(source: str, other: str, width: int) -> list[int]:
    """Mismatch and gap columns per `width` bp of the source sequence."""
    import edlib

    result = edlib.align(source, other, mode="NW", task="path")
    nice = edlib.getNiceAlignment(result, source, other)
    bins = [0] * (len(source) // width + 1)
    position = 0
    for base, match in zip(nice["query_aligned"], nice["matched_aligned"]):
        if base != "-":
            position += 1
        if match != "|":
            bins[min(max(position - 1, 0) // width, len(bins) - 1)] += 1
    return bins


def _run(command: Sequence[str]) -> str:
    return subprocess.run(list(command), check=True, capture_output=True, text=True).stdout


def write_tsv(path: Path, rows: list[dict], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--attribution", type=Path, required=True, help="06b_dq_attribution")
    parser.add_argument("--assembly-dir", type=Path, required=True, help="06_local_assembly/dq_dr")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=3000)
    parser.add_argument("--window-padding", type=int, default=20000)
    parser.add_argument("--min-internal-identity", type=float, default=90.0)
    parser.add_argument("--min-internal-length", type=int, default=1000)
    parser.add_argument("--bin", type=int, default=500)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--blastn", default="blastn")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    status = json.loads((args.attribution / "dq_attribution_status.json").read_text(encoding="utf-8"))
    selected = status["selected_deletion"]
    source = (selected["sample_id"], selected["contig"])
    boundaries = (int(selected["reference_left"]), int(selected["reference_right"]))
    with (args.attribution / "dq_contig_events.tsv").open(encoding="utf-8", newline="") as handle:
        events = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if row["tier"] == "1" and row["event"] == "deletion" and row["read_supported"] == "True"
            and (int(row["reference_left"]), int(row["reference_right"])) == boundaries
        ]
    if (source[0], source[1]) not in {(r["sample_id"], r["contig"]) for r in events}:
        raise AttributionError(f"{':'.join(source)} is not among the read-supported junctions")
    events.sort(key=lambda r: ((r["sample_id"], r["contig"]) != source, r["sample_id"]))

    window_low = boundaries[0] - args.window_padding
    window_high = boundaries[1] + args.window_padding
    rows, sequences = [], {}
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        window = workdir / "window.fa"
        window.write_text(_run([args.samtools, "faidx", str(args.reference), f"{CONTIG}:{window_low}-{window_high}"]),
                          encoding="utf-8")
        for event in events:
            sample, contig = event["sample_id"], event["contig"]
            low = int(event["contig_left"]) - args.flank + 1
            high = int(event["contig_right"]) + args.flank - 1
            sequence = fasta_sequence(_run([args.samtools, "faidx", str(args.assembly_dir / sample / "assembly.fasta"),
                                            f"{contig}:{low}-{high}"]))
            query = workdir / f"{sample}.fa"
            query.write_text(f">q\n{sequence}\n", encoding="utf-8")
            hits = parse_hits(_run([args.blastn, "-query", str(query), "-subject", str(window), "-evalue", "1e-20",
                                    "-outfmt", "6 qstart qend sstart send pident length"]), window_low - 1)
            longest = max(hits, key=lambda h: h["length"])
            orientation = longest["strand"]
            if orientation == "-":
                sequence = reverse_complement(sequence)
                query.write_text(f">q\n{sequence}\n", encoding="utf-8")
                hits = parse_hits(_run([args.blastn, "-query", str(query), "-subject", str(window), "-evalue",
                                        "1e-20", "-outfmt", "6 qstart qend sstart send pident length"]),
                                  window_low - 1)
            structure = junction_structure(hits, len(sequence), args.flank, *boundaries,
                                           args.min_internal_identity, args.min_internal_length)
            structure["left_overlap_loci"] = overlap_loci(hits, structure.pop("left_overlap_query", (1, 0)))
            structure["right_overlap_loci"] = overlap_loci(hits, structure.pop("right_overlap_query", (1, 0)))
            sequences[sample] = sequence[args.flank:len(sequence) - args.flank]
            rows.append({"sample_id": sample, "contig": contig, "contig_orientation": orientation,
                         "contig_gap_bp": event["contig_gap_bp"], "spanning_reads": event["spanning_reads"],
                         **structure})

    pairs = []
    for a, b in itertools.combinations(sequences, 2):
        edits, columns, identity = global_identity(sequences[a], sequences[b])
        pairs.append({"sample_a": a, "sample_b": b, "length_a": len(sequences[a]), "length_b": len(sequences[b]),
                      "edits": edits, "aligned_columns": columns, "identity": f"{identity:.4f}"})
    bins = []
    for other in list(sequences)[1:]:
        for k, count in enumerate(difference_bins(sequences[source[0]], sequences[other], args.bin)):
            bins.append({"source": source[0], "other": other, "source_bin_start": k * args.bin + 1,
                         "differences": count})

    args.output.mkdir(parents=True, exist_ok=True)
    structure_fields = ["sample_id", "contig", "contig_orientation", "contig_gap_bp", "spanning_reads",
                        "left_anchor_reference_end", "right_anchor_reference_start", "internal_hits",
                        "internal_reference_start", "internal_reference_end", "internal_min_identity",
                        "left_boundary_overlap_bp", "left_overlap_loci", "right_boundary_overlap_bp",
                        "right_overlap_loci"]
    write_tsv(args.output / "dq_junction_structure.tsv", rows, structure_fields)
    write_tsv(args.output / "dq_junction_pairwise.tsv", pairs,
              ["sample_a", "sample_b", "length_a", "length_b", "edits", "aligned_columns", "identity"])
    write_tsv(args.output / "dq_junction_difference_bins.tsv", bins,
              ["source", "other", "source_bin_start", "differences"])
    parameters = {
        "selected_deletion": f"{source[0]}:{source[1]}",
        "reference_boundaries": list(boundaries),
        "flank_bp": args.flank,
        "window": f"{CONTIG}:{window_low}-{window_high}",
        "internal_hits": f"blastn identity >= {args.min_internal_identity}%, length >= {args.min_internal_length} bp",
        "pairwise": "Edlib NW over the between-anchor contig sequence",
        "bin_bp": args.bin,
    }
    (args.output / "dq_junction_status.json").write_text(
        json.dumps({"status": "passed", "junctions": len(rows), "parameters": parameters}, indent=2,
                   sort_keys=True) + "\n", encoding="utf-8")
    write_provenance(args.output, {
        "inputs": {
            "events_sha256": sha256_file(args.attribution / "dq_contig_events.tsv"),
            "attribution_status_sha256": sha256_file(args.attribution / "dq_attribution_status.json"),
            "reference_sha256": sha256_file(args.reference),
        },
        "parameters": parameters,
    })
    print(f"DQ junction comparison written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
