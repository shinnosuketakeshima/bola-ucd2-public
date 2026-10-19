"""Look for the cohort's DQ configurations in published cattle assemblies.

Stage 06b describes two contig-anchored DQ configurations (the deletion
junction on the `--deletion-source` contig, and the DQA2-retaining
configuration with a non-reference segment on the `--insertion-source`
contig) and genotypes them, and the reference configuration, from reads.
This stage compares those configurations with haplotypes sequenced
elsewhere, so that a configuration is not called new, or rare, only because
ARS-UCD2.0 lacks it:

1. In each published assembly's BTA23, the DQ region is located by aligning
   the ARS-UCD2.0 interval `--locate` (`minimap2 -x asm20`); the hit span
   plus `--padding` bp on each side is cut out.
2. That region is aligned (`minimap2 -x asm20`) to each of three targets:
   the two source contigs and the ARS-UCD2.0 DQ window. Each alignment is
   scored over the 06b diagnostic windows of that target with the same
   continuity rule used for the cohort reads
   (`attribute_dq_contigs.supporting_reads`: end points enclosing the window
   with `--flank` bp each side, no long gap, aligned bases over the flanks
   and the interval), plus an identity floor: matching bases over alignment
   columns in the window >= `--min-window-identity`. Each target is aligned
   separately, not competitively, so the identity floor is what keeps a
   continuous paralogous block from counting. Primary and supplementary
   alignments both count, because one long region aligns to a contig in
   several pieces.
3. For every alignment the target span, strand and identity (matching bases
   over alignment block length) are written, so the synteny of the region
   with each contig can be read directly.
4. Each assembly is classified by a fixed rule (`classify`): configuration
   A when the deletion junction is present; configuration B when the left
   edge, the whole non-reference segment and the right edge are all
   present; the reference configuration when the reference deletion point,
   reference DQA2 and the reference BOLA-DQB-DQA5 interval are all present;
   otherwise not classified. Partial matches are not evidence for any
   configuration, and failing a reference diagnostic means only that the
   assembly is not colinear with ARS-UCD2.0 there under this rule.

A collapsed (pseudo-haploid) assembly carries one configuration per
position; for haplotype-resolved assemblies each haplotype is listed
separately. Absence of a configuration from an assembly is a statement
about that assembly only. Only coordinates, identities and accessions are
written; no sequence leaves the scratch directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.attribute_dq_contigs import (
    CONTIG,
    _run,
    load_features,
    parse_sam_alignments,
    parse_source,
    supporting_reads,
)
from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file


def read_fasta(path: Path, wanted: str) -> str:
    """The sequence of one record of a (gzipped) FASTA."""
    opener = gzip.open if path.suffix == ".gz" else open
    chunks: list[str] = []
    found = False
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == wanted
            elif found:
                chunks.append(line.strip())
    if not found:
        raise ValueError(f"{path}: no record {wanted}")
    return "".join(chunks).upper()


def locate_span(paf: str, min_block: int) -> tuple[int, int, str] | None:
    """Target span (1-based) and majority strand of the locate query's alignment blocks >= min_block."""
    blocks = []
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) >= 12 and int(f[10]) >= min_block:
            blocks.append((int(f[7]) + 1, int(f[8]), f[4], int(f[10])))
    if not blocks:
        return None
    # keep blocks near the longest one so a distant paralogous hit does not stretch the span
    anchor = max(blocks, key=lambda b: b[3])
    near = [b for b in blocks if abs(b[0] - anchor[0]) < 1_000_000]
    strand = "+" if sum(b[3] for b in near if b[2] == "+") >= sum(b[3] for b in near if b[2] == "-") else "-"
    return min(b[0] for b in near), max(b[1] for b in near), strand


def paf_blocks(paf: str, min_block: int) -> list[dict]:
    """Alignment blocks >= min_block on the target, with identity = matching bases / block length."""
    blocks = []
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12 or int(f[10]) < min_block:
            continue
        blocks.append({"target_start": int(f[7]) + 1, "target_end": int(f[8]), "strand": f[4],
                       "query_start": int(f[2]) + 1, "query_end": int(f[3]),
                       "identity": round(int(f[9]) / int(f[10]), 4), "block_bp": int(f[10])})
    return sorted(blocks, key=lambda b: b["target_start"])


def diagnostic_windows(status: dict, hit_rows: Sequence[dict], features: dict, reference_offset: int) -> dict:
    """Diagnostic windows per target, in that target's coordinates."""
    deletion, insertion = status["selected_deletion"], status["selected_insertion"]
    ins_sample, ins_contig = parse_source(status["parameters"]["insertion_source"])
    dqa2 = next(r for r in hit_rows if (r["sample_id"], r["contig"], r["family"], r["attributed_locus"])
                == (ins_sample, ins_contig, "DQA", "DQA2"))
    dqa2_start, dqa2_end = features["BOLA-DQA2"]
    dqb_start = features["BOLA-DQB"][0]
    dqa5_end = features["BOLA-DQA5"][1]

    def ref(position: int) -> int:
        return position - reference_offset

    return {
        "deletion_contig": [
            ("deletion_junction", deletion["contig_left"], deletion["contig_right"]),
        ],
        "insertion_contig": [
            ("insertion_left_edge", insertion["contig_left"], insertion["contig_left"]),
            ("insertion_whole_segment", insertion["contig_left"], insertion["contig_right"]),
            ("insertion_right_edge", insertion["contig_right"], insertion["contig_right"]),
            ("insertion_config_DQA2", int(dqa2["contig_start"]), int(dqa2["contig_end"])),
        ],
        "reference": [
            ("reference_insertion_point", ref(insertion["reference_left"]), ref(insertion["reference_right"])),
            ("reference_deletion_point", ref(deletion["reference_left"]), ref(deletion["reference_left"])),
            ("reference_DQA2", ref(dqa2_start), ref(dqa2_end)),
            ("reference_DQB_DQA5", ref(dqb_start), ref(dqa5_end)),
        ],
    }


def classify(present: dict[str, bool]) -> str:
    calls = []
    if present.get("deletion_junction"):
        calls.append("A")
    if all(present.get(d) for d in ("insertion_left_edge", "insertion_whole_segment", "insertion_right_edge")):
        calls.append("B")
    if all(present.get(d) for d in ("reference_deletion_point", "reference_DQA2", "reference_DQB_DQA5")):
        calls.append("reference")
    return "+".join(calls) if calls else "not_classified"


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--attribution", type=Path, required=True, help="06b_dq_attribution")
    parser.add_argument("--assembly-dir", type=Path, required=True, help="06_local_assembly/dq_dr")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--assemblies", type=Path, required=True,
                        help="TSV: assembly_accession, assembly_name, breed, haplotype, source, "
                             "sequence_accession, url, fasta (holding BTA23; the local path is not written out)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--locate", default=f"{CONTIG}:25550001-25750000")
    parser.add_argument("--reference-window", default=f"{CONTIG}:25450001-25800000")
    parser.add_argument("--padding", type=int, default=150000)
    parser.add_argument("--min-block", type=int, default=1000)
    parser.add_argument("--min-window-identity", type=float, default=0.95,
                        help="matching bases over alignment columns within a diagnostic window")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimap2", default="minimap2")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    status = json.loads((args.attribution / "dq_attribution_status.json").read_text(encoding="utf-8"))
    p = status["parameters"]
    with (args.attribution / "dq_family_attribution.tsv").open(encoding="utf-8", newline="") as handle:
        hit_rows = list(csv.DictReader(handle, delimiter="\t"))
    features = load_features(args.gff)
    reference_offset = int(args.reference_window.split(":")[1].split("-")[0]) - 1
    windows = diagnostic_windows(status, hit_rows, features, reference_offset)
    gaps = (p["max_gap_bp"], p["max_gap_sensitivity_bp"])
    primary, strict = f"continuous_gap{p['max_gap_bp']}", f"continuous_gap{p['max_gap_sensitivity_bp']}"
    with args.assemblies.open(encoding="utf-8", newline="") as handle:
        assemblies = list(csv.DictReader(handle, delimiter="\t"))

    region_rows, diagnostic_rows, block_rows, sources = [], [], [], {}
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        targets = {}
        for name, source in (("deletion_contig", p["deletion_source"]), ("insertion_contig", p["insertion_source"])):
            sample, contig = parse_source(source)
            sequence = _run([args.samtools, "faidx", str(args.assembly_dir / sample / "assembly.fasta"), contig]).stdout
            targets[name] = workdir / f"{name}.fa"
            targets[name].write_text(f">{name}\n" + "".join(sequence.splitlines(True)[1:]), encoding="utf-8")
        reference = _run([args.samtools, "faidx", str(args.reference), args.reference_window]).stdout
        targets["reference"] = workdir / "reference.fa"
        targets["reference"].write_text(">reference\n" + "".join(reference.splitlines(True)[1:]), encoding="utf-8")
        locate = workdir / "locate.fa"
        locate.write_text(_run([args.samtools, "faidx", str(args.reference), args.locate]).stdout, encoding="utf-8")

        for assembly in assemblies:
            fasta = Path(assembly["fasta"])
            sources[f"{assembly['assembly_accession']}:{assembly['sequence_accession']}"] = sha256_file(fasta)
            chromosome = read_fasta(fasta, assembly["sequence_accession"])
            whole = workdir / "chromosome.fa"
            whole.write_text(f">chr\n{chromosome}\n", encoding="utf-8")
            span = locate_span(_run([args.minimap2, "-x", "asm20", "-t", str(args.threads), str(whole),
                                     str(locate)]).stdout, 5000)
            label = {k: assembly[k] for k in ("assembly_accession", "assembly_name", "breed", "haplotype",
                                              "sequence_accession")}
            if span is None:
                region_rows.append({**label, "region_start": "", "region_end": "", "strand": ""})
                continue
            low, high = max(span[0] - args.padding, 1), min(span[1] + args.padding, len(chromosome))
            region = workdir / "region.fa"
            region.write_text(f">region\n{chromosome[low - 1:high]}\n", encoding="utf-8")
            region_rows.append({**label, "region_start": low, "region_end": high, "strand": span[2]})
            for name, target in targets.items():
                sam = _run([args.minimap2, "-a", "-x", "asm20", "--eqx", "-t", str(args.threads), str(target),
                            str(region)]).stdout
                # one region aligns in pieces: score supplementary pieces like primary ones
                alignments = [replace(a, flag=a.flag & ~0x800) for a in parse_sam_alignments(sam)
                              if not a.flag & 0x100]
                for diagnostic, left, right in windows[name]:
                    support = supporting_reads(alignments, name, left, right, p["spanning_flank_bp"], gaps,
                                               p["min_aligned_fraction"], 0, p["boundary_max_gap_bp"],
                                               p["boundary_margin_bp"])
                    similar = [r for r in support
                               if r["window_identity"] != "" and r["window_identity"] >= args.min_window_identity]
                    diagnostic_rows.append({
                        **label, "target": name, "diagnostic": diagnostic, "window_left": left - p["spanning_flank_bp"],
                        "window_right": right + p["spanning_flank_bp"],
                        "enclosing_alignments": len(support),
                        "best_window_identity": max((r["window_identity"] for r in support), default=""),
                        "present": any(r[primary] for r in similar),
                        "present_strict": any(r[strict] for r in similar),
                        "longest_deletion": min((r["longest_deletion"] for r in support), default=""),
                        "longest_insertion": min((r["longest_insertion"] for r in support), default=""),
                    })
                paf = _run([args.minimap2, "-c", "-x", "asm20", "-t", str(args.threads), str(target),
                            str(region)]).stdout
                for block in paf_blocks(paf, args.min_block):
                    block_rows.append({**label, "target": name, **block,
                                       "region_start": block["query_start"] + low - 1,
                                       "region_end": block["query_end"] + low - 1})

    classification_rows = []
    for region in region_rows:
        present = {r["diagnostic"]: r["present"] for r in diagnostic_rows
                   if (r["assembly_accession"], r["sequence_accession"])
                   == (region["assembly_accession"], region["sequence_accession"])}
        classification_rows.append({**region, "configuration": classify(present) if present else "region_not_found"})
    manifest_fields = ["assembly_accession", "assembly_name", "breed", "haplotype", "source", "sequence_accession",
                       "url"]
    args.output.mkdir(parents=True, exist_ok=True)
    label_fields = ["assembly_accession", "assembly_name", "breed", "haplotype", "sequence_accession"]
    tables = {
        "public_dq_regions.tsv": (region_rows, label_fields + ["region_start", "region_end", "strand"]),
        "public_dq_configurations.tsv": (diagnostic_rows, label_fields + [
            "target", "diagnostic", "window_left", "window_right", "enclosing_alignments", "best_window_identity",
            "present", "present_strict", "longest_deletion", "longest_insertion"]),
        "public_dq_classification.tsv": (classification_rows, label_fields + ["configuration"]),
        "public_dq_assemblies.tsv": ([{k: a.get(k, "") for k in manifest_fields} for a in assemblies],
                                     manifest_fields),
        "public_dq_synteny_blocks.tsv": (block_rows, label_fields + [
            "target", "target_start", "target_end", "strand", "region_start", "region_end", "identity",
            "block_bp"]),
    }
    for filename, (rows, fields) in tables.items():
        with (args.output / filename).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    parameters = {
        "locate": args.locate, "reference_window": args.reference_window, "padding_bp": args.padding,
        "alignment": "minimap2 -x asm20 (region as query; primary and supplementary pieces)",
        "support_rule": p["support_rule"], "flank_bp": p["spanning_flank_bp"],
        "targets": {"deletion_contig": p["deletion_source"], "insertion_contig": p["insertion_source"],
                    "reference": args.reference_window},
        "min_synteny_block_bp": args.min_block,
        "min_window_identity": args.min_window_identity,
        "classification": "A: deletion_junction; B: insertion left edge + whole segment + right edge; "
                          "reference: reference deletion point + DQA2 + BOLA-DQB-DQA5; else not_classified",
    }
    (args.output / "public_dq_status.json").write_text(
        json.dumps({"status": "passed", "assemblies": len(assemblies), "parameters": parameters}, indent=2,
                   sort_keys=True) + "\n", encoding="utf-8")
    write_provenance(args.output, {
        "inputs": {"assembly_chromosome_sha256": sources,
                   "attribution_status_sha256": sha256_file(args.attribution / "dq_attribution_status.json"),
                   "reference_sha256": sha256_file(args.reference)},
        "parameters": parameters,
    })
    print(f"public DQ comparison written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
