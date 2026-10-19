"""Build Additional file 2: the DQ configuration sequences and read-support records.

The manuscript's DQ result rests on two contigs, a three-sequence panel built
from them, diagnostic coordinates on that panel, and per-read support
records. This script collects them into one archive for submission with the
manuscript (never for git: it contains FASTA, see
docs/REPOSITORY_DATA_POLICY.md):

- `dq_configuration_panel.fasta`: the panel exactly as stage 06b builds it;
- `dq_representative_contigs.fasta`: the full configuration A (ID11) and
  configuration B (ID06) contigs;
- the 06b diagnostics, contig events, support comparison and per-read
  support table, and the 06c junction comparison;
- the support comparison and per-read table of the panel-sensitivity run,
  in which another animal's configuration A contig replaced the ID11 contig;
- `parameters.json` and `SHA256SUMS.txt`.

Sample names are replaced by the manuscript's anonymized IDs (ID01-ID13)
from the internal anonymization key, and the script refuses to write an
archive in which an original sample name survives.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import tempfile
import zipfile
from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.attribute_dq_contigs import (
    _run,
    _select_dqa2_hit,
    build_panel,
    load_features,
    parse_source,
)

KEY_ROW = re.compile(r"^\|\s*(ID\d{2})\s*\|\s*([^|]+?)\s*\|\s*(Site [A-D])\s*\|")
FASTA_WIDTH = 80
README = """Additional file 2. DQ configuration sequences and read-support records

Coordinates are 1-based and closed. ARS-UCD2.0 = RefSeq GCF_002263795.3 (BTA23 = NC_037350.1).
Animals are identified by the anonymized IDs used in the manuscript (ID01-ID13).

dq_configuration_panel.fasta
    The three-sequence panel to which every animal's reads were mapped competitively:
    configuration A (ID11 contig across the junction), configuration B (ID06 contig across the
    25.8-kb non-reference segment and its DQA2-placed sequence), and ARS-UCD2.0 over the same
    interval. Headers give the source contig or chromosome interval.
dq_representative_contigs.fasta
    The full ID11 (configuration A) and ID06 (configuration B) Flye contigs (R9 consensus).
dq_panel_diagnostics.tsv
    Diagnostic windows on the panel (target, left, right, meaning).
dq_contig_events.tsv
    Junctions (reference sequence absent from a contig) and non-reference segments on every
    animal's DQ/DR contigs, with own-read support at the primary rule, the 50-bp rule, and by
    alignment end points.
dq_configuration_support_comparison.tsv
    Panel read counts per animal and diagnostic: end-point count, primary rule (gap <= 200 bp,
    <= 50 bp within 200 bp of a diagnostic boundary), and strict rule (gap <= 50 bp).
dq_support_reads.tsv
    Every read alignment whose end points enclosed a diagnostic window: read identifier (as in
    the deposited FASTQ), strand, MAPQ, aligned bases in each flank and in the diagnostic
    interval, longest deletion/insertion, longest gap near a boundary, the CIGAR clipped to the
    window, and whether the read met each rule.
dq_junction_structure.tsv, dq_junction_pairwise.tsv, dq_junction_difference_bins.tsv
    Configuration A junction compared across the five contigs that carry it: anchors, retained
    reference sequence, breakpoint overlap (repeat / microhomology), and pairwise identity.
dq_panel_sensitivity_support_comparison.tsv, dq_panel_sensitivity_support_reads.tsv
    The same counts and per-read records with the ID13 configuration A contig, which lacks the 91 bp
    carried by the ID11 contig, in place of the ID11 contig in the panel.
parameters.json
    Thresholds and sources used to build the panel and count support.
"""


def anonymization_key(path: Path) -> dict[str, str]:
    key = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = KEY_ROW.match(line)
        if match:
            key[match.group(2)] = match.group(1)
    if len(key) != 13:
        raise ValueError(f"expected 13 anonymization rows in {path}, found {len(key)}")
    return key


def anonymize(text: str, key: dict[str, str]) -> str:
    for original in sorted(key, key=len, reverse=True):
        text = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(original)}(?![A-Za-z0-9_])", key[original], text)
    leaked = [original for original in key if original in text]
    if leaked:
        raise ValueError(f"original sample names survive anonymization: {', '.join(leaked)}")
    return text


def fasta(records: Sequence[tuple[str, str]]) -> str:
    out = []
    for name, sequence in records:
        out.append(f">{name}")
        out.extend(sequence[i:i + FASTA_WIDTH] for i in range(0, len(sequence), FASTA_WIDTH))
    return "\n".join(out) + "\n"


def sequence_of(text: str) -> str:
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--attribution", type=Path, required=True, help="06b_dq_attribution")
    parser.add_argument("--junctions", type=Path, required=True, help="06c_dq_junction_comparison")
    parser.add_argument("--sensitivity", type=Path, required=True,
                        help="06b run with another configuration A contig in the panel")
    parser.add_argument("--assembly-dir", type=Path, required=True, help="06_local_assembly/dq_dr")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--key", type=Path, required=True, help="docs/ANONYMIZATION_KEY_INTERNAL.md")
    parser.add_argument("--output", type=Path, required=True, help="the .zip to write")
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    key = anonymization_key(args.key)
    status = json.loads((args.attribution / "dq_attribution_status.json").read_text(encoding="utf-8"))
    p = status["parameters"]
    with (args.attribution / "dq_family_attribution.tsv").open(encoding="utf-8", newline="") as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    panel_args = SimpleNamespace(deletion_source=p["deletion_source"], insertion_source=p["insertion_source"],
                                 assembly_dir=args.assembly_dir, reference=args.reference, samtools=args.samtools)
    files: dict[str, str] = {}
    with tempfile.TemporaryDirectory() as tmp:
        panel, diagnostics = build_panel(panel_args, load_features(args.gff), status["selected_deletion"],
                                         status["selected_insertion"], _select_dqa2_hit(hits, p["insertion_source"]),
                                         Path(tmp))
        panel_records = []
        name = None
        for line in panel.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                name = line[1:]
            else:
                panel_records.append((name, line))
    sources = {
        "deletion_config": ("configuration_A_junction", p["deletion_source"]),
        "insertion_config": ("configuration_B_segment_DQA2", p["insertion_source"]),
    }
    contigs = []
    for target, (label, source) in sources.items():
        sample, contig = parse_source(source)
        sequence = sequence_of(_run([args.samtools, "faidx", str(args.assembly_dir / sample / "assembly.fasta"),
                                     contig]).stdout)
        contigs.append((f"{label} {key[sample]}_{contig} length={len(sequence)}", sequence))
    named_panel = []
    for target, sequence in panel_records:
        if target in sources:
            label, source = sources[target]
            sample, contig = parse_source(source)
            named_panel.append((f"{target} {label} from {key[sample]}_{contig} length={len(sequence)}", sequence))
        else:
            named_panel.append((f"{target} ARS-UCD2.0 NC_037350.1 interval length={len(sequence)}", sequence))
    files["dq_configuration_panel.fasta"] = fasta(named_panel)
    files["dq_representative_contigs.fasta"] = fasta(contigs)
    copies = {
        "dq_panel_diagnostics.tsv": args.attribution / "dq_configuration_diagnostics.tsv",
        "dq_contig_events.tsv": args.attribution / "dq_contig_events.tsv",
        "dq_configuration_support_comparison.tsv": args.attribution / "dq_configuration_support_comparison.tsv",
        "dq_support_reads.tsv": args.attribution / "dq_support_reads.tsv",
        "dq_junction_structure.tsv": args.junctions / "dq_junction_structure.tsv",
        "dq_junction_pairwise.tsv": args.junctions / "dq_junction_pairwise.tsv",
        "dq_junction_difference_bins.tsv": args.junctions / "dq_junction_difference_bins.tsv",
        "dq_panel_sensitivity_support_comparison.tsv":
            args.sensitivity / "dq_configuration_support_comparison.tsv",
        "dq_panel_sensitivity_support_reads.tsv": args.sensitivity / "dq_support_reads.tsv",
    }
    for name, path in copies.items():
        files[name] = path.read_text(encoding="utf-8")
    files["parameters.json"] = json.dumps(
        {"attribution": p,
         "panel_sensitivity": json.loads(
             (args.sensitivity / "dq_attribution_status.json").read_text(encoding="utf-8"))["parameters"],
         "junction_comparison": json.loads(
            (args.junctions / "dq_junction_status.json").read_text(encoding="utf-8"))["parameters"]},
        indent=2, sort_keys=True) + "\n"
    files["README.txt"] = README
    files = {name: anonymize(text, key) for name, text in files.items()}
    files["SHA256SUMS.txt"] = "".join(
        f"{hashlib.sha256(text.encode('utf-8')).hexdigest()}  {name}\n" for name, text in sorted(files.items()))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(files):
            info = zipfile.ZipInfo(f"Additional_file_2/{name}", date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, files[name])
    args.output.write_bytes(buffer.getvalue())
    print(f"{args.output}\tsha256={hashlib.sha256(buffer.getvalue()).hexdigest()}\tbytes={len(buffer.getvalue())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
