"""Generate validated analytical targets from a fixed RefSeq GFF."""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from dataclasses import dataclass, replace
from dataclasses import fields as dataclass_fields
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from urllib.parse import unquote


class CoordinateError(ValueError):
    """Raised when a coordinate cannot be resolved without guessing."""


@dataclass(frozen=True)
class GeneFeature:
    seqid: str
    start: int
    end: int
    strand: str
    feature_id: str
    gene_symbol: str


@dataclass(frozen=True)
class TargetSpec:
    name: str
    gene_symbols: tuple[str, ...]
    category: str
    absent_policy: str = "error"
    probe_id: str = ""
    probe_loci: tuple[str, ...] = ()


@dataclass(frozen=True)
class TargetRecord:
    assembly_accession: str
    annotation_release: str
    contig_accession: str
    display_alias: str
    feature_id: str
    gene_symbol: str
    target_name: str
    category: str
    target_type: str
    gff_start: int | None
    gff_end: int | None
    bed_start: int | None
    bed_end: int | None
    strand: str
    sequence_sha256: str
    probe_id: str
    probe_loci: str
    probe_alignment_status: str


def gff_to_bed(start: int, end: int) -> tuple[int, int]:
    """Convert a 1-based closed GFF interval to 0-based half-open BED."""
    if start < 1 or end < start:
        raise CoordinateError(f"invalid GFF interval: {start}-{end}")
    return start - 1, end


def interval_length(start: int, end: int) -> int:
    """Return the length of a 0-based half-open interval."""
    if start < 0 or end < start:
        raise CoordinateError(f"invalid BED interval: {start}-{end}")
    return end - start


def resolve_contig(name: str, aliases: Mapping[str, Sequence[str]]) -> str:
    """Resolve an alias only when it identifies one canonical accession."""
    candidates: set[str] = set()
    if name in aliases:
        candidates.update(str(value) for value in aliases[name])
    for values in aliases.values():
        normalized = [str(value) for value in values]
        if name in normalized:
            candidates.add(name)
    if not candidates:
        raise CoordinateError(f"unknown contig alias: {name}")
    if len(candidates) != 1:
        raise CoordinateError(f"ambiguous contig alias {name!r}: {sorted(candidates)}")
    return next(iter(candidates))


def parse_gene_features(gff_path: Path) -> list[GeneFeature]:
    """Read gene features and the RefSeq ID/gene attributes from GFF3."""
    features: list[GeneFeature] = []
    with Path(gff_path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise CoordinateError(f"invalid GFF row at line {line_number}")
            if fields[2] != "gene":
                continue
            attributes = _parse_attributes(fields[8])
            feature_id = attributes.get("ID", "")
            gene_symbol = attributes.get("gene") or attributes.get("Name", "")
            if not feature_id or not gene_symbol:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise CoordinateError(f"invalid GFF coordinates at line {line_number}") from exc
            gff_to_bed(start, end)
            if fields[6] not in {"+", "-", ".", "?"}:
                raise CoordinateError(f"invalid GFF strand at line {line_number}")
            features.append(
                GeneFeature(
                    seqid=fields[0],
                    start=start,
                    end=end,
                    strand=fields[6],
                    feature_id=feature_id,
                    gene_symbol=gene_symbol,
                )
            )
    return features


def load_target_specs(path: Path) -> list[TargetSpec]:
    """Load a coordinate-free target specification."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        forbidden = {
            "chrom",
            "chromosome",
            "contig",
            "seqid",
            "start",
            "end",
            "position",
            "bed_start",
            "bed_end",
            "gff_start",
            "gff_end",
        }
        if fieldnames & forbidden:
            raise CoordinateError("target specification must not contain coordinate columns")
        required = {"name", "gene_symbols", "category", "absent_policy", "probe_id"}
        if not required.issubset(fieldnames):
            missing = ", ".join(sorted(required - fieldnames))
            raise CoordinateError(f"target specification is missing columns: {missing}")
        specs: list[TargetSpec] = []
        names: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            symbols = tuple(
                value.strip() for value in row["gene_symbols"].split("|") if value.strip()
            )
            name = row["name"].strip()
            if not name or not symbols or not row["category"].strip():
                raise CoordinateError(f"incomplete target specification at line {line_number}")
            if name in names:
                raise CoordinateError(f"duplicate target name at line {line_number}: {name}")
            names.add(name)
            absent_policy = row["absent_policy"].strip() or "error"
            if absent_policy not in {"error", "probe_only"}:
                raise CoordinateError(f"invalid absent policy at line {line_number}")
            probe_id = (row.get("probe_id") or "").strip()
            if probe_id == "-":
                probe_id = ""
            probe_loci = tuple(
                value.strip()
                for value in (row.get("probe_loci") or "").split("|")
                if value.strip() and value.strip() != "-"
            )
            if absent_policy == "probe_only" and not probe_id:
                raise CoordinateError(f"probe-only target lacks a probe at line {line_number}")
            specs.append(
                TargetSpec(
                    name=name,
                    gene_symbols=symbols,
                    category=row["category"].strip(),
                    absent_policy=absent_policy,
                    probe_id=probe_id,
                    probe_loci=probe_loci,
                )
            )
    if not specs:
        raise CoordinateError("target specification is empty")
    return specs


def load_fasta_sequences(path: Path, names: set[str] | None = None) -> dict[str, str]:
    """Stream a FASTA and retain only requested sequence records."""
    requested = None if names is None else set(names)
    sequences: dict[str, str] = {}
    current_name = ""
    chunks: list[str] = []
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                _store_fasta_record(sequences, current_name, chunks, requested)
                current_name = line[1:].split(maxsplit=1)[0]
                if not current_name:
                    raise CoordinateError(f"empty FASTA name at line {line_number}")
                chunks = []
            elif not current_name:
                raise CoordinateError(f"FASTA sequence precedes header at line {line_number}")
            elif requested is None or current_name in requested:
                chunks.append(line)
    _store_fasta_record(sequences, current_name, chunks, requested)
    if requested is not None:
        missing = requested - sequences.keys()
        if missing:
            raise CoordinateError(f"reference FASTA is missing sequences: {', '.join(sorted(missing))}")
    return sequences


#: Flank added to every gene interval when building the mapping retention
#: filter. Roughly twice the cohort's mean read length, so a read anchored just
#: outside a gene, and the supplementary alignment of a read crossing a
#: breakpoint, are still retained.
#: The analytical gene intervals; per-gene depth is measured over these.
ANALYTICAL_BED_NAME = "targets.bed"
#: The wider retention filter mapping runs with, and the file the reference
#: seal must cover. Defined once so a caller cannot seal the other one.
MAPPING_BED_NAME = "mapping_targets.bed"

MAPPING_FLANK = 50_000

#: Categories whose members are joined into one contiguous interval. Gene
#: bodies alone leave most of each cluster unmapped: measured on the pilot,
#: 30% of the class I span and 24% of the class II span had zero coverage,
#: which makes local assembly and same-contig determination impossible.
CONTIGUOUS_CATEGORIES = frozenset(
    {"class_I", "class_II", "class_III", "marker_host_gene", "btn_btnl"}
)

#: Categories joined into a single span rather than one span each, because they
#: interleave along the same stretch of BTA23.
_MHC_GROUP = frozenset({"class_I", "class_II", "class_III", "marker_host_gene"})


def _field(row: object, name: str) -> object:
    if isinstance(row, Mapping):
        return row.get(name)
    return getattr(row, name, None)


def build_mapping_intervals(
    records: Iterable[object],
    flank: int,
    contiguous_categories: Iterable[str],
    contig_lengths: Mapping[str, int],
) -> list[tuple[str, int, int, str]]:
    """Build the merged retention filter used to keep alignments during mapping.

    This is deliberately separate from the analytical gene intervals: widening
    those would silently change every per-gene depth measurement.
    """
    if flank < 0:
        raise CoordinateError("mapping flank must not be negative")
    contiguous = set(contiguous_categories)
    rows = [row for row in records if _field(row, "target_type") == "reference_feature"]

    raw: dict[str, list[tuple[int, int, str]]] = {}
    groups: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for row in rows:
        contig = str(_field(row, "contig_accession"))
        start, end = _field(row, "bed_start"), _field(row, "bed_end")
        if start is None or end is None:
            raise CoordinateError(f"reference feature without coordinates: {_field(row, 'target_name')}")
        if contig not in contig_lengths:
            raise CoordinateError(f"contig length is unknown for {contig}")
        limit = int(contig_lengths[contig])
        raw.setdefault(contig, []).append(
            (max(0, int(start) - flank), min(limit, int(end) + flank), str(_field(row, "target_name")))
        )
        category = str(_field(row, "category"))
        if category in contiguous:
            key = (contig, "MHC" if category in _MHC_GROUP else category)
            groups.setdefault(key, []).append((int(start), int(end)))

    for (contig, label), spans in groups.items():
        limit = int(contig_lengths[contig])
        raw[contig].append(
            (
                max(0, min(start for start, _ in spans) - flank),
                min(limit, max(end for _, end in spans) + flank),
                f"span_{label}",
            )
        )

    merged: list[tuple[str, int, int, str]] = []
    for contig in sorted(raw):
        current: tuple[int, int, list[str]] | None = None
        for start, end, name in sorted(raw[contig]):
            if current is not None and start <= current[1]:
                current = (current[0], max(current[1], end), [*current[2], name])
            else:
                if current is not None:
                    merged.append((contig, current[0], current[1], "+".join(current[2])))
                current = (start, end, [name])
        if current is not None:
            merged.append((contig, current[0], current[1], "+".join(current[2])))
    return merged


def write_coordinate_outputs(
    records: Iterable[TargetRecord],
    *,
    sequences: Mapping[str, str],
    output_dir: Path,
) -> dict[str, Path]:
    """Write the complete registry and only validated reference targets to BED/FASTA.

    Two BED files are produced and they are not interchangeable.
    ``targets.bed`` holds the analytical gene intervals and is what per-gene
    depth is measured over. ``mapping_targets.bed`` is the wider retention
    filter used while mapping; widening the analytical intervals instead would
    silently change every depth measurement.
    """
    rows = list(records)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    registry_path = destination / "coordinate_registry.tsv"
    bed_path = destination / ANALYTICAL_BED_NAME
    mapping_bed_path = destination / MAPPING_BED_NAME
    fasta_path = destination / "target_sequences.fa"
    registry_partial = registry_path.with_suffix(f"{registry_path.suffix}.partial")
    bed_partial = bed_path.with_suffix(f"{bed_path.suffix}.partial")
    fasta_partial = fasta_path.with_suffix(f"{fasta_path.suffix}.partial")
    columns = [field.name for field in dataclass_fields(TargetRecord)]
    with registry_partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: getattr(row, column) for column in columns})
    with bed_partial.open("w", encoding="utf-8", newline="") as bed_handle, fasta_partial.open(
        "w", encoding="ascii", newline=""
    ) as fasta_handle:
        for row in rows:
            if row.target_type != "reference_feature":
                continue
            assert row.bed_start is not None and row.bed_end is not None
            sequence = sequences.get(row.contig_accession)
            if sequence is None:
                raise CoordinateError(f"reference sequence is absent for {row.contig_accession}")
            extracted = sequence[row.bed_start : row.bed_end].upper()
            if sha256(extracted.encode("ascii")).hexdigest() != row.sequence_sha256:
                raise CoordinateError(f"sequence checksum changed for {row.target_name}")
            bed_handle.write(
                f"{row.contig_accession}\t{row.bed_start}\t{row.bed_end}\t{row.target_name}\n"
            )
            fasta_handle.write(
                f">{row.target_name}|{row.contig_accession}:{row.gff_start}-{row.gff_end}\n"
                f"{extracted}\n"
            )
    mapping_partial = mapping_bed_path.with_suffix(f"{mapping_bed_path.suffix}.partial")
    intervals = build_mapping_intervals(
        rows,
        MAPPING_FLANK,
        CONTIGUOUS_CATEGORIES,
        {contig: len(sequence) for contig, sequence in sequences.items()},
    )
    with mapping_partial.open("w", encoding="utf-8", newline="") as handle:
        for contig, start, end, name in intervals:
            handle.write(f"{contig}\t{start}\t{end}\t{name}\n")

    registry_partial.replace(registry_path)
    bed_partial.replace(bed_path)
    mapping_partial.replace(mapping_bed_path)
    fasta_partial.replace(fasta_path)
    return {
        "registry": registry_path,
        "bed": bed_path,
        "mapping_bed": mapping_bed_path,
        "fasta": fasta_path,
    }


def apply_probe_alignment_statuses(
    records: Iterable[TargetRecord],
    locus_statuses: Mapping[str, str],
) -> list[TargetRecord]:
    """Attach panel-level evidence without converting probes into reference genes."""
    rows = list(records)
    locus_usage = Counter(
        locus.casefold()
        for record in rows
        for locus in record.probe_loci.split("|")
        if locus
    )
    updated: list[TargetRecord] = []
    for record in rows:
        if not record.probe_id:
            updated.append(record)
            continue
        expected_loci = {
            value.casefold() for value in record.probe_loci.split("|") if value
        }
        matching = [
            status
            for locus, status in locus_statuses.items()
            if locus.casefold() in expected_loci
        ]
        aligned_loci = {
            locus.casefold()
            for locus, status in locus_statuses.items()
            if status == "aligned" and locus.casefold() in expected_loci
        }
        if aligned_loci:
            probe_status = (
                "panel_family_aligned"
                if any(locus_usage[locus] > 1 for locus in aligned_loci)
                else "panel_aligned"
            )
        elif matching:
            probe_status = "panel_no_alignment"
        else:
            probe_status = "probe_unavailable"
        updated.append(replace(record, probe_alignment_status=probe_status))
    return updated


def generate_target_records(
    features: Iterable[GeneFeature],
    target_specs: Iterable[TargetSpec],
    *,
    assembly_accession: str,
    annotation_release: str,
    aliases: Mapping[str, Sequence[str]],
    sequences: Mapping[str, str],
) -> list[TargetRecord]:
    """Resolve target specifications without importing historical coordinates."""
    feature_rows = list(features)
    records: list[TargetRecord] = []
    for spec in target_specs:
        matches = [row for row in feature_rows if row.gene_symbol in spec.gene_symbols]
        if not matches:
            if spec.absent_policy != "probe_only" or not spec.probe_id:
                raise CoordinateError(f"target feature is absent: {spec.name}")
            records.append(_probe_only_record(spec, assembly_accession, annotation_release))
            continue
        if len(matches) != 1:
            raise CoordinateError(f"non-unique target feature for {spec.name}: {len(matches)}")

        feature = matches[0]
        contig = resolve_contig(feature.seqid, aliases)
        sequence = sequences.get(contig)
        if sequence is None:
            raise CoordinateError(f"reference sequence is absent for {contig}")
        bed_start, bed_end = gff_to_bed(feature.start, feature.end)
        if bed_end > len(sequence):
            raise CoordinateError(f"target {spec.name} extends beyond {contig}")
        target_sequence = sequence[bed_start:bed_end].upper()
        if len(target_sequence) != interval_length(bed_start, bed_end):
            raise CoordinateError(f"sequence extraction failed for {spec.name}")
        records.append(
            TargetRecord(
                assembly_accession=assembly_accession,
                annotation_release=annotation_release,
                contig_accession=contig,
                display_alias=_display_alias(contig, aliases),
                feature_id=feature.feature_id,
                gene_symbol=feature.gene_symbol,
                target_name=spec.name,
                category=spec.category,
                target_type="reference_feature",
                gff_start=feature.start,
                gff_end=feature.end,
                bed_start=bed_start,
                bed_end=bed_end,
                strand=feature.strand,
                sequence_sha256=sha256(target_sequence.encode("ascii")).hexdigest(),
                probe_id=spec.probe_id,
                probe_loci="|".join(spec.probe_loci),
                probe_alignment_status="pending" if spec.probe_id else "not_available",
            )
        )
    return records


def _parse_attributes(raw: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in raw.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            parsed[unquote(key)] = unquote(value)
    return parsed


def _store_fasta_record(
    sequences: dict[str, str],
    name: str,
    chunks: list[str],
    requested: set[str] | None,
) -> None:
    if not name or (requested is not None and name not in requested):
        return
    if name in sequences:
        raise CoordinateError(f"duplicate FASTA sequence: {name}")
    sequence = "".join(chunks).upper()
    if not sequence:
        raise CoordinateError(f"empty FASTA sequence: {name}")
    sequences[name] = sequence


def _display_alias(contig: str, aliases: Mapping[str, Sequence[str]]) -> str:
    matching = [alias for alias, values in aliases.items() if contig in values]
    numeric = [alias for alias in matching if alias.isdigit()]
    if len(numeric) == 1:
        return numeric[0]
    non_accession = [alias for alias in matching if alias != contig and not alias.startswith("chr")]
    if len(non_accession) == 1:
        return non_accession[0]
    if len(matching) == 1:
        return matching[0]
    raise CoordinateError(f"contig {contig!r} does not have one preferred display alias")


def _probe_only_record(
    spec: TargetSpec,
    assembly_accession: str,
    annotation_release: str,
) -> TargetRecord:
    return TargetRecord(
        assembly_accession=assembly_accession,
        annotation_release=annotation_release,
        contig_accession="",
        display_alias="",
        feature_id="",
        gene_symbol=spec.name,
        target_name=spec.name,
        category=spec.category,
        target_type="haplotype_probe",
        gff_start=None,
        gff_end=None,
        bed_start=None,
        bed_end=None,
        strand="",
        sequence_sha256="",
        probe_id=spec.probe_id,
        probe_loci="|".join(spec.probe_loci),
        probe_alignment_status="pending",
    )
