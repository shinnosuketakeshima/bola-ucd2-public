"""Freeze and validate external reference resources for the reanalysis."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
from datetime import date
from dataclasses import dataclass
from hashlib import md5, sha256
from pathlib import Path, PurePosixPath
from tempfile import NamedTemporaryFile
from typing import Callable, Iterable, Mapping
from urllib.parse import urlparse
from urllib.request import urlopen

import yaml

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.coordinate_registry import (
    apply_probe_alignment_statuses,
    generate_target_records,
    load_fasta_sequences,
    load_target_specs,
    parse_gene_features,
    write_coordinate_outputs,
)
from scripts.codex_reanalysis.provenance import write_provenance


class ReferenceRegistryError(RuntimeError):
    """Raised when a reference package cannot be verified safely."""


@dataclass(frozen=True)
class SequenceComparison:
    contig: str
    old_accession: str
    new_accession: str
    old_length: int
    new_length: int
    old_sha256: str
    new_sha256: str
    status: str


CommandRunner = Callable[[list[str]], None]
AlignmentRunner = Callable[[list[str], Path], None]


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_tool_version(output: str, expected: str, tool: str) -> str:
    """Require the configured pinned version as an exact first-line token."""
    first_line = output.strip().splitlines()[0] if output.strip() else ""
    if expected not in first_line.replace(",", " ").split():
        raise ReferenceRegistryError(
            f"{tool} version mismatch: expected {expected!r}, observed {first_line!r}"
        )
    return first_line


def parse_official_checksums(manifest_path: Path) -> dict[str, str]:
    """Parse NCBI md5checksums.txt while rejecting escaping paths."""
    checksums: dict[str, str] = {}
    with Path(manifest_path).open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split(maxsplit=1)
            if len(parts) != 2 or len(parts[0]) != 32:
                raise ReferenceRegistryError(f"invalid checksum row at line {line_number}")
            checksum, raw_name = parts
            normalized = raw_name.removeprefix("*").removeprefix("./")
            package_path = PurePosixPath(normalized)
            if package_path.is_absolute() or ".." in package_path.parts or not normalized:
                raise ReferenceRegistryError(f"unsafe checksum path: {raw_name}")
            if any(character not in "0123456789abcdefABCDEF" for character in checksum):
                raise ReferenceRegistryError(f"invalid checksum at line {line_number}")
            checksums[package_path.as_posix()] = checksum.lower()
    return checksums


def parse_assembly_report_aliases(report_path: Path) -> dict[str, list[str]]:
    """Build unambiguous chromosome aliases from an NCBI assembly report."""
    aliases: dict[str, list[str]] = {}
    header: list[str] | None = None
    with Path(report_path).open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith("# Sequence-Name\t"):
                header = line.removeprefix("# ").split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise ReferenceRegistryError("assembly report header is absent")
            values = line.split("\t")
            if len(values) != len(header):
                raise ReferenceRegistryError("invalid assembly report row")
            row = dict(zip(header, values, strict=True))
            if row["Sequence-Role"] != "assembled-molecule":
                continue
            accession = row["RefSeq-Accn"]
            if accession in {"", "na"}:
                continue
            for alias in (
                row["Sequence-Name"],
                row["Assigned-Molecule"],
                row["UCSC-style-name"],
                accession,
            ):
                if alias in {"", "na"}:
                    continue
                existing = aliases.setdefault(alias, [])
                if accession not in existing:
                    existing.append(accession)
                if len(existing) > 1:
                    raise ReferenceRegistryError(
                        f"ambiguous chromosome alias {alias!r}: {existing}"
                    )
    if not aliases:
        raise ReferenceRegistryError("assembly report has no assembled chromosomes")
    return aliases


def verify_official_checksums(
    package_root: Path,
    expected: Mapping[str, str],
    required_names: Iterable[str],
) -> dict[str, str]:
    verified: dict[str, str] = {}
    root = Path(package_root).resolve()
    for raw_name in required_names:
        name = PurePosixPath(raw_name).as_posix()
        if name not in expected:
            raise ReferenceRegistryError(f"official checksum is absent for {name}")
        path = (root / Path(*PurePosixPath(name).parts)).resolve()
        if path != root and root not in path.parents:
            raise ReferenceRegistryError(f"unsafe package path: {name}")
        if not path.is_file():
            raise ReferenceRegistryError(f"reference package file is absent: {name}")
        observed = _md5_file(path)
        if observed != expected[name].lower():
            raise ReferenceRegistryError(
                f"checksum mismatch for {name}: expected {expected[name]}, observed {observed}"
            )
        verified[name] = observed
    return verified


def build_index_commands(
    samtools: str,
    minimap2: str,
    fasta_path: str | Path,
) -> list[list[str]]:
    fasta = str(fasta_path)
    return [
        [samtools, "faidx", "-o", f"{fasta}.fai.partial", fasta],
        [minimap2, "-d", f"{fasta}.mmi.partial", fasta],
    ]


def compare_sequence_records(
    *,
    contig: str,
    old_accession: str,
    old_sequence: str,
    new_accession: str,
    new_sequence: str,
) -> SequenceComparison:
    old_bytes = old_sequence.upper().encode("ascii")
    new_bytes = new_sequence.upper().encode("ascii")
    old_digest = sha256(old_bytes).hexdigest()
    new_digest = sha256(new_bytes).hexdigest()
    identical = (
        old_accession == new_accession
        and len(old_bytes) == len(new_bytes)
        and old_digest == new_digest
    )
    return SequenceComparison(
        contig=contig,
        old_accession=old_accession,
        new_accession=new_accession,
        old_length=len(old_bytes),
        new_length=len(new_bytes),
        old_sha256=old_digest,
        new_sha256=new_digest,
        status="identical" if identical else "decision_required",
    )


def compare_fasta_contigs(
    *,
    old_fasta: Path,
    new_fasta: Path,
    old_aliases: Mapping[str, list[str]],
    aliases: Mapping[str, list[str]],
    contigs: Iterable[str],
    output_dir: Path,
    minimap2: str,
    alignment_runner: AlignmentRunner | None = None,
) -> list[dict[str, object]]:
    """Compare named assembly records and retain an alignment for every record."""
    requested_aliases = [str(contig) for contig in contigs]
    accessions: dict[str, str] = {}
    old_accessions: dict[str, str] = {}
    old_fasta_names: dict[str, str] = {}
    for alias in requested_aliases:
        values = aliases.get(alias, [])
        if len(values) != 1:
            raise ReferenceRegistryError(f"comparison contig {alias!r} is not unambiguous")
        accessions[alias] = values[0]
        old_values = old_aliases.get(alias, [])
        if len(old_values) != 1:
            raise ReferenceRegistryError(
                f"old comparison contig {alias!r} is not unambiguous"
            )
        old_accessions[alias] = old_values[0]
        ucsc_alias = f"chr{alias}"
        old_fasta_names[alias] = (
            ucsc_alias if old_aliases.get(ucsc_alias) == old_values else old_values[0]
        )
    old_sequences = load_fasta_sequences(old_fasta, set(old_fasta_names.values()))
    new_sequences = load_fasta_sequences(new_fasta, set(accessions.values()))
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    run_alignment = alignment_runner or _run_alignment
    rows: list[dict[str, object]] = []
    for alias in requested_aliases:
        accession = accessions[alias]
        old_fasta_name = old_fasta_names[alias]
        comparison = compare_sequence_records(
            contig=alias,
            old_accession=old_accessions[alias],
            old_sequence=old_sequences[old_fasta_name],
            new_accession=accession,
            new_sequence=new_sequences[accession],
        )
        old_path = destination / f"BTA{alias}.old.fa"
        new_path = destination / f"BTA{alias}.new.fa"
        alignment_path = destination / f"BTA{alias}.old_to_new.paf"
        _write_single_fasta(old_path, old_accessions[alias], old_sequences[old_fasta_name])
        _write_single_fasta(new_path, accession, new_sequences[accession])
        command = [minimap2, "-x", "asm5", "-c", str(new_path), str(old_path)]
        run_alignment(command, alignment_path)
        if not alignment_path.is_file() or alignment_path.stat().st_size == 0:
            raise ReferenceRegistryError(f"alignment output is absent for BTA{alias}")
        rows.append(
            {
                "contig": alias,
                "old_accession": comparison.old_accession,
                "new_accession": comparison.new_accession,
                "old_fasta_name": old_fasta_name,
                "old_length": comparison.old_length,
                "new_length": comparison.new_length,
                "old_sha256": comparison.old_sha256,
                "new_sha256": comparison.new_sha256,
                "status": comparison.status,
                "alignment_status": "completed",
                "alignment_path": str(alignment_path),
            }
        )
    table_path = destination / "sequence_comparison.tsv"
    columns = [
        "contig",
        "old_accession",
        "new_accession",
        "old_fasta_name",
        "old_length",
        "new_length",
        "old_sha256",
        "new_sha256",
        "status",
        "alignment_status",
        "alignment_path",
    ]
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def summarize_ipd_alignments(
    manifest_path: Path,
    paf_path: Path,
    output_path: Path,
    *,
    minimum_identity: float = 0.95,
    minimum_coverage: float = 0.80,
    minimum_mapq: int = 20,
) -> dict[str, dict[str, object]]:
    """Summarize panel alignment evidence by the manifest's declared locus."""
    allele_loci: dict[str, str] = {}
    summaries: dict[str, dict[str, object]] = {}
    with Path(manifest_path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"allele_name", "locus"}
        if not required.issubset(reader.fieldnames or ()):
            raise ReferenceRegistryError("IPD manifest lacks allele_name or locus")
        for row in reader:
            allele = row["allele_name"]
            locus = row["locus"]
            if allele in allele_loci and allele_loci[allele] != locus:
                raise ReferenceRegistryError(f"IPD allele maps to multiple loci: {allele}")
            allele_loci[allele] = locus
            summaries.setdefault(
                locus,
                {
                    "locus": locus,
                    "probe_count": 0,
                    "alignment_count": 0,
                    "best_identity": 0.0,
                    "best_coverage": 0.0,
                    "best_mapq": 0,
                    "status": "no_alignment",
                    "_best_rank": None,
                },
            )["probe_count"] += 1
    with Path(paf_path).open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            fields = raw_line.rstrip("\n\r").split("\t")
            if len(fields) < 12:
                raise ReferenceRegistryError(f"invalid PAF row at line {line_number}")
            allele = fields[0]
            if allele not in allele_loci:
                raise ReferenceRegistryError(f"PAF query is absent from IPD manifest: {allele}")
            try:
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                matches = int(fields[9])
                block_length = int(fields[10])
                mapq = int(fields[11])
            except ValueError as exc:
                raise ReferenceRegistryError(f"invalid numeric PAF field at line {line_number}") from exc
            if (
                query_length <= 0
                or block_length <= 0
                or not 0 <= query_start <= query_end <= query_length
                or not 0 <= matches <= block_length
                or not 0 <= mapq <= 255
            ):
                raise ReferenceRegistryError(f"invalid PAF length at line {line_number}")
            summary = summaries[allele_loci[allele]]
            summary["alignment_count"] += 1
            identity = matches / block_length
            coverage = (query_end - query_start) / query_length
            qualified = (
                identity >= minimum_identity
                and coverage >= minimum_coverage
                and mapq >= minimum_mapq
            )
            rank = (qualified, coverage, identity, mapq)
            if summary["_best_rank"] is None or rank > summary["_best_rank"]:
                summary["_best_rank"] = rank
                summary["best_identity"] = identity
                summary["best_coverage"] = coverage
                summary["best_mapq"] = mapq
            if qualified:
                summary["status"] = "aligned"
            elif summary["status"] == "no_alignment":
                summary["status"] = "partial_alignment"
    columns = [
        "locus",
        "probe_count",
        "alignment_count",
        "best_identity",
        "best_coverage",
        "best_mapq",
        "status",
    ]
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for locus in sorted(summaries):
            writer.writerow({column: summaries[locus][column] for column in columns})
    for summary in summaries.values():
        summary.pop("_best_rank", None)
    return summaries


def align_ipd_panel(
    *,
    reference_fasta: Path,
    panel_path: Path,
    manifest_path: Path,
    output_dir: Path,
    minimap2: str,
    alignment_runner: AlignmentRunner | None = None,
) -> dict[str, dict[str, object]]:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    paf_path = destination / "ipd_mhc_panel.paf"
    command = [
        minimap2,
        "-x",
        "asm5",
        "-c",
        "--secondary=yes",
        "-N",
        "100",
        str(reference_fasta),
        str(panel_path),
    ]
    (alignment_runner or _run_alignment)(command, paf_path)
    return summarize_ipd_alignments(
        manifest_path,
        paf_path,
        destination / "ipd_mhc_locus_summary.tsv",
    )


def write_reference_manifest(
    output_path: Path,
    *,
    resources: Iterable[Mapping[str, object]],
    assembly_accession: str,
    annotation_release: str,
    retrieval_date: str,
) -> Path:
    """Write traceable checksums for each frozen external resource."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "assembly_accession",
        "annotation_release",
        "retrieval_date",
        "resource",
        "url",
        "size_bytes",
        "official_md5",
        "sha256",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for resource in resources:
            resource_path = Path(resource["path"])
            if not resource_path.is_file():
                raise ReferenceRegistryError(f"manifest resource is absent: {resource_path}")
            writer.writerow(
                {
                    "assembly_accession": assembly_accession,
                    "annotation_release": annotation_release,
                    "retrieval_date": retrieval_date,
                    "resource": str(resource["name"]),
                    "url": str(resource["url"]),
                    "size_bytes": resource_path.stat().st_size,
                    "official_md5": str(resource.get("official_md5", "")),
                    "sha256": sha256_file(resource_path),
                }
            )
    return path


def build_reference_registry(
    config_path: Path,
    output_dir: Path,
    *,
    command_runner: CommandRunner | None = None,
    retrieval_date: str | None = None,
    legacy_reference: Path | None = None,
    legacy_assembly_report: Path | None = None,
    legacy_checksum_manifest: Path | None = None,
    ipd_manifest: Path | None = None,
    ipd_panel: Path | None = None,
) -> dict[str, object]:
    """Build a verified reference and GFF-derived coordinate registry."""
    config_file = Path(config_path).resolve()
    project_root = _find_project_root(config_file)
    output = Path(output_dir).resolve()
    safe_root = (project_root / "results" / "codex_reanalysis_v1").resolve()
    if output != safe_root and safe_root not in output.parents:
        raise ReferenceRegistryError("output must be inside results/codex_reanalysis_v1")
    parsed = yaml.safe_load(config_file.read_text(encoding="utf-8"))
    if not isinstance(parsed, dict):
        raise ReferenceRegistryError("configuration must contain a mapping")
    reference = _mapping(parsed, "reference")
    tools = _mapping(parsed, "tools")
    accession = _string(reference, "accession")
    if accession != "GCF_002263795.3" or _string(reference, "assembly") != "ARS-UCD2.0":
        raise ReferenceRegistryError("reference must be ARS-UCD2.0 GCF_002263795.3")
    annotation_release = _string(reference, "annotation_release")
    base_url = _string(reference, "base_url").rstrip("/")
    parsed_url = urlparse(base_url)
    if parsed_url.scheme != "https" or parsed_url.hostname != "ftp.ncbi.nlm.nih.gov":
        raise ReferenceRegistryError("reference base_url must use the official NCBI HTTPS host")

    checksum_name = _safe_resource_name(_string(reference, "checksum_file"))
    resource_names = [
        _safe_resource_name(_string(reference, key))
        for key in (
            "fasta_file",
            "gff_file",
            "assembly_report_file",
            "annotation_report_file",
        )
    ]
    target_spec = _resolve_project_file(
        project_root,
        _string(reference, "target_spec_relative_path"),
    )
    if command_runner is None:
        samtools_version = _capture_tool_version(
            _string(tools, "samtools"),
            ["--version"],
            _string(tools, "samtools_version"),
        )
        minimap2_version = _capture_tool_version(
            _string(tools, "minimap2"),
            ["--version"],
            _string(tools, "minimap2_version"),
        )
    else:
        samtools_version = str(tools.get("samtools_version", "injected-test-runner"))
        minimap2_version = str(tools.get("minimap2_version", "injected-test-runner"))
    output.mkdir(parents=True, exist_ok=True)
    _write_json_atomic(output / "reference_status.json", {"status": "running"})
    downloads = output / "downloads"
    downloads.mkdir(exist_ok=True)
    checksum_path = downloads / checksum_name
    _ensure_download(f"{base_url}/{checksum_name}", checksum_path)
    expected = parse_official_checksums(checksum_path)
    for resource_name in resource_names:
        _ensure_download(f"{base_url}/{resource_name}", downloads / resource_name)
    verified = verify_official_checksums(downloads, expected, resource_names)

    fasta_gz, gff_gz, assembly_report, annotation_report = (
        downloads / name for name in resource_names
    )
    fasta_path = output / _without_gzip_suffix(resource_names[0])
    gff_path = output / _without_gzip_suffix(resource_names[1])
    _decompress_gzip(fasta_gz, fasta_path)
    _decompress_gzip(gff_gz, gff_path)

    runner = command_runner or _run_command
    commands = build_index_commands(
        _string(tools, "samtools"),
        _string(tools, "minimap2"),
        fasta_path,
    )
    index_paths = [Path(f"{fasta_path}.fai"), Path(f"{fasta_path}.mmi")]
    for command, index_path in zip(commands, index_paths, strict=True):
        partial_path = Path(command[-2] if command[1] == "faidx" else command[2])
        partial_path.unlink(missing_ok=True)
        runner(command)
        if not partial_path.is_file():
            raise ReferenceRegistryError(
                f"index command did not create {partial_path.name}"
            )
        partial_path.replace(index_path)

    aliases = parse_assembly_report_aliases(assembly_report)
    features = parse_gene_features(gff_path)
    specs = load_target_specs(target_spec)
    requested_symbols = {symbol for spec in specs for symbol in spec.gene_symbols}
    target_contigs = {feature.seqid for feature in features if feature.gene_symbol in requested_symbols}
    sequences = load_fasta_sequences(fasta_path, target_contigs)
    records = generate_target_records(
        features,
        specs,
        assembly_accession=accession,
        annotation_release=annotation_release,
        aliases=aliases,
        sequences=sequences,
    )
    probe_panel_status = "not_run"
    if (ipd_manifest is None) != (ipd_panel is None):
        raise ReferenceRegistryError("IPD manifest and panel must be supplied together")
    if ipd_manifest is not None and ipd_panel is not None:
        manifest_path = Path(ipd_manifest).resolve()
        panel_path = Path(ipd_panel).resolve()
        if not manifest_path.is_file() or not panel_path.is_file():
            raise ReferenceRegistryError("IPD manifest or panel is absent")
        locus_summaries = align_ipd_panel(
            reference_fasta=fasta_path,
            panel_path=panel_path,
            manifest_path=manifest_path,
            output_dir=output / "probe_alignments",
            minimap2=_string(tools, "minimap2"),
        )
        records = apply_probe_alignment_statuses(
            records,
            {locus: str(summary["status"]) for locus, summary in locus_summaries.items()},
        )
        probe_panel_status = "completed"
    coordinate_outputs = write_coordinate_outputs(records, sequences=sequences, output_dir=output)

    comparison_status = "not_run"
    comparisons: list[dict[str, object]] = []
    if legacy_reference is not None:
        if legacy_assembly_report is None or legacy_checksum_manifest is None:
            raise ReferenceRegistryError(
                "legacy comparison requires an assembly report and checksum manifest"
            )
        legacy_path = Path(legacy_reference).resolve()
        legacy_report_path = Path(legacy_assembly_report).resolve()
        legacy_checksums_path = Path(legacy_checksum_manifest).resolve()
        if not legacy_path.is_file():
            raise ReferenceRegistryError(f"legacy comparison reference is absent: {legacy_path}")
        if not legacy_report_path.is_file() or not legacy_checksums_path.is_file():
            raise ReferenceRegistryError("legacy assembly report or checksum manifest is absent")
        legacy_expected = parse_official_checksums(legacy_checksums_path)
        legacy_verified = verify_official_checksums(
            legacy_report_path.parent,
            legacy_expected,
            [legacy_report_path.name],
        )
        old_aliases = parse_assembly_report_aliases(legacy_report_path)
        configured_contigs = reference.get("comparison_contigs")
        if not isinstance(configured_contigs, list) or not configured_contigs:
            raise ReferenceRegistryError("reference.comparison_contigs must be a non-empty list")
        comparisons = compare_fasta_contigs(
            old_fasta=legacy_path,
            new_fasta=fasta_path,
            old_aliases=old_aliases,
            aliases=aliases,
            contigs=[str(value) for value in configured_contigs],
            output_dir=output / "assembly_comparison",
            minimap2=_string(tools, "minimap2"),
        )
        comparison_status = (
            "identical"
            if all(row["status"] == "identical" for row in comparisons)
            else "decision_required"
        )

    retrieved = retrieval_date or date.today().isoformat()
    manifest_resources = [
        {
            "name": checksum_name,
            "path": checksum_path,
            "url": f"{base_url}/{checksum_name}",
            "official_md5": "",
        },
        *[
        {
            "name": name,
            "path": downloads / name,
            "url": f"{base_url}/{name}",
            "official_md5": verified[name],
        }
        for name in resource_names
        ],
    ]
    manifest_resources.extend(
        {
            "name": path.name,
            "path": path,
            "url": f"generated-from:{source.name}",
            "official_md5": "",
        }
        for path, source in (
            (fasta_path, fasta_gz),
            (gff_path, gff_gz),
            (index_paths[0], fasta_path),
            (index_paths[1], fasta_path),
        )
    )
    manifest_resources.extend(
        {
            "name": path.name,
            "path": path,
            "url": "generated-coordinate-registry",
            "official_md5": "",
        }
        for path in coordinate_outputs.values()
    )
    if legacy_reference is not None:
        manifest_resources.extend(
            [
                {
                    "name": Path(legacy_checksum_manifest).name,
                    "path": Path(legacy_checksum_manifest).resolve(),
                    "url": "local-official-checksum-manifest",
                    "official_md5": "",
                },
                {
                    "name": Path(legacy_assembly_report).name,
                    "path": Path(legacy_assembly_report).resolve(),
                    "url": "official-ncbi-legacy-assembly-report",
                    "official_md5": legacy_verified[Path(legacy_assembly_report).name],
                },
            ]
        )
    manifest = write_reference_manifest(
        output / "reference_manifest.tsv",
        resources=manifest_resources,
        assembly_accession=accession,
        annotation_release=annotation_release,
        retrieval_date=retrieved,
    )
    overall_status = "incomplete"
    if comparison_status == "decision_required":
        overall_status = "decision_required"
    elif comparison_status == "identical" and probe_panel_status == "completed":
        overall_status = "passed"
    status: dict[str, object] = {
        "status": overall_status,
        "assembly_accession": accession,
        "annotation_release": annotation_release,
        "reference_manifest": str(manifest),
        "reference_manifest_sha256": sha256_file(manifest),
        "target_count": len(records),
        "reference_feature_count": sum(row.target_type == "reference_feature" for row in records),
        "probe_only_count": sum(row.target_type == "haplotype_probe" for row in records),
        "comparison_status": comparison_status,
        "comparison_count": len(comparisons),
        "probe_panel_status": probe_panel_status,
    }
    provenance_records = {
            "versions": {
                "python": __import__("platform").python_version(),
                "samtools": samtools_version,
                "minimap2": minimap2_version,
                "annotation_pipeline": str(reference.get("annotation_pipeline_version", "unknown")),
            },
            "input_hashes": {
                "target_spec": sha256_file(target_spec),
                **(
                    {
                        "ipd_manifest": sha256_file(Path(ipd_manifest).resolve()),
                        "ipd_panel": sha256_file(Path(ipd_panel).resolve()),
                    }
                    if ipd_manifest is not None and ipd_panel is not None
                    else {}
                ),
            },
            "reference_hashes": {
                **{name: sha256_file(downloads / name) for name in resource_names},
                **(
                    {
                        "legacy_reference": sha256_file(Path(legacy_reference).resolve()),
                        "legacy_assembly_report": sha256_file(
                            Path(legacy_assembly_report).resolve()
                        ),
                        "legacy_checksum_manifest": sha256_file(
                            Path(legacy_checksum_manifest).resolve()
                        ),
                    }
                    if legacy_reference is not None
                    else {}
                ),
            },
            "config_hash": sha256_file(config_file),
        }
    _publish_reference_completion(output, status, provenance_records)
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--legacy-reference", type=Path)
    parser.add_argument("--legacy-assembly-report", type=Path)
    parser.add_argument("--legacy-checksum-manifest", type=Path)
    parser.add_argument("--ipd-manifest", type=Path)
    parser.add_argument("--ipd-panel", type=Path)
    args = parser.parse_args(argv)
    build_reference_registry(
        args.config,
        args.output,
        legacy_reference=args.legacy_reference,
        legacy_assembly_report=args.legacy_assembly_report,
        legacy_checksum_manifest=args.legacy_checksum_manifest,
        ipd_manifest=args.ipd_manifest,
        ipd_panel=args.ipd_panel,
    )
    return 0


def _md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    destination = Path(path)
    temporary = destination.with_suffix(f"{destination.suffix}.partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(destination)


def _publish_reference_completion(
    output: Path,
    status: Mapping[str, object],
    provenance_records: dict[str, object],
) -> None:
    """Publish passed status only after provenance is safely available."""
    destination = Path(output)
    write_provenance(destination, provenance_records)
    _write_json_atomic(destination / "reference_status.json", status)


def _find_project_root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise ReferenceRegistryError("configuration must be inside a Git repository")


def _mapping(config: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = config.get(key)
    if not isinstance(value, dict):
        raise ReferenceRegistryError(f"configuration section {key!r} must be a mapping")
    return value


def _string(config: Mapping[str, object], key: str) -> str:
    value = config.get(key)
    if not isinstance(value, str) or not value:
        raise ReferenceRegistryError(f"configuration value {key!r} must be a string")
    return value


def _safe_resource_name(name: str) -> str:
    path = PurePosixPath(name)
    if (
        path.is_absolute()
        or len(path.parts) != 1
        or path.name != name
        or name in {".", ".."}
        or any(character in name for character in ("/", "\\", ":", "?", "#"))
        or any(ord(character) < 32 for character in name)
    ):
        raise ReferenceRegistryError(f"unsafe resource name: {name}")
    return name


def _resolve_project_file(project_root: Path, relative: str) -> Path:
    candidate = Path(relative)
    if candidate.is_absolute():
        raise ReferenceRegistryError("tracked input path must be project-relative")
    resolved = (project_root / candidate).resolve()
    if resolved != project_root and project_root not in resolved.parents:
        raise ReferenceRegistryError("tracked input path escapes the project")
    if not resolved.is_file():
        raise ReferenceRegistryError(f"tracked input is absent: {relative}")
    return resolved


def _ensure_download(url: str, destination: Path) -> None:
    if destination.is_file():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with urlopen(url) as response, NamedTemporaryFile(
            mode="wb", dir=destination.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            shutil.copyfileobj(response, handle, length=1024 * 1024)
        temporary.replace(destination)
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def _without_gzip_suffix(name: str) -> str:
    if not name.endswith(".gz"):
        raise ReferenceRegistryError(f"expected a gzip resource: {name}")
    return name[:-3]


def _decompress_gzip(source: Path, destination: Path) -> None:
    temporary = destination.with_name(f".{destination.name}.partial")
    try:
        with gzip.open(source, "rb") as input_handle, temporary.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _run_command(command: list[str]) -> None:
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise ReferenceRegistryError(
            f"reference command failed with exit code {result.returncode}: {command}"
        )


def _capture_tool_version(
    command: str,
    arguments: list[str],
    expected: str,
) -> str:
    result = subprocess.run(
        [command, *arguments],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise ReferenceRegistryError(f"cannot query {command} version")
    combined = result.stdout or result.stderr
    return validate_tool_version(combined, expected, command)


def _run_alignment(command: list[str], output_path: Path) -> None:
    temporary = output_path.with_suffix(f"{output_path.suffix}.partial")
    with temporary.open("w", encoding="ascii", newline="") as handle:
        result = subprocess.run(command, stdout=handle, check=False, text=True)
    if result.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise ReferenceRegistryError(
            f"reference alignment failed with exit code {result.returncode}: {command}"
        )
    temporary.replace(output_path)


def _write_single_fasta(path: Path, name: str, sequence: str) -> None:
    with path.open("w", encoding="ascii", newline="") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
