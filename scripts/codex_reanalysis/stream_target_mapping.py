"""Stream minimap2 alignments into a validated target-only BAM."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import yaml

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file


#: minimap2 keeps 5 secondary alignments by default, which truncates the
#: multi-hit evidence BoLA class I paralogue analysis relies on.
DEFAULT_SECONDARY_MAX = 20


class MappingError(RuntimeError):
    """Raised when mapping cannot produce a validated target alignment."""


def build_mapping_argv(
    fastq: str | Path,
    reference: str | Path,
    *,
    threads: int,
    sample_id: str,
    deliverable_id: str,
    minimap2: str = "minimap2",
    secondary_max: int = DEFAULT_SECONDARY_MAX,
) -> list[str]:
    """Build the mapper argv while retaining secondary/supplementary records."""
    if threads < 1:
        raise MappingError("threads must be positive")
    if secondary_max < 1:
        raise MappingError("secondary alignment cap must be positive")
    _validate_identifier(sample_id, "sample ID")
    _validate_identifier(deliverable_id, "deliverable ID")
    read_group = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ONT\\tPU:{deliverable_id}"
    return [
        minimap2,
        "-ax",
        "map-ont",
        "-L",
        "--secondary=yes",
        "-N",
        str(secondary_max),
        "-t",
        str(threads),
        "-R",
        read_group,
        str(reference),
        str(fastq),
    ]


def build_mapping_command(
    fastq: str | Path,
    reference: str | Path,
    *,
    threads: int,
    sample_id: str = "sample",
    deliverable_id: str = "deliverable",
    minimap2: str = "minimap2",
    secondary_max: int = DEFAULT_SECONDARY_MAX,
) -> str:
    """Return a shell-escaped display command; execution uses the argv builder."""
    return shlex.join(
        build_mapping_argv(
            fastq,
            reference,
            threads=threads,
            sample_id=sample_id,
            deliverable_id=deliverable_id,
            minimap2=minimap2,
            secondary_max=secondary_max,
        )
    )


def build_redacted_mapping_command(
    fastq: str | Path,
    reference: str | Path,
    *,
    threads: int,
    sample_id: str,
    deliverable_id: str,
    minimap2: str = "minimap2",
    secondary_max: int = DEFAULT_SECONDARY_MAX,
) -> str:
    """Build a reproducible display command without exposing delivery paths."""
    del fastq
    return build_mapping_command(
        f"<audited-fastq:{sample_id}>",
        Path(reference).name,
        threads=threads,
        sample_id=sample_id,
        deliverable_id=deliverable_id,
        minimap2=minimap2,
        secondary_max=secondary_max,
    )


def build_target_filter_argv(
    targets_bed: str | Path,
    output_bam: str | Path,
    *,
    temp_prefix: str | Path,
    threads: int,
    samtools: str = "samtools",
) -> tuple[list[str], list[str]]:
    """Build the target filter and filtered-only coordinate-sort argv vectors."""
    if threads < 1:
        raise MappingError("threads must be positive")
    view = [samtools, "view", "-h", "-b", "-L", str(targets_bed), "-"]
    sort = [
        samtools,
        "sort",
        "-@",
        str(threads),
        "-T",
        str(temp_prefix),
        "-o",
        str(output_bam),
        "-",
    ]
    return view, sort


def build_target_filter_command(
    targets_bed: str | Path,
    output_bam: str | Path,
    *,
    temp_prefix: str | Path = "target-sort",
    threads: int = 1,
    samtools: str = "samtools",
) -> str:
    """Return the logged filter pipeline without executing a shell string."""
    view, sort = build_target_filter_argv(
        targets_bed,
        output_bam,
        temp_prefix=temp_prefix,
        threads=threads,
        samtools=samtools,
    )
    return f"{shlex.join(view)} | {shlex.join(sort)}"


def build_redacted_filter_commands(
    *, threads: int, samtools: str = "samtools"
) -> tuple[str, str]:
    """Describe target filtering without recording local filesystem paths."""
    view, sort = build_target_filter_argv(
        "<validated-target-bed>",
        "<staged-target-bam>",
        temp_prefix="<sample-temp-prefix>",
        threads=threads,
        samtools=samtools,
    )
    return shlex.join(view), shlex.join(sort)


def validate_stage_gates(
    audit_status_path: Path,
    input_manifest_path: Path,
    sample_registry_path: Path,
    reference_status_path: Path,
    reference_manifest_path: Path,
) -> None:
    """Refuse mapping unless input and reference stages are explicitly complete."""
    audit = _read_json_object(audit_status_path, "input audit status")
    if audit.get("status") != "pass" or audit.get("analytical_sample_count") != 13:
        raise MappingError("input audit is not passing for exactly 13 analytical samples")
    if (
        audit.get("input_manifest_sha256") != sha256_file(input_manifest_path)
        or audit.get("sample_registry_sha256") != sha256_file(sample_registry_path)
    ):
        raise MappingError("audit status is not bound to the active manifest/registry")
    reference = _read_json_object(reference_status_path, "reference status")
    if reference.get("status") != "passed":
        raise MappingError("reference registry is not passing")
    manifest = Path(reference_manifest_path)
    if not manifest.is_file() or manifest.stat().st_size == 0:
        raise MappingError("reference manifest is absent or empty")
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if not {"resource", "sha256"}.issubset(fields) or next(reader, None) is None:
            raise MappingError("reference manifest is malformed or empty")
    if reference.get("reference_manifest_sha256") != sha256_file(manifest):
        raise MappingError("reference status is not bound to the active manifest")


def resolve_audited_sample(
    sample_id: str,
    sample_registry_path: Path,
    input_manifest_path: Path,
    deliver_root: Path,
    normalized_manifest_path: Path | None = None,
    input_seal: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Resolve one analytical FASTQ only through matching audited records.

    When a normalisation register is supplied, a sample it marks as filtered is
    mapped from the registered derivative instead of the delivered file. The
    derivative is accepted only if it is bound to the same delivered SHA-256 the
    audit recorded, so a derivative built from some other input cannot enter.

    When an input seal is supplied, the delivered SHA-256 it records was
    obtained by hashing the delivery on this host, so the record can state the
    hash as verified here rather than carried over from the audit. A seal that
    does not cover the sample, or disagrees with the audit, stops the mapping
    instead of quietly falling back to the weaker claim.
    """
    _validate_identifier(sample_id, "sample ID")
    root = Path(deliver_root).resolve()
    if not root.is_dir():
        raise MappingError("delivery root is absent")
    registry_rows = _read_tsv(sample_registry_path, "sample registry")
    registry_matches = [
        row
        for row in registry_rows
        if row.get("sample_id") == sample_id and row.get("analysis_status") == "analytical"
    ]
    if len(registry_matches) != 1:
        raise MappingError(f"sample registry does not contain one analytical row for {sample_id}")
    registry = registry_matches[0]
    manifest_rows = _read_tsv(input_manifest_path, "input manifest")
    manifest_matches = [row for row in manifest_rows if row.get("sample_id") == sample_id]
    if len(manifest_matches) != 1:
        raise MappingError(f"input manifest does not contain one row for {sample_id}")
    manifest = manifest_matches[0]
    for field in ("deliverable_id", "input_relative_path"):
        if not registry.get(field) or registry.get(field) != manifest.get(field):
            raise MappingError(f"audit tables disagree on {field} for {sample_id}")
    relative = Path(str(registry["input_relative_path"]))
    if relative.is_absolute():
        raise MappingError("audited FASTQ path must be delivery-root-relative")
    fastq = (root / relative).resolve()
    if fastq != root and root not in fastq.parents:
        raise MappingError("audited FASTQ path escapes the delivery root")
    if not fastq.is_file():
        raise MappingError(f"audited FASTQ is absent: {fastq}")
    digest = str(manifest.get("sha256", "")).strip()
    if not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
        raise MappingError("input manifest lacks a valid audited FASTQ SHA-256")
    try:
        recorded_size = int(str(manifest.get("byte_size", "")))
    except ValueError as exc:
        raise MappingError("input manifest has an invalid FASTQ byte size") from exc
    if recorded_size != fastq.stat().st_size:
        raise MappingError("audited FASTQ byte size does not match the input manifest")
    resolved: dict[str, object] = {
        "deliverable_id": str(registry["deliverable_id"]),
        "fastq": fastq,
        "input_sha256": digest,
        "source_sha256": digest,
        "input_normalization": "delivered",
        # The audit hashed this file; mapping re-checks only its size, so the
        # record must not imply the hash was recomputed here.
        "input_hash_verified_at_mapping": False,
        "input_hash_verification": "not_verified",
    }
    if input_seal is not None:
        _apply_input_seal(resolved, sample_id, input_seal, digest, recorded_size)
    if normalized_manifest_path is not None:
        resolved = _apply_normalized_input(resolved, sample_id, Path(normalized_manifest_path))
    return resolved


def _apply_input_seal(
    resolved: dict[str, object],
    sample_id: str,
    seal: Mapping[str, object],
    digest: str,
    recorded_size: int,
) -> None:
    """Upgrade the delivered hash from carried over to verified on this host."""
    if seal.get("status") != "verified":
        raise MappingError(f"input seal is not verified: {seal.get('status')}")
    samples = seal.get("samples")
    if not isinstance(samples, Mapping) or sample_id not in samples:
        raise MappingError(f"input seal does not cover {sample_id}")
    record = samples[sample_id]
    if not isinstance(record, Mapping):
        raise MappingError(f"input seal record is malformed for {sample_id}")
    if str(record.get("sha256", "")).lower() != digest.lower():
        raise MappingError(f"input seal disagrees with the audited SHA-256 for {sample_id}")
    if int(record.get("byte_size", -1)) != recorded_size:
        raise MappingError(f"input seal disagrees with the audited byte size for {sample_id}")
    resolved["input_hash_verified_at_mapping"] = True
    resolved["input_hash_verification"] = "input_seal"


def _apply_normalized_input(
    resolved: dict[str, object], sample_id: str, manifest_path: Path
) -> dict[str, object]:
    """Redirect one sample to its registered length-filtered derivative."""
    rows = [
        row
        for row in _read_tsv(manifest_path, "normalization manifest")
        if row.get("sample_id") == sample_id
    ]
    if len(rows) != 1:
        raise MappingError(f"normalization manifest does not contain one row for {sample_id}")
    row = rows[0]
    if str(row.get("source_sha256", "")).strip().lower() != str(resolved["source_sha256"]).lower():
        raise MappingError(
            f"normalized derivative for {sample_id} was built from a different delivered source"
        )
    action = str(row.get("action", "")).strip()
    if action == "not_required":
        return resolved
    if action != "length_filter":
        raise MappingError(f"unknown normalization action for {sample_id}: {action!r}")

    relative = Path(str(row.get("normalized_relative_path", "")).strip())
    if not relative.name or relative.is_absolute() or len(relative.parts) != 1:
        raise MappingError(f"normalized derivative path is not a plain file name for {sample_id}")
    derived = (manifest_path.parent / relative).resolve()
    if derived.parent != manifest_path.parent.resolve() or not derived.is_file():
        raise MappingError(f"normalized derivative is absent for {sample_id}: {derived}")
    expected = str(row.get("normalized_sha256", "")).strip()
    if not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        raise MappingError(f"normalization manifest lacks a valid derivative SHA-256 for {sample_id}")
    try:
        recorded_size = int(str(row.get("normalized_byte_size", "")))
    except ValueError as exc:
        raise MappingError(f"normalization manifest has an invalid derivative size for {sample_id}") from exc
    if recorded_size != derived.stat().st_size or sha256_file(derived).lower() != expected.lower():
        raise MappingError(f"normalized derivative does not match its registered checksum for {sample_id}")
    resolved.update(
        {
            "fastq": derived,
            "input_sha256": expected.lower(),
            "input_normalization": "length_filtered",
            "input_hash_verified_at_mapping": True,
            # The derivative's own hash was recomputed here; whether the
            # delivered source was re-read is a separate question the seal
            # answers, and the two must not be conflated.
            "input_hash_verification": (
                "input_seal+recomputed_derivative"
                if resolved.get("input_hash_verification") == "input_seal"
                else "recomputed_derivative"
            ),
        }
    )
    return resolved


def reference_resource_sha256(manifest_path: Path, resource_name: str) -> str:
    """Reuse the verified resource digest instead of re-reading a multi-GB FASTA."""
    rows = _read_tsv(manifest_path, "reference manifest")
    matches = [row for row in rows if row.get("resource") == resource_name]
    if len(matches) != 1 or not str(matches[0].get("sha256", "")).strip():
        raise MappingError(f"reference manifest lacks one SHA-256 for {resource_name}")
    return str(matches[0]["sha256"]).strip()


def create_reference_seal(
    reference: Path,
    targets_bed: Path,
    manifest_path: Path,
    config_path: Path,
    output_path: Path,
) -> dict[str, object]:
    """Hash live mapping resources once and seal them to manifest/config inputs."""
    resources = {}
    for path in (Path(reference), Path(targets_bed)):
        expected = reference_resource_sha256(manifest_path, path.name)
        actual = sha256_file(path)
        if actual.lower() != expected.lower():
            raise MappingError(f"reference manifest SHA-256 mismatch for {path.name}")
        stat = path.stat()
        resources[path.name] = {
            "sha256": actual,
            "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    seal: dict[str, object] = {
        "status": "passed",
        "manifest_sha256": sha256_file(manifest_path),
        "config_sha256": sha256_file(config_path),
        "resources": resources,
    }
    _write_json_atomic(Path(output_path), seal)
    return seal


def validate_reference_seal(
    seal_path: Path,
    reference: Path,
    targets_bed: Path,
    manifest_path: Path,
    config_path: Path,
) -> dict[str, object]:
    """Bind each map job to the exact files hashed by the shared seal job."""
    seal = _read_json_object(seal_path, "reference seal")
    if (
        seal.get("status") != "passed"
        or seal.get("manifest_sha256") != sha256_file(manifest_path)
        or seal.get("config_sha256") != sha256_file(config_path)
    ):
        raise MappingError("reference seal does not match the active manifest/config")
    resources = seal.get("resources")
    if not isinstance(resources, dict):
        raise MappingError("reference seal lacks resource records")
    for path in (Path(reference), Path(targets_bed)):
        record = resources.get(path.name)
        if not isinstance(record, dict):
            raise MappingError(f"reference seal lacks {path.name}")
        stat = path.stat()
        if record.get("size_bytes") != stat.st_size or record.get("mtime_ns") != stat.st_mtime_ns:
            raise MappingError(f"reference seal no longer matches {path.name}")
        # size/mtime match is a cheap proxy, not proof: a same-size overwrite
        # can land inside one filesystem timestamp tick (exFAT's is coarse,
        # and even local filesystems are not immune), so the recorded
        # SHA-256 -- which the size/mtime pair exists only to let us usually
        # skip recomputing -- is re-checked on every call rather than trusted
        # implicitly.
        if record.get("sha256") != sha256_file(path):
            raise MappingError(f"reference seal no longer matches {path.name}")
    return seal


def validate_mapping_config(
    config_path: Path,
    threads: int,
    minimap2_version: str,
    samtools_version: str,
    secondary_max: int | None = None,
) -> dict[str, object]:
    """Enforce the approved mapping semantics and pinned tool versions."""
    try:
        payload = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise MappingError("mapping config is unreadable") from exc
    if not isinstance(payload, dict):
        raise MappingError("mapping config must be a mapping")
    mapping = payload.get("mapping")
    tools = payload.get("tools")
    if not isinstance(mapping, dict) or not isinstance(tools, dict):
        raise MappingError("mapping config lacks mapping/tools sections")
    if mapping.get("preset") != "map-ont":
        raise MappingError("mapping preset must be map-ont")
    if mapping.get("retain_secondary") is not True or mapping.get("retain_supplementary") is not True:
        raise MappingError("secondary and supplementary retention must be enabled")
    if threads < 1 or threads != mapping.get("threads"):
        raise MappingError("runtime threads must equal the positive configured thread count")
    configured_secondary = mapping.get("secondary_max")
    if not isinstance(configured_secondary, int) or configured_secondary < 1:
        raise MappingError("mapping config must set a positive secondary_max")
    if secondary_max is not None and secondary_max != configured_secondary:
        raise MappingError(
            f"runtime secondary cap {secondary_max} does not match the configured "
            f"secondary_max {configured_secondary}"
        )
    expected_minimap2 = str(tools.get("minimap2_version", ""))
    expected_samtools = str(tools.get("samtools_version", ""))
    if minimap2_version.strip() != expected_minimap2:
        raise MappingError("minimap2 version does not match the configured pin")
    samtools_tokens = samtools_version.strip().split()
    if len(samtools_tokens) < 2 or samtools_tokens[:2] != ["samtools", expected_samtools]:
        raise MappingError("samtools version does not match the configured pin")
    return payload


def summarize_sam_alignments(lines: Iterable[str]) -> dict[str, int]:
    """Summarize MAPQ and alignment class from headerless SAM rows."""
    result = {
        "alignment_count": 0,
        "primary_count": 0,
        "secondary_count": 0,
        "supplementary_count": 0,
        "mapq_0": 0,
        "mapq_1_9": 0,
        "mapq_10_19": 0,
        "mapq_20_59": 0,
        "mapq_60_plus": 0,
    }
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\n\r")
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise MappingError(f"invalid SAM row at line {line_number}")
        try:
            flag = int(fields[1])
            mapq = int(fields[4])
        except ValueError as exc:
            raise MappingError(f"invalid SAM flag or MAPQ at line {line_number}") from exc
        if flag < 0 or not 0 <= mapq <= 255:
            raise MappingError(f"out-of-range SAM flag or MAPQ at line {line_number}")
        result["alignment_count"] += 1
        if flag & 0x100:
            result["secondary_count"] += 1
        elif flag & 0x800:
            result["supplementary_count"] += 1
        else:
            result["primary_count"] += 1
        if mapq == 0:
            result["mapq_0"] += 1
        elif mapq < 10:
            result["mapq_1_9"] += 1
        elif mapq < 20:
            result["mapq_10_19"] += 1
        elif mapq < 60:
            result["mapq_20_59"] += 1
        else:
            result["mapq_60_plus"] += 1
    return result


def validate_mapping_summary(summary: Mapping[str, int]) -> None:
    """Reject an empty target result instead of publishing a biological zero."""
    if int(summary.get("alignment_count", 0)) < 1:
        raise MappingError("zero target alignments; mapping output is not usable")


def parse_flagstat_total(path: Path) -> int:
    """Read the exact total alignment-record count from samtools flagstat."""
    source = Path(path)
    if not source.is_file():
        raise MappingError("flagstat output is absent")
    lines = source.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise MappingError("flagstat output is empty")
    match = re.fullmatch(
        r"(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)",
        lines[0].strip(),
    )
    if match is None:
        raise MappingError("flagstat total line is malformed")
    return int(match.group(1)) + int(match.group(2))


def publish_staged_artifacts(artifacts: Mapping[Path, Path]) -> None:
    """Publish a complete validated artifact set; never publish an incomplete set."""
    missing = [str(source) for source in artifacts if not Path(source).is_file()]
    if missing:
        raise MappingError(f"staged mapping artifact is absent: {', '.join(missing)}")
    for destination in artifacts.values():
        Path(destination).parent.mkdir(parents=True, exist_ok=True)
    for source, destination in artifacts.items():
        Path(source).replace(destination)


def invalidate_status_before_cleanup(
    status_path: Path,
    sample_id: str,
    stale_paths: Iterable[Path],
    staging_dir: Path | None = None,
) -> None:
    """Invalidate any prior success marker before touching published artifacts.

    Published artifacts are removed by name because the destination holds other
    samples. The staging directory belongs to this sample alone and is cleared
    wholesale: enumerating files leaves behind anything not on the list, and a
    leftover from an earlier run can satisfy a later existence check.
    """
    _write_json_atomic(status_path, {"status": "running", "sample_id": sample_id})
    for stale in stale_paths:
        Path(stale).unlink(missing_ok=True)
    if staging_dir is not None:
        staging = Path(staging_dir)
        shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True, exist_ok=True)


def run_stream_mapping(
    *,
    sample_id: str,
    deliverable_id: str,
    fastq: Path,
    input_manifest: Path,
    sample_registry: Path,
    reference: Path,
    targets_bed: Path,
    audit_status: Path,
    reference_status: Path,
    reference_manifest: Path,
    reference_seal: Path,
    config_path: Path,
    output_dir: Path,
    threads: int,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    input_sha256: str,
    source_sha256: str | None = None,
    input_normalization: str = "delivered",
    input_hash_verified: bool = False,
    input_hash_verification: str = "not_verified",
    secondary_max: int = DEFAULT_SECONDARY_MAX,
) -> dict[str, object]:
    """Map one FASTQ and persist only validated target/control alignments."""
    _validate_identifier(sample_id, "sample ID")
    _validate_identifier(deliverable_id, "deliverable ID")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", input_sha256):
        raise MappingError("mapping requires an audited FASTQ SHA-256")
    validate_stage_gates(
        audit_status,
        input_manifest,
        sample_registry,
        reference_status,
        reference_manifest,
    )
    for path, label in (
        (Path(fastq), "FASTQ"),
        (Path(reference), "reference FASTA"),
        (Path(targets_bed), "target BED"),
        (Path(config_path), "workflow config"),
    ):
        if not path.is_file():
            raise MappingError(f"{label} is absent: {path}")
    reference_seal_payload = validate_reference_seal(
        reference_seal, reference, targets_bed, reference_manifest, config_path
    )
    minimap2_version = _capture_version([minimap2, "--version"])
    samtools_version = _capture_version([samtools, "--version"])
    validate_mapping_config(
        config_path, threads, minimap2_version, samtools_version, secondary_max=secondary_max
    )
    destination = Path(output_dir).resolve()
    _validate_mapping_output(destination)
    destination.mkdir(parents=True, exist_ok=True)
    temporary_dir = destination / "tmp" / sample_id
    temporary_dir.mkdir(parents=True, exist_ok=True)
    status_path = destination / f"{sample_id}.mapping_status.json"

    staging_dir = temporary_dir / "publish"
    staging_dir.mkdir(parents=True, exist_ok=True)
    bam = destination / f"{sample_id}.targets.bam"
    bam_partial = staging_dir / f"{sample_id}.targets.bam"
    bai = destination / f"{sample_id}.targets.bam.bai"
    bai_partial = staging_dir / f"{sample_id}.targets.bam.bai"
    flagstat = destination / f"{sample_id}.flagstat.txt"
    flagstat_partial = staging_dir / flagstat.name
    mapq_table = destination / f"{sample_id}.mapq_summary.tsv"
    mapq_partial = staging_dir / mapq_table.name
    command_log = destination / f"{sample_id}.commands.json"
    command_partial = staging_dir / command_log.name
    mapper_log = destination / f"{sample_id}.minimap2.log"
    samtools_log = destination / f"{sample_id}.samtools.log"
    provenance = destination / "provenance" / sample_id / "provenance.json"
    invalidate_status_before_cleanup(
        status_path,
        sample_id,
        (
            bam,
            bai,
            flagstat,
            mapq_table,
            command_log,
            provenance,
        ),
        staging_dir=staging_dir,
    )

    mapping_argv = build_mapping_argv(
        fastq,
        reference,
        threads=threads,
        sample_id=sample_id,
        deliverable_id=deliverable_id,
        minimap2=minimap2,
        secondary_max=secondary_max,
    )
    view_argv, sort_argv = build_target_filter_argv(
        targets_bed,
        bam_partial,
        temp_prefix=temporary_dir / "sort",
        threads=threads,
        samtools=samtools,
    )
    redacted_view, redacted_sort = build_redacted_filter_commands(
        threads=threads, samtools=samtools
    )
    _write_json_atomic(
        command_partial,
        {
            "mapping": build_redacted_mapping_command(
                fastq,
                reference,
                threads=threads,
                sample_id=sample_id,
                deliverable_id=deliverable_id,
                minimap2=minimap2,
                secondary_max=secondary_max,
            ),
            "filter": redacted_view,
            "sort": redacted_sort,
        },
    )
    try:
        _run_pipeline(mapping_argv, view_argv, sort_argv, mapper_log, samtools_log)
        if not bam_partial.is_file() or bam_partial.stat().st_size == 0:
            raise MappingError("target BAM is absent or empty")
        _run_checked([samtools, "quickcheck", "-v", str(bam_partial)], samtools_log)
        _run_checked(
            [samtools, "index", "-@", str(threads), str(bam_partial), str(bai_partial)],
            samtools_log,
        )
        if not bai_partial.is_file() or bai_partial.stat().st_size == 0:
            raise MappingError("target BAM index is absent or empty")
        idxstats_partial = staging_dir / f"{sample_id}.idxstats.txt"
        _capture_command([samtools, "idxstats", str(bam_partial)], idxstats_partial)
        if idxstats_partial.stat().st_size == 0:
            raise MappingError("target BAM index cannot be read by samtools idxstats")
        _capture_command(
            [samtools, "flagstat", "-@", str(threads), str(bam_partial)], flagstat_partial
        )
        summary = _summarize_bam(samtools, bam_partial, samtools_log)
        validate_mapping_summary(summary)
        flagstat_count = parse_flagstat_total(flagstat_partial)
        if flagstat_count != summary["alignment_count"]:
            raise MappingError(
                "flagstat total does not match the independent alignment count: "
                f"{flagstat_count} != {summary['alignment_count']}"
            )
        _write_summary(mapq_partial, sample_id, summary)
        output_hashes = {
            "target_bam": sha256_file(bam_partial),
            "target_bai": sha256_file(bai_partial),
            "flagstat": sha256_file(flagstat_partial),
            "mapq_summary": sha256_file(mapq_partial),
            "commands": sha256_file(command_partial),
        }
        provenance_partial_dir = staging_dir / "provenance" / sample_id
        provenance_partial = write_provenance(
            provenance_partial_dir,
            {
                "command": f"stream_target_mapping.py --sample-id {sample_id} [paths recorded separately]",
                "versions": {
                    "python": sys.version.split()[0],
                    "minimap2": minimap2_version,
                    "samtools": samtools_version,
                },
                "input_hashes": {
                    "fastq": input_sha256,
                    "delivered_source": source_sha256 or input_sha256,
                    "targets_bed": sha256_file(targets_bed),
                },
                "input_normalization": input_normalization,
                "input_hash_verified_at_mapping": input_hash_verified,
                "input_hash_verification": input_hash_verification,
                "reference_hashes": {
                    "reference_manifest": sha256_file(reference_manifest),
                    "reference_seal": sha256_file(reference_seal),
                    "reference_fasta": reference_seal_payload["resources"][Path(reference).name]["sha256"],
                },
                "config_hash": sha256_file(config_path),
                "output_hashes": output_hashes,
                "alignment_retention_scope": "all alignment classes whose records overlap the target/control BED",
            },
        )
        provenance_hash = sha256_file(provenance_partial)
        publish_staged_artifacts(
            {
                bam_partial: bam,
                bai_partial: bai,
                flagstat_partial: flagstat,
                mapq_partial: mapq_table,
                command_partial: command_log,
                provenance_partial: provenance,
            }
        )
        status: dict[str, object] = {
            "status": "passed",
            "sample_id": sample_id,
            "deliverable_id": deliverable_id,
            "input_normalization": input_normalization,
            "flagstat_alignment_count": flagstat_count,
            "artifact_hashes": {**output_hashes, "provenance": provenance_hash},
            **summary,
        }
        _write_json_atomic(status_path, status)
        return status
    except Exception as exc:
        _write_json_atomic(
            status_path,
            {"status": "failed", "sample_id": sample_id, "error": str(exc)},
        )
        raise


def _run_pipeline(
    mapping: Sequence[str],
    view: Sequence[str],
    sort: Sequence[str],
    mapper_log: Path,
    samtools_log: Path,
) -> None:
    processes: list[subprocess.Popen[bytes]] = []
    with mapper_log.open("w", encoding="utf-8") as map_stderr, samtools_log.open(
        "w", encoding="utf-8"
    ) as sam_stderr:
        try:
            mapper = subprocess.Popen(mapping, stdout=subprocess.PIPE, stderr=map_stderr)
            processes.append(mapper)
            assert mapper.stdout is not None
            viewer = subprocess.Popen(
                view, stdin=mapper.stdout, stdout=subprocess.PIPE, stderr=sam_stderr
            )
            processes.append(viewer)
            mapper.stdout.close()
            assert viewer.stdout is not None
            sorter = subprocess.Popen(sort, stdin=viewer.stdout, stderr=sam_stderr)
            processes.append(sorter)
            viewer.stdout.close()
            sort_code = sorter.wait()
            view_code = viewer.wait()
            map_code = mapper.wait()
        except BaseException:
            for process in reversed(processes):
                if process.poll() is None:
                    process.terminate()
            for process in reversed(processes):
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise
    if (map_code, view_code, sort_code) != (0, 0, 0):
        raise MappingError(
            f"mapping pipeline failed: minimap2={map_code}, view={view_code}, sort={sort_code}"
        )


def _summarize_bam(samtools: str, bam: Path, log_path: Path) -> dict[str, int]:
    with log_path.open("a", encoding="utf-8") as stderr_log:
        process = subprocess.Popen(
            [samtools, "view", str(bam)],
            stdout=subprocess.PIPE,
            stderr=stderr_log,
            text=True,
        )
        assert process.stdout is not None
        try:
            summary = summarize_sam_alignments(process.stdout)
            process.stdout.close()
            process.wait()
        except BaseException:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise
    if process.returncode != 0:
        raise MappingError("samtools view summary failed; see samtools log")
    return summary


def _run_checked(command: Sequence[str], log_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=log, check=False, text=True)
    if result.returncode != 0:
        raise MappingError(f"command failed with exit code {result.returncode}: {shlex.join(command)}")


def _capture_command(command: Sequence[str], output: Path) -> None:
    temporary = output.with_suffix(f"{output.suffix}.partial")
    with temporary.open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, check=False, text=True)
    if result.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise MappingError(f"command failed: {shlex.join(command)}: {result.stderr.strip()}")
    temporary.replace(output)


def _capture_version(command: Sequence[str]) -> str:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise MappingError(f"cannot query tool version: {command[0]}")
    output = result.stdout or result.stderr
    return output.strip().splitlines()[0] if output.strip() else "unknown"


def _write_summary(path: Path, sample_id: str, summary: Mapping[str, int]) -> None:
    temporary = path.with_suffix(f"{path.suffix}.partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", *summary.keys()],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"sample_id": sample_id, **summary})
    temporary.replace(path)


def _read_json_object(path: Path, label: str) -> dict[str, object]:
    source = Path(path)
    if not source.is_file():
        raise MappingError(f"{label} is absent: {source}")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MappingError(f"{label} is not valid JSON") from exc
    if not isinstance(payload, dict):
        raise MappingError(f"{label} must contain a JSON object")
    return payload


def _read_tsv(path: Path, label: str) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file():
        raise MappingError(f"{label} is absent: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise MappingError(f"{label} has no header")
        return [dict(row) for row in reader]


def _validate_identifier(value: str, label: str) -> None:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise MappingError(f"invalid {label}: {value!r}")


def _validate_mapping_output(path: Path) -> None:
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            safe_root = (candidate / "results" / "codex_reanalysis_v1" / "02_mapping").resolve()
            if path != safe_root and safe_root not in path.parents:
                raise MappingError("mapping output must be inside 02_mapping")
            return
    raise MappingError("mapping output must be inside the project repository")


def _write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_suffix(f"{path.suffix}.partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    seal = subparsers.add_parser("seal-reference")
    seal.add_argument("--reference", required=True, type=Path)
    seal.add_argument("--targets-bed", required=True, type=Path)
    seal.add_argument("--reference-manifest", required=True, type=Path)
    seal.add_argument("--config", required=True, type=Path)
    seal.add_argument("--output", required=True, type=Path)
    mapping = subparsers.add_parser("map")
    mapping.add_argument("--sample-id", required=True)
    mapping.add_argument("--sample-registry", required=True, type=Path)
    mapping.add_argument("--input-manifest", required=True, type=Path)
    mapping.add_argument("--deliver-root", required=True, type=Path)
    mapping.add_argument("--reference", required=True, type=Path)
    mapping.add_argument("--targets-bed", required=True, type=Path)
    mapping.add_argument("--audit-status", required=True, type=Path)
    mapping.add_argument("--reference-status", required=True, type=Path)
    mapping.add_argument("--reference-manifest", required=True, type=Path)
    mapping.add_argument("--reference-seal", required=True, type=Path)
    mapping.add_argument("--config", required=True, type=Path)
    mapping.add_argument("--normalized-manifest", type=Path, default=None)
    mapping.add_argument(
        "--input-seal",
        type=Path,
        default=None,
        help="input_seal.json, so provenance can state the hash was verified on this host",
    )
    mapping.add_argument("--output-dir", required=True, type=Path)
    mapping.add_argument("--threads", type=int, default=8)
    mapping.add_argument("--minimap2", default="minimap2")
    mapping.add_argument("--samtools", default="samtools")
    mapping.add_argument("--secondary-max", type=int, default=DEFAULT_SECONDARY_MAX)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.mode == "seal-reference":
        create_reference_seal(
            args.reference, args.targets_bed, args.reference_manifest, args.config, args.output
        )
        return 0
    seal = None
    if args.input_seal is not None:
        try:
            seal = json.loads(Path(args.input_seal).read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise MappingError(f"input seal is missing or unreadable: {args.input_seal}") from error
    audited = resolve_audited_sample(
        args.sample_id,
        args.sample_registry,
        args.input_manifest,
        args.deliver_root,
        normalized_manifest_path=args.normalized_manifest,
        input_seal=seal,
    )
    fastq = Path(audited["fastq"])
    deliverable_id = str(audited["deliverable_id"])
    input_sha256 = str(audited["input_sha256"])
    run_stream_mapping(
        sample_id=args.sample_id,
        deliverable_id=deliverable_id,
        fastq=fastq,
        input_manifest=args.input_manifest,
        sample_registry=args.sample_registry,
        reference=args.reference,
        targets_bed=args.targets_bed,
        audit_status=args.audit_status,
        reference_status=args.reference_status,
        reference_manifest=args.reference_manifest,
        reference_seal=args.reference_seal,
        config_path=args.config,
        output_dir=args.output_dir,
        threads=args.threads,
        minimap2=args.minimap2,
        samtools=args.samtools,
        input_sha256=input_sha256,
        source_sha256=str(audited["source_sha256"]),
        input_normalization=str(audited["input_normalization"]),
        input_hash_verified=bool(audited["input_hash_verified_at_mapping"]),
        input_hash_verification=str(audited["input_hash_verification"]),
        secondary_max=args.secondary_max,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
