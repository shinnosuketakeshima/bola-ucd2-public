"""Check a remote host is ready to map the cohort, before any long job starts.

The delivery drive is exFAT and is moved between machines, so the reference
seal's cheap staleness guard (size and modification time) is expected to fail
even when the bytes are identical. This tool distinguishes the two cases: a
changed modification time is re-sealed, and re-sealing can only succeed if the
content still matches the SHA-256 the reference stage recorded. Changed content
is a hard stop.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.coordinate_registry import MAPPING_BED_NAME
from scripts.codex_reanalysis.input_seal import InputSealError, validate_input_seal
from scripts.codex_reanalysis.stream_target_mapping import (
    MappingError,
    create_reference_seal,
    validate_reference_seal,
    validate_stage_gates,
)


class PreflightError(RuntimeError):
    """Raised when the host cannot be trusted to run the cohort."""


#: Observed peak resident memory of one minimap2 job against ARS-UCD2.0 with
#: the configured secondary-alignment cap of 20. Measured at 14.0 GB and
#: still rising slowly, so the figure is rounded up rather than fitted.
MAPPING_MEMORY_GB = 14.5
#: Headroom for the operating system, page cache, samtools and sort buffers.
#: Deliberately generous: an out-of-memory kill ten hours into a twelve-hour
#: job costs far more than one fewer concurrent sample.
RESERVED_MEMORY_GB = 8.0
THREADS_PER_JOB = 8


def locate_delivery_root(candidates: Iterable[Path]) -> Path:
    """Return the single readable delivery root among the given candidates."""
    found = []
    for candidate in candidates:
        path = Path(candidate)
        try:
            is_delivery_root = path.is_dir() and (path / "Filtered").is_dir()
        except PermissionError:
            # Another account's automount can sit alongside the real one
            # (e.g. /media/<other-user> on a shared desktop) with no
            # traverse permission for us. That is not a candidate, not a
            # crash.
            continue
        if is_delivery_root:
            found.append(path.resolve())
    unique = sorted(set(found))
    if not unique:
        raise PreflightError(
            "no delivery root found; mount the drive and pass --deliver-root explicitly"
        )
    if len(unique) > 1:
        raise PreflightError(
            "several delivery roots are visible, refusing to guess: "
            + ", ".join(str(path) for path in unique)
        )
    return unique[0]


def default_candidates(explicit: Path | None) -> list[Path]:
    """Build the search list for an auto-mounted or hand-mounted drive."""
    if explicit is not None:
        return [Path(explicit)]
    roots = [Path("/media"), Path("/mnt"), Path("/run/media")]
    candidates: list[Path] = []
    for root in roots:
        if not root.is_dir():
            continue
        for entry in _safe_iterdir(root):
            candidates.append(entry / "Deliver_ADD")
            for nested in _safe_iterdir(entry):
                candidates.append(nested / "Deliver_ADD")
    return candidates


def _safe_iterdir(path: Path) -> list[Path]:
    try:
        return [entry for entry in path.iterdir() if entry.is_dir()]
    except (PermissionError, OSError):
        return []


def mount_options(path: Path, findmnt: str = "findmnt") -> str:
    """Return the mount options for the filesystem holding a path."""
    result = subprocess.run(
        [findmnt, "-no", "OPTIONS", "--target", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise PreflightError(f"cannot determine mount options for {path}")
    return result.stdout.strip().splitlines()[0]


def is_read_only(options: str) -> bool:
    """Decide whether a mount option string marks the filesystem read-only."""
    return "ro" in {token.strip() for token in options.split(",")}


def recommended_jobs(available_memory_gb: float, cpu_count: int) -> int:
    """Choose how many samples may run at once without risking the machine."""
    by_memory = int((available_memory_gb - RESERVED_MEMORY_GB) // MAPPING_MEMORY_GB)
    by_cpu = cpu_count // THREADS_PER_JOB
    return max(1, min(by_memory, by_cpu))


def evaluate_seal(
    seal_path: Path,
    reference: Path,
    targets_bed: Path,
    manifest_path: Path,
    config_path: Path,
) -> str:
    """Return 'valid', or re-seal and return 'resealed'; raise if content moved."""
    try:
        validate_reference_seal(seal_path, reference, targets_bed, manifest_path, config_path)
        return "valid"
    except MappingError:
        pass
    try:
        # Only succeeds when every byte still matches the recorded SHA-256.
        create_reference_seal(reference, targets_bed, manifest_path, config_path, seal_path)
    except MappingError as error:
        raise PreflightError(
            f"reference content no longer matches the recorded checksums: {error}"
        ) from error
    return "resealed"


def read_json(path: Path, label: str) -> dict[str, object]:
    source = Path(path)
    if not source.is_file():
        raise PreflightError(f"{label} is absent: {source}")
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise PreflightError(f"{label} must contain a JSON object")
    return payload


def check_tool(command: Sequence[str], expected: str, label: str) -> str:
    if shutil.which(command[0]) is None:
        raise PreflightError(f"{label} is not on PATH: {command[0]}")
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    output = (result.stdout or result.stderr).strip().splitlines()
    observed = output[0] if output else ""
    if expected not in observed:
        raise PreflightError(f"{label} version mismatch: expected {expected}, found {observed!r}")
    return observed


def available_memory_gb() -> float:
    """Read memory the kernel says is actually available, not merely free."""
    meminfo = Path("/proc/meminfo")
    if not meminfo.is_file():
        return 0.0
    for line in meminfo.read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / (1024 * 1024)
    return 0.0


def run_preflight(
    project_root: Path,
    deliver_root: Path | None,
    config_rel: str,
    allow_writable: bool = False,
) -> dict[str, object]:
    root = Path(project_root).resolve()
    if not (root / ".git").exists():
        raise PreflightError(f"project root has no .git, provenance cannot record a commit: {root}")
    output = root / "results" / "codex_reanalysis_v1"
    audit = output / "00_input_audit"
    normalized = output / "00b_normalized"
    reference = output / "01_reference"
    config = root / config_rel

    checks: list[dict[str, object]] = []

    delivery = locate_delivery_root(default_candidates(deliver_root))
    options = mount_options(delivery)
    if not is_read_only(options) and not allow_writable:
        raise PreflightError(
            f"the delivery mount is writable ({options}); remount it read-only, or pass "
            "--allow-writable-delivery to accept the risk to irreplaceable data"
        )
    checks.append(
        {
            "name": "delivery_root",
            "status": "ok",
            "detail": f"{delivery} [{'ro' if is_read_only(options) else 'RW ACCEPTED'}]",
        }
    )

    audit_status = read_json(audit / "audit_status.json", "input audit status")
    normalization = read_json(normalized / "normalization_status.json", "normalization status")
    validate_stage_gates(
        audit / "audit_status.json",
        audit / "input_manifest.tsv",
        audit / "sample_registry.tsv",
        reference / "reference_status.json",
        reference / "reference_manifest.tsv",
    )
    checks.append(
        {
            "name": "stage_gates",
            "status": "ok",
            "detail": (
                f"audit={audit_status.get('status')}, "
                f"analytical={audit_status.get('analytical_sample_count')}, "
                f"normalised={normalization.get('filtered_sample_count')}"
            ),
        }
    )

    # Content verification of the delivered FASTQ is expensive and lives in its
    # own long-lived seal; this check only confirms it exists and still
    # describes the current audit manifest.
    seal_file = audit / "input_seal.json"
    if not seal_file.is_file():
        raise PreflightError(
            "delivered inputs have not been verified on this host; run "
            "run_remote/02_verify_inputs.sh first"
        )
    try:
        validate_input_seal(read_json(seal_file, "input seal"), audit / "input_manifest.tsv")
    except InputSealError as error:
        raise PreflightError(f"input seal is not usable: {error}") from error
    checks.append(
        {
            "name": "input_seal",
            "status": "ok",
            "detail": f"{read_json(seal_file, 'input seal').get('sample_count')} FASTQ verified by content",
        }
    )

    seal_state = evaluate_seal(
        reference / "mapping_reference_seal.json",
        reference / "GCF_002263795.3_ARS-UCD2.0_genomic.fna",
        reference / MAPPING_BED_NAME,
        reference / "reference_manifest.tsv",
        config,
    )
    checks.append({"name": "reference_seal", "status": "ok", "detail": seal_state})

    versions = {
        "minimap2": check_tool(["minimap2", "--version"], "2.31-r1302", "minimap2"),
        "samtools": check_tool(["samtools", "--version"], "1.23.1", "samtools"),
        "pigz": check_tool(["pigz", "--version"], "2.8", "pigz"),
    }
    checks.append({"name": "tool_versions", "status": "ok", "detail": json.dumps(versions)})

    memory = available_memory_gb()
    cpus = os.cpu_count() or 1
    jobs = recommended_jobs(memory, cpus)
    free_gb = shutil.disk_usage(root).free / (1024 ** 3)
    checks.append(
        {
            "name": "resources",
            "status": "ok",
            "detail": (
                f"available_memory={memory:.1f} GB, cpus={cpus}, "
                f"project_disk_free={free_gb:.1f} GB, recommended_jobs={jobs}"
            ),
        }
    )
    if free_gb < 30:
        raise PreflightError(f"project disk has only {free_gb:.1f} GB free; 30 GB is the minimum")

    status = {
        "status": "ready",
        "delivery_mount_options": options,
        "project_root": str(root),
        "delivery_root": str(delivery),
        "reference_seal": seal_state,
        "recommended_jobs": jobs,
        "cores_for_snakemake": jobs * THREADS_PER_JOB,
        "available_memory_gb": round(memory, 1),
        "cpu_count": cpus,
        "checks": checks,
    }
    write_status(output / "preflight_status.json", status)
    return status


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--deliver-root", type=Path, default=None)
    parser.add_argument("--config", default="configs/codex_reanalysis.example.yaml")
    parser.add_argument(
        "--allow-writable-delivery",
        action="store_true",
        help="proceed even if the delivery mount is writable",
    )
    return parser.parse_args(argv)


def status_path_for(project_root: Path) -> Path:
    return Path(project_root) / "results" / "codex_reanalysis_v1" / "preflight_status.json"


def write_status(path: Path, payload: dict[str, object]) -> None:
    """Replace the status atomically so no reader sees a half-written verdict."""
    path.parent.mkdir(parents=True, exist_ok=True)
    stamped = {**payload, "checked_at_utc": datetime.now(timezone.utc).isoformat()}
    temporary = path.with_suffix(f"{path.suffix}.partial")
    temporary.write_text(json.dumps(stamped, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    status_file = status_path_for(args.project_root)
    # Invalidate any earlier verdict before checking anything. A crash, a kill
    # or a power loss part-way through must not leave a previous "ready" behind
    # for the cohort runner to act on.
    write_status(status_file, {"status": "running"})
    try:
        status = run_preflight(
            args.project_root, args.deliver_root, args.config, args.allow_writable_delivery
        )
    except (PreflightError, MappingError) as error:
        failed_status: dict[str, object] = {"status": "failed", "error": str(error)}
        # A fresh host can never reach "ready" before 02_verify_inputs.sh has
        # written the input seal, because that seal is itself one of the
        # checks below. Without this, that script has no way to learn the
        # delivery root and the two scripts deadlock on each other. A
        # read-only root located here is real evidence, not a guess, so it
        # is handed over even though the overall verdict is a failure.
        try:
            delivery = locate_delivery_root(default_candidates(args.deliver_root))
            if is_read_only(mount_options(delivery)):
                failed_status["delivery_root"] = str(delivery)
        except (PreflightError, OSError):
            pass
        write_status(status_file, failed_status)
        sys.stderr.write(f"PREFLIGHT FAILED: {error}\n")
        return 1
    for check in status["checks"]:
        sys.stderr.write(f"  {check['name']:18s} {check['status']:6s} {check['detail']}\n")
    sys.stderr.write(
        "\nPREFLIGHT READY\n"
        f"  delivery root   : {status['delivery_root']}\n"
        f"  reference seal  : {status['reference_seal']}\n"
        f"  recommended     : {status['recommended_jobs']} concurrent samples "
        f"(--cores {status['cores_for_snakemake']})\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
