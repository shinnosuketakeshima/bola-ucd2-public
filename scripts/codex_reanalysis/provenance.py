"""Create complete, machine-readable provenance records for analysis stages."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REQUIRED_FIELDS = (
    "command",
    "versions",
    "input_hashes",
    "reference_hashes",
    "config_hash",
    "git_sha",
    "git_tracked_dirty",
    "git_untracked_present",
    "timestamp_utc",
)


def write_provenance(output_dir: Path, records: dict[str, Any]) -> Path:
    """Write one sorted JSON record with mandatory rerun metadata."""
    output_path = Path(output_dir)
    _validate_output_path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    payload = _default_record()
    payload.update(records)
    for field in REQUIRED_FIELDS:
        if field not in payload:
            raise ValueError(f"provenance field is required: {field}")

    provenance_path = output_path / "provenance.json"
    temporary_path = provenance_path.with_suffix(f"{provenance_path.suffix}.partial")
    temporary_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_path.replace(provenance_path)
    return provenance_path


def _validate_output_path(output_path: Path) -> None:
    resolved = output_path.resolve()
    project_root = _project_root_for(resolved)
    if project_root is None:
        raise ValueError("provenance output must be inside a safe reanalysis output directory")
    safe_root = (project_root / "results" / "codex_reanalysis_v1").resolve()
    if resolved != safe_root and safe_root not in resolved.parents:
        raise ValueError("provenance output must be inside a safe reanalysis output directory")


def _project_root_for(path: Path) -> Path | None:
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    return None


def _default_record() -> dict[str, Any]:
    return {
        "command": " ".join(sys.argv),
        "versions": {"python": sys.version.split()[0]},
        "input_hashes": {},
        "reference_hashes": {},
        "config_hash": os.environ.get("CODEX_REANALYSIS_CONFIG_HASH", "unknown"),
        "git_sha": _git_output("rev-parse", "HEAD"),
        "git_tracked_dirty": _git_is_dirty(include_untracked=False),
        "git_untracked_present": _git_is_dirty(include_untracked=True, only_untracked=True),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def _git_output(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        capture_output=True,
        check=False,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def _git_is_dirty(*, include_untracked: bool, only_untracked: bool = False) -> bool | str:
    """Report tracked and untracked change separately.

    Release eligibility is about tracked content: this workspace deliberately
    keeps a large untracked legacy tree, so folding it into one "dirty" flag
    would mark every run ineligible for ever.
    """
    command = ["git", "status", "--porcelain"]
    command.append("--untracked-files=all" if include_untracked else "--untracked-files=no")
    result = subprocess.run(command, capture_output=True, check=False, text=True)
    if result.returncode != 0:
        return "unknown"
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if only_untracked:
        return any(line.startswith("??") for line in lines)
    return bool(lines)
