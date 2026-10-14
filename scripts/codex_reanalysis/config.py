"""Load the portable configuration without exposing delivery data paths."""

from __future__ import annotations

from dataclasses import dataclass
import csv
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

import yaml


PROTECTED_REPOSITORY_PATHS = (Path("data/raw"), Path("data/meta/private"))


@dataclass(frozen=True)
class RuntimeConfig:
    """Validated settings for one reanalysis run.

    Delivery data are always supplied at runtime and never written into the
    portable configuration file.
    """

    project_root: Path
    deliver_root: Path
    config_path: Path
    output_root: Path
    sample_sheet_path: Path
    reference_fasta_path: Path
    tools: Mapping[str, str]
    reference: Mapping[str, Any]
    sample_sheet_schema: tuple[str, ...]
    thresholds: Mapping[str, Any]
    config_sha256: str


def load_runtime_config(deliver_root: Path | None, config_path: Path) -> RuntimeConfig:
    """Load immutable settings while requiring an external delivery root."""
    if deliver_root is None:
        raise ValueError("deliver root must be supplied explicitly")

    resolved_config = Path(config_path).resolve()
    project_root = _project_root_for(resolved_config)
    _validate_config_path(resolved_config, project_root)
    if not resolved_config.is_file():
        raise FileNotFoundError(f"configuration file does not exist: {resolved_config}")
    resolved_deliver_root = Path(deliver_root).resolve()
    _validate_deliver_root(resolved_deliver_root, project_root)

    raw_text = resolved_config.read_text(encoding="utf-8")
    parsed = yaml.safe_load(raw_text)
    if not isinstance(parsed, dict):
        raise ValueError("configuration must contain a mapping")

    tools = _required_mapping(parsed, "tools")
    reference = _required_mapping(parsed, "reference")
    samples = _required_mapping(parsed, "samples")
    outputs = _required_mapping(parsed, "outputs")
    thresholds = _required_mapping(parsed, "thresholds")
    _require_keys(reference, "assembly", "accession", "fasta_relative_path")
    _require_keys(samples, "sheet_relative_path", "sample_id_column")

    sample_schema = tuple(
        str(value) for key, value in samples.items() if key.endswith("_column") and value
    )
    sample_sheet_path = _resolve_sample_sheet_path(
        project_root,
        samples["sheet_relative_path"],
        sample_schema,
    )
    return RuntimeConfig(
        project_root=project_root,
        deliver_root=resolved_deliver_root,
        config_path=resolved_config,
        output_root=_resolve_output_root(project_root, outputs.get("root")),
        sample_sheet_path=sample_sheet_path,
        reference_fasta_path=_resolve_reference_fasta_path(
            project_root,
            reference["fasta_relative_path"],
        ),
        tools=_freeze_mapping(tools),
        reference=_freeze_mapping(reference),
        sample_sheet_schema=sample_schema,
        thresholds=_freeze_mapping(thresholds),
        config_sha256=sha256(raw_text.encode("utf-8")).hexdigest(),
    )


def _project_root_for(config_path: Path) -> Path:
    for candidate in (config_path.parent, *config_path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise ValueError("configuration must be located inside a Git repository")


def _validate_deliver_root(deliver_root: Path, project_root: Path) -> None:
    if not deliver_root.is_dir():
        raise ValueError("deliver root must be an existing directory")
    if deliver_root == project_root or project_root in deliver_root.parents:
        protected = any(
            deliver_root == (project_root / relative).resolve()
            or (project_root / relative).resolve() in deliver_root.parents
            for relative in PROTECTED_REPOSITORY_PATHS
        )
        if protected:
            raise ValueError("deliver root cannot resolve to repository protected data")
        raise ValueError("deliver root must be outside the repository")


def _validate_config_path(config_path: Path, project_root: Path) -> None:
    if not _is_within(config_path, project_root):
        raise ValueError("configuration path must remain inside the repository")
    if _is_protected_path(config_path, project_root):
        raise ValueError("configuration path cannot resolve to repository protected data")


def _required_mapping(config: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = config.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"configuration section {key!r} must be a mapping")
    return value


def _require_keys(section: Mapping[str, Any], *keys: str) -> None:
    missing = [key for key in keys if not section.get(key)]
    if missing:
        raise ValueError(f"missing required configuration values: {', '.join(missing)}")


def _resolve_output_root(project_root: Path, configured_root: Any) -> Path:
    if not isinstance(configured_root, str) or not configured_root:
        raise ValueError("outputs.root must be a non-empty relative path")
    candidate = Path(configured_root)
    if candidate.is_absolute():
        raise ValueError("outputs.root must be project-root-relative")
    resolved = (project_root / candidate).resolve()
    if not _is_within(resolved, project_root):
        raise ValueError("outputs.root must remain inside the repository")
    if _is_protected_path(resolved, project_root):
        raise ValueError("outputs.root cannot resolve to repository protected data")
    safe_root = (project_root / "results" / "codex_reanalysis_v1").resolve()
    if not _is_within(resolved, safe_root):
        raise ValueError("outputs.root must be inside results/codex_reanalysis_v1")
    return resolved


def _resolve_sample_sheet_path(
    project_root: Path,
    configured_path: Any,
    required_columns: tuple[str, ...],
) -> Path:
    if not isinstance(configured_path, str) or not configured_path:
        raise ValueError("samples.sheet_relative_path must be a non-empty relative path")
    candidate = Path(configured_path)
    if candidate.is_absolute():
        raise ValueError("samples.sheet_relative_path must be project-root-relative")
    resolved = (project_root / candidate).resolve()
    if not _is_within(resolved, project_root):
        raise ValueError("samples.sheet_relative_path must remain inside the repository")
    if _is_protected_path(resolved, project_root):
        raise ValueError("samples.sheet_relative_path cannot resolve to repository protected data")
    if not resolved.is_file():
        raise ValueError("tracked sample-sheet schema does not exist")
    with resolved.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    missing = [column for column in required_columns if column not in header]
    if missing:
        raise ValueError(f"tracked sample-sheet schema is missing columns: {', '.join(missing)}")
    return resolved


def _resolve_reference_fasta_path(project_root: Path, configured_path: Any) -> Path:
    if not isinstance(configured_path, str) or not configured_path:
        raise ValueError("reference.fasta_relative_path must be a non-empty relative path")
    candidate = Path(configured_path)
    if candidate.is_absolute():
        raise ValueError("reference.fasta_relative_path must be project-root-relative")
    resolved = (project_root / candidate).resolve()
    if not _is_within(resolved, project_root):
        raise ValueError("reference.fasta_relative_path must remain inside the repository")
    if _is_protected_path(resolved, project_root):
        raise ValueError("reference.fasta_relative_path cannot resolve to repository protected data")
    return resolved


def _is_within(path: Path, root: Path) -> bool:
    return path == root or root in path.parents


def _is_protected_path(path: Path, project_root: Path) -> bool:
    return any(_is_within(path, (project_root / relative).resolve()) for relative in PROTECTED_REPOSITORY_PATHS)


def _freeze_mapping(value: Mapping[str, Any]) -> Mapping[str, Any]:
    frozen = {
        str(key): _freeze_mapping(item) if isinstance(item, dict) else tuple(item) if isinstance(item, list) else item
        for key, item in value.items()
    }
    return MappingProxyType(frozen)
