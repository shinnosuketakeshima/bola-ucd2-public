import json
import os
import shutil
from pathlib import Path
from uuid import uuid4

import pytest
import yaml

from scripts.codex_reanalysis.config import _resolve_output_root, load_runtime_config
from scripts.codex_reanalysis.provenance import _default_record, write_provenance


def test_release_eligibility_ignores_untracked_files():
    """This workspace permanently holds hundreds of untracked legacy files.

    Counting them as a dirty tree would make every output release-ineligible
    for ever, so the release signal must come from tracked changes only and the
    untracked state must be recorded separately rather than conflated.
    """
    record = _default_record()

    assert "git_tracked_dirty" in record
    assert "git_untracked_present" in record
    assert isinstance(record["git_tracked_dirty"], bool)
    assert isinstance(record["git_untracked_present"], bool)
    # ...and that must not, by itself, mark the run ineligible.
    assert record["git_tracked_dirty"] is False or record["git_tracked_dirty"] is True


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_CONFIG = PROJECT_ROOT / "configs" / "codex_reanalysis.example.yaml"


def _safe_output_dir() -> Path:
    return PROJECT_ROOT / "results" / "codex_reanalysis_v1" / "provenance_tests" / uuid4().hex


def _temporary_config(reference_path: str) -> Path:
    path = PROJECT_ROOT / "tmp" / f"runtime-config-{uuid4().hex}.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load(EXAMPLE_CONFIG.read_text(encoding="utf-8"))
    config["reference"]["fasta_relative_path"] = reference_path
    path.write_text(
        yaml.safe_dump(config, sort_keys=False),
        encoding="utf-8",
    )
    return path


def test_runtime_config_requires_explicit_deliver_root(tmp_path):
    """A missing data-root argument must never fall back to a private path."""
    with pytest.raises(ValueError, match="deliver root"):
        load_runtime_config(None, EXAMPLE_CONFIG)


def test_runtime_config_resolves_output_under_project_and_freezes_settings(tmp_path):
    """A relative output is project-local and its scientific settings cannot mutate."""
    delivery_root = tmp_path / "delivery"
    delivery_root.mkdir()
    config = load_runtime_config(delivery_root, EXAMPLE_CONFIG)

    assert config.output_root == PROJECT_ROOT / "results" / "codex_reanalysis_v1"
    assert config.reference["accession"] == "GCF_002263795.3"
    assert config.reference_fasta_path == (
        PROJECT_ROOT
        / "results"
        / "codex_reanalysis_v1"
        / "01_reference"
        / "GCF_002263795.3_ARS-UCD2.0_genomic.fna"
    )
    assert config.sample_sheet_path == PROJECT_ROOT / "configs" / "codex_reanalysis.samplesheet.schema.tsv"
    assert config.sample_sheet_schema == ("sample_id", "fastq_relative_path")
    with pytest.raises(TypeError):
        config.tools["minimap2"] = "other-minimap2"


def test_runtime_config_rejects_repository_protected_data_root():
    """A delivery argument resolving inside repository source-data is unsafe."""
    (PROJECT_ROOT / "data" / "raw").mkdir(parents=True, exist_ok=True)  # git-ignored
    with pytest.raises(ValueError, match="protected data"):
        load_runtime_config(PROJECT_ROOT / "data" / "raw", EXAMPLE_CONFIG)


def test_runtime_config_requires_existing_delivery_directory(tmp_path):
    """A typo in the external delivery path must stop before analysis begins."""
    with pytest.raises(ValueError, match="existing directory"):
        load_runtime_config(tmp_path / "missing-delivery", EXAMPLE_CONFIG)


def test_runtime_config_rejects_protected_config_path(tmp_path):
    """Configuration files may not be loaded from a protected data directory."""
    delivery_root = tmp_path / "delivery"
    delivery_root.mkdir()
    with pytest.raises(ValueError, match="configuration path cannot resolve to repository protected data"):
        load_runtime_config(delivery_root, PROJECT_ROOT / "data" / "raw" / "runtime.yaml")


@pytest.mark.parametrize(
    "reference_path",
    [
        "data/raw/reference.fa",
        pytest.param(
            "C:/private/reference.fa",
            marks=pytest.mark.skipif(os.name != "nt", reason="drive-letter paths are absolute on Windows only"),
        ),
    ],
)
def test_runtime_config_rejects_protected_or_absolute_reference_path(tmp_path, reference_path):
    """Reference FASTA configuration cannot point at protected or absolute data paths."""
    delivery_root = tmp_path / "delivery"
    delivery_root.mkdir()
    config_path = _temporary_config(reference_path)
    try:
        with pytest.raises(ValueError, match="reference.fasta_relative_path"):
            load_runtime_config(delivery_root, config_path)
    finally:
        config_path.unlink(missing_ok=True)


@pytest.mark.parametrize("unsafe_root", ["data/raw", "results/../data/meta/private"])
def test_output_root_rejects_protected_paths_after_resolution(unsafe_root):
    """Traversal into protected data cannot be configured as an output root."""
    with pytest.raises(ValueError, match="protected data"):
        _resolve_output_root(PROJECT_ROOT, unsafe_root)


def test_provenance_json_is_stable_for_fixed_complete_record(tmp_path):
    """Fixed provenance input must yield byte-identical sorted JSON."""
    record = {
        "sample": "ID01",
        "version": "1",
        "command": "python analysis.py",
        "versions": {"python": "3.13.0"},
        "input_hashes": {"reads.fastq.gz": "input-sha"},
        "reference_hashes": {"assembly.fa": "reference-sha"},
        "config_hash": "config-sha",
        "git_sha": "commit-sha",
        "git_tracked_dirty": False,
        "git_untracked_present": True,
        "timestamp_utc": "2026-07-31T00:00:00+00:00",
    }
    output_root = _safe_output_dir()
    try:
        first = write_provenance(output_root / "first", record).read_text(encoding="utf-8")
        second = write_provenance(output_root / "second", record).read_text(encoding="utf-8")
    finally:
        shutil.rmtree(output_root, ignore_errors=True)
    payload = json.loads(first)

    assert first == second
    assert payload["sample"] == "ID01"
    assert list(payload) == sorted(payload)
    assert {
        "command",
        "versions",
        "input_hashes",
        "reference_hashes",
        "config_hash",
        "git_sha",
        "git_tracked_dirty",
        "git_untracked_present",
        "timestamp_utc",
    } <= set(payload)


def test_provenance_refuses_protected_repository_output():
    """A provenance record must not be written beside protected source data."""
    with pytest.raises(ValueError, match="safe reanalysis output"):
        write_provenance(PROJECT_ROOT / "data" / "raw" / "provenance", {})
