import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.codex_reanalysis.remote_preflight import (
    PreflightError,
    evaluate_seal,
    locate_delivery_root,
    recommended_jobs,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def make_delivery(tmp_path, name):
    root = tmp_path / name
    (root / "Filtered").mkdir(parents=True)
    return root


def test_delivery_root_is_found_by_its_filtered_directory(tmp_path):
    root = make_delivery(tmp_path, "mybook")

    assert locate_delivery_root([tmp_path / "absent", root]) == root.resolve()


def test_missing_delivery_root_is_a_clear_failure(tmp_path):
    with pytest.raises(PreflightError, match="no delivery root"):
        locate_delivery_root([tmp_path / "absent"])


def test_two_visible_delivery_roots_refuse_to_be_guessed(tmp_path):
    """Two mounts of the same drive would silently pick one and map from it."""
    first = make_delivery(tmp_path, "media_copy")
    second = make_delivery(tmp_path, "mnt_copy")

    with pytest.raises(PreflightError, match="refusing to guess"):
        locate_delivery_root([first, second])


def test_a_directory_without_filtered_is_not_a_delivery_root(tmp_path):
    (tmp_path / "not_delivery").mkdir()

    with pytest.raises(PreflightError):
        locate_delivery_root([tmp_path / "not_delivery"])


def test_an_unreadable_sibling_mount_is_skipped_not_a_crash(tmp_path):
    """A shared host can have another account's automount alongside the real
    one (e.g. /media/<other-user>) with no traverse permission for us. That
    must be treated as absent, not raise PermissionError."""
    blocked = tmp_path / "other_user"
    blocked.mkdir(mode=0o000)
    root = make_delivery(tmp_path, "mybook")

    try:
        assert (
            locate_delivery_root([blocked / "Deliver_ADD", root]) == root.resolve()
        )
    finally:
        blocked.chmod(0o755)


@pytest.mark.parametrize(
    "memory_gb,cpus,expected",
    [
        (64.0, 32, 3),   # the remote host: memory is the limit, not the cpus
        (64.0, 16, 2),   # cpu becomes the limit
        (128.0, 32, 4),  # more memory, but the cpu budget caps it
        (15.0, 16, 1),   # this workstation: only one job fits
        (8.0, 16, 1),    # never returns zero
        (256.0, 8, 1),   # abundant memory cannot exceed the cpu budget
    ],
)
def test_job_count_respects_both_memory_and_cpu(memory_gb, cpus, expected):
    assert recommended_jobs(memory_gb, cpus) == expected


def _seal_fixture(tmp_path):
    reference = tmp_path / "reference.fa"
    targets = tmp_path / "targets.bed"
    manifest = tmp_path / "reference_manifest.tsv"
    config = tmp_path / "config.yaml"
    seal = tmp_path / "seal.json"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    targets.write_text("chr1\t0\t4\tt\n", encoding="utf-8")
    config.write_text("mapping: {}\n", encoding="utf-8")
    manifest.write_text(
        "resource\tsha256\n"
        f"reference.fa\t{hashlib.sha256(reference.read_bytes()).hexdigest()}\n"
        f"targets.bed\t{hashlib.sha256(targets.read_bytes()).hexdigest()}\n",
        encoding="utf-8",
    )
    return reference, targets, manifest, config, seal


def test_an_untouched_seal_is_reported_valid(tmp_path):
    reference, targets, manifest, config, seal = _seal_fixture(tmp_path)
    from scripts.codex_reanalysis.stream_target_mapping import create_reference_seal

    create_reference_seal(reference, targets, manifest, config, seal)

    assert evaluate_seal(seal, reference, targets, manifest, config) == "valid"


def test_a_changed_modification_time_is_resealed_not_rejected(tmp_path):
    """Moving the exFAT drive between hosts changes mtime but not content."""
    reference, targets, manifest, config, seal = _seal_fixture(tmp_path)
    from scripts.codex_reanalysis.stream_target_mapping import create_reference_seal

    create_reference_seal(reference, targets, manifest, config, seal)
    payload = json.loads(seal.read_text(encoding="utf-8"))
    payload["resources"]["reference.fa"]["mtime_ns"] += 10_000_000_000
    seal.write_text(json.dumps(payload), encoding="utf-8")

    assert evaluate_seal(seal, reference, targets, manifest, config) == "resealed"
    assert evaluate_seal(seal, reference, targets, manifest, config) == "valid"


def test_changed_reference_content_is_a_hard_stop(tmp_path):
    reference, targets, manifest, config, seal = _seal_fixture(tmp_path)
    from scripts.codex_reanalysis.stream_target_mapping import create_reference_seal

    create_reference_seal(reference, targets, manifest, config, seal)
    reference.write_text(">chr1\nTTTT\n", encoding="utf-8")

    with pytest.raises(PreflightError, match="no longer matches the recorded checksums"):
        evaluate_seal(seal, reference, targets, manifest, config)


@pytest.mark.parametrize(
    "options,expected",
    [
        ("ro,relatime,uid=1000", True),
        ("rw,relatime,uid=1000", False),
        ("ro", True),
        ("rw,errors=remount-ro", False),   # 'remount-ro' must not read as read-only
        ("relatime,ro,nosuid", True),
    ],
)
def test_read_only_detection_reads_whole_options_not_substrings(options, expected):
    from scripts.codex_reanalysis.remote_preflight import is_read_only

    assert is_read_only(options) is expected


def test_preflight_script_offers_an_explicit_writable_override():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "remote_preflight.py"),
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert "--allow-writable-delivery" in result.stdout


def test_preflight_script_runs_as_a_direct_command():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "remote_preflight.py"),
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--project-root" in result.stdout
    assert "--deliver-root" in result.stdout


def test_preflight_seals_the_same_bed_the_workflow_maps_with():
    """A drift here re-seals the wrong file and then every mapping job fails.

    The preflight's re-seal succeeds against any BED listed in the reference
    manifest, so sealing the analytical intervals would silently destroy the
    valid seal and only surface as a failure once the cohort started.
    """
    from scripts.codex_reanalysis.coordinate_registry import (
        ANALYTICAL_BED_NAME,
        MAPPING_BED_NAME,
    )

    preflight = (
        PROJECT_ROOT / "scripts" / "codex_reanalysis" / "remote_preflight.py"
    ).read_text(encoding="utf-8")
    rules = (
        PROJECT_ROOT / "workflows" / "codex_reanalysis" / "rules" / "mapping.smk"
    ).read_text(encoding="utf-8")

    assert MAPPING_BED_NAME != ANALYTICAL_BED_NAME
    # The preflight must name the mapping BED through the shared constant.
    assert "MAPPING_BED_NAME" in preflight
    assert f'"{ANALYTICAL_BED_NAME}"' not in preflight
    # Both the seal rule and the map rule must filter on the mapping BED.
    assert rules.count("output.mapping_bed") >= 2
    assert "output.bed" not in rules


def test_preflight_invalidates_an_earlier_verdict_before_checking(tmp_path):
    """A crash mid-check must not leave a previous 'ready' for the runner.

    The runner starts a multi-day cohort on the strength of this file, so the
    window where a stale pass is readable has to be closed at the start, not
    only in the failure handler.
    """
    from scripts.codex_reanalysis.remote_preflight import main, status_path_for

    project = tmp_path / "project"
    (project / ".git").mkdir(parents=True)
    status = status_path_for(project)
    status.parent.mkdir(parents=True, exist_ok=True)
    status.write_text(json.dumps({"status": "ready", "recommended_jobs": 4}), encoding="utf-8")

    # No delivery root exists here, so the preflight must fail.
    code = main(["--project-root", str(project), "--deliver-root", str(tmp_path / "absent")])

    assert code == 1
    payload = json.loads(status.read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert "checked_at_utc" in payload
    assert "recommended_jobs" not in payload


def test_a_failed_preflight_still_hands_over_a_readonly_delivery_root(tmp_path, monkeypatch):
    """A fresh host can never reach a 'ready' preflight before
    02_verify_inputs.sh has written the input seal, because the seal is
    itself one of the checks run_preflight performs. Without this, that
    script has no way to learn the delivery root and the two scripts
    deadlock on each other. A failed verdict must still hand over a root it
    already located and confirmed read-only."""
    from scripts.codex_reanalysis import remote_preflight
    from scripts.codex_reanalysis.remote_preflight import main, status_path_for

    project = tmp_path / "project"
    (project / ".git").mkdir(parents=True)
    delivery = make_delivery(tmp_path, "mybook")
    monkeypatch.setattr(remote_preflight, "mount_options", lambda path: "ro,relatime")

    status = status_path_for(project)
    code = main(["--project-root", str(project), "--deliver-root", str(delivery)])

    assert code == 1
    payload = json.loads(status.read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert payload["delivery_root"] == str(delivery.resolve())


def test_a_failed_preflight_withholds_a_writable_delivery_root(tmp_path, monkeypatch):
    """The same bootstrap allowance must not become a way to skip the
    read-only check: a writable mount must not be handed to
    02_verify_inputs.sh even when the overall verdict already failed for an
    unrelated reason."""
    from scripts.codex_reanalysis import remote_preflight
    from scripts.codex_reanalysis.remote_preflight import main, status_path_for

    project = tmp_path / "project"
    (project / ".git").mkdir(parents=True)
    delivery = make_delivery(tmp_path, "mybook")
    monkeypatch.setattr(remote_preflight, "mount_options", lambda path: "rw,relatime")

    status = status_path_for(project)
    code = main(["--project-root", str(project), "--deliver-root", str(delivery)])

    assert code == 1
    payload = json.loads(status.read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert "delivery_root" not in payload
