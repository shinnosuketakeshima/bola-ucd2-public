import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.codex_reanalysis.pilot_gate import (
    CHECK_NAMES,
    PilotGateError,
    SCOPE_LIMITATIONS,
    SMOKE_CHECK_NAMES,
    artefact_bindings,
    build_report,
    carrier_contrast,
    evaluate_pilot,
    parse_bedcov,
    smoke_verdicts,
    summarize_verdict,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def complete_observations(**overrides):
    observations = {
        "input_audit": True,
        "reference": True,
        "bam": True,
        "index": True,
        "coverage": 1.0,
        "control_depth": 24.5,
        "secondary_saturation": True,
        "split_reads": True,
        "assembly": True,
    }
    observations.update(overrides)
    return observations


def test_pilot_gate_accepts_complete_qc():
    assert evaluate_pilot(complete_observations()) is True


def test_pilot_gate_rejects_missing_index():
    with pytest.raises(PilotGateError, match="index"):
        evaluate_pilot(complete_observations(index=False))


def test_pilot_gate_rejects_an_absent_observation_instead_of_assuming_it_passed():
    observations = complete_observations()
    del observations["assembly"]

    with pytest.raises(PilotGateError, match="assembly"):
        evaluate_pilot(observations)


def test_pilot_gate_rejects_incomplete_target_coverage():
    with pytest.raises(PilotGateError, match="coverage"):
        evaluate_pilot(complete_observations(coverage=0.98))


def test_pilot_gate_rejects_control_depth_below_the_configured_minimum():
    with pytest.raises(PilotGateError, match="control depth"):
        evaluate_pilot(complete_observations(control_depth=9.9), minimum_depth=10)


def test_pilot_gate_accepts_control_depth_at_the_configured_minimum():
    assert evaluate_pilot(complete_observations(control_depth=10.0), minimum_depth=10) is True


def test_bedcov_mean_depth_is_summed_depth_divided_by_region_length():
    rows = parse_bedcov(
        [
            "NC_037350.1\t100\t200\tDRB3\t2500",
            "NC_037333.1\t0\t50\tcontrol_ALB\t0",
        ]
    )

    assert rows["DRB3"] == pytest.approx(25.0)
    assert rows["control_ALB"] == pytest.approx(0.0)


def test_bedcov_parser_rejects_a_zero_length_region():
    with pytest.raises(PilotGateError, match="region length"):
        parse_bedcov(["NC_037350.1\t100\t100\tDRB3\t0"])


def test_bedcov_parser_rejects_a_truncated_row():
    with pytest.raises(PilotGateError, match="bedcov"):
        parse_bedcov(["NC_037350.1\t100\t200\tDRB3"])


def test_verdict_is_passed_only_when_every_check_passed():
    assert summarize_verdict([("a", "passed"), ("b", "passed")]) == "passed"


def test_verdict_is_incomplete_when_a_stage_has_not_run():
    assert summarize_verdict([("a", "passed"), ("b", "not_run")]) == "incomplete"


def test_verdict_is_failed_when_any_check_failed_even_if_another_has_not_run():
    assert summarize_verdict([("a", "failed"), ("b", "not_run")]) == "failed"


def test_carrier_contrast_requires_one_016_carrier_and_one_non_carrier():
    registry = [
        {"sample_id": "ID02", "pcr_sbt_drb3_hap1": "*015:01", "pcr_sbt_drb3_hap2": "*016:01"},
        {"sample_id": "ID01", "pcr_sbt_drb3_hap1": "*005:02", "pcr_sbt_drb3_hap2": "*005:03"},
    ]

    contrast = carrier_contrast(registry, ["ID02", "ID01"])

    assert contrast == {"carriers": ["ID02"], "non_carriers": ["ID01"]}


def test_carrier_contrast_rejects_a_pilot_without_the_planned_contrast():
    registry = [
        {"sample_id": "ID02", "pcr_sbt_drb3_hap1": "*015:01", "pcr_sbt_drb3_hap2": "*016:01"},
        {"sample_id": "ID04", "pcr_sbt_drb3_hap1": "*015:01", "pcr_sbt_drb3_hap2": "*016:01"},
    ]

    with pytest.raises(PilotGateError, match="carrier"):
        carrier_contrast(registry, ["ID02", "ID04"])


def test_scope_limitation_is_declared_against_a_real_check_name():
    """A typo in a scope limitation would silently widen the acceptance."""
    assert set(SCOPE_LIMITATIONS) <= set(CHECK_NAMES)
    assert all(reason.strip() for reason in SCOPE_LIMITATIONS.values())


def test_verdict_is_accepted_with_scope_limitation_when_a_check_is_excluded():
    verdict = summarize_verdict(
        [("mapping", "passed"), ("assembly_smoke", "not_run")],
        excluded=("assembly_smoke",),
    )

    assert verdict == "accepted_with_scope_limitation"


def test_verdict_is_plain_passed_when_nothing_was_excluded():
    assert summarize_verdict([("mapping", "passed")], excluded=()) == "passed"


def test_an_excluded_check_cannot_mask_a_failure_in_a_required_check():
    verdict = summarize_verdict(
        [("mapping", "failed"), ("assembly_smoke", "not_run")],
        excluded=("assembly_smoke",),
    )

    assert verdict == "failed"


def test_a_required_check_that_has_not_run_still_blocks_acceptance():
    verdict = summarize_verdict(
        [("control_depth", "not_run"), ("assembly_smoke", "not_run")],
        excluded=("assembly_smoke",),
    )

    assert verdict == "incomplete"


def test_evaluate_pilot_accepts_when_the_only_unmet_check_is_excluded():
    observations = complete_observations(assembly=False)

    assert evaluate_pilot(observations, excluded=("assembly_smoke",)) is True


def test_evaluate_pilot_still_enforces_checks_outside_the_scope_limitation():
    with pytest.raises(PilotGateError, match="index"):
        evaluate_pilot(complete_observations(index=False), excluded=("assembly_smoke",))


def test_report_states_the_scope_limitation_and_its_recorded_reason():
    report = build_report(
        verdict="accepted_with_scope_limitation",
        checks=[("mapping", "passed", "2 samples"), ("assembly_smoke", "not_run", "stage absent")],
        qc_rows=[],
        contrast={"carriers": ["ID02"], "non_carriers": ["ID01"]},
        excluded={"assembly_smoke": "deferred to the assembly stage"},
    )

    assert "accepted_with_scope_limitation" in report
    assert "Scope limitation" in report
    assert "assembly_smoke" in report
    assert "deferred to the assembly stage" in report


def test_report_names_every_check_and_lists_zero_depth_targets():
    report = build_report(
        verdict="incomplete",
        checks=[("mapping", "passed", "2 samples"), ("assembly_smoke", "not_run", "stage absent")],
        qc_rows=[
            {"sample_id": "ID02", "target_name": "DRB3", "mean_depth": 31.2, "category": "class_II"},
            {"sample_id": "ID02", "target_name": "DQA5", "mean_depth": 0.0, "category": "class_II"},
        ],
        contrast={"carriers": ["ID02"], "non_carriers": ["ID01"]},
    )

    assert "incomplete" in report
    assert "assembly_smoke" in report and "not_run" in report
    assert "DQA5" in report
    # A zero-depth target must be visible as an explicit observation, because a
    # divergent haplotype and a broken pipeline look identical in a summary mean.
    assert "zero" in report.lower()


def smoke_document(bam_sha, **statuses):
    defaults = {"secondary_saturation": "passed", "sv_smoke": "passed", "assembly_smoke": "passed"}
    defaults.update(statuses)
    return {
        "bam_sha256": bam_sha,
        "checks": [{"name": n, "status": s, "detail": ""} for n, s in defaults.items()],
    }


def test_smoke_checks_pass_only_when_every_pilot_sample_passed():
    verdicts = dict(
        (name, status)
        for name, status, _ in smoke_verdicts(
            {"a": smoke_document("aa"), "b": smoke_document("bb")},
            ["a", "b"],
            {"a": "aa", "b": "bb"},
        )
    )

    assert verdicts == {
        "secondary_saturation": "passed",
        "sv_smoke": "passed",
        "assembly_smoke": "passed",
    }


def test_a_sample_without_a_smoke_document_leaves_the_check_not_run():
    """Absent evidence is not a pass, even when the other animal passed."""
    verdicts = dict(
        (name, status)
        for name, status, _ in smoke_verdicts(
            {"a": smoke_document("aa")}, ["a", "b"], {"a": "aa", "b": "bb"}
        )
    )

    assert verdicts["assembly_smoke"] == "not_run"


def test_a_smoke_document_describing_another_bam_is_a_failure_not_evidence():
    """A re-map must not inherit the smoke result of the BAM it replaced."""
    verdicts = dict(
        (name, status)
        for name, status, _ in smoke_verdicts(
            {"a": smoke_document("old-digest")}, ["a"], {"a": "new-digest"}
        )
    )

    assert set(verdicts.values()) == {"failed"}


def test_one_failing_sample_fails_the_check_for_the_pilot():
    verdicts = dict(
        (name, status)
        for name, status, _ in smoke_verdicts(
            {"a": smoke_document("aa"), "b": smoke_document("bb", sv_smoke="failed")},
            ["a", "b"],
            {"a": "aa", "b": "bb"},
        )
    )

    assert verdicts["sv_smoke"] == "failed"
    assert verdicts["assembly_smoke"] == "passed"


def test_no_mapped_samples_yields_not_run_rather_than_a_vacuous_pass():
    verdicts = dict((name, status) for name, status, _ in smoke_verdicts({}, [], {}))

    assert set(verdicts.values()) == {"not_run"}


def test_the_smoke_checks_are_declared_in_the_gate():
    assert set(SMOKE_CHECK_NAMES) <= set(CHECK_NAMES)


def test_nothing_is_excluded_from_acceptance_by_default():
    """The assembly exclusion was retired when the smoke stage landed; leaving
    it in place would be a permanent allowance rather than a scope note."""
    assert SCOPE_LIMITATIONS == {}


def test_the_verdict_records_the_artefacts_it_was_reached_on(tmp_path):
    root = tmp_path / "results"
    (root / "02_mapping").mkdir(parents=True)
    (root / "01_reference").mkdir(parents=True)
    (root / "02_mapping" / "a.targets.bam.bai").write_bytes(b"index")
    (root / "01_reference" / "targets.bed").write_text("x\t1\t2\tt\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("outputs: {}\n", encoding="utf-8")

    bindings = artefact_bindings(root, {"a": "digest-of-a"}, config)

    assert bindings["bam"] == {"a": "digest-of-a"}
    assert len(bindings["bai"]["a"]) == 64
    assert "targets.bed" in bindings["reference"]
    assert bindings["config_sha256"] != "not_supplied"


def test_an_unsupplied_config_is_recorded_as_such_not_as_a_hash(tmp_path):
    root = tmp_path / "results"
    (root / "02_mapping").mkdir(parents=True)
    (root / "01_reference").mkdir(parents=True)

    assert artefact_bindings(root, {}, None)["config_sha256"] == "not_supplied"


def test_pilot_gate_runs_as_a_direct_command():
    result = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "pilot_gate.py"), "--help"],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--root" in result.stdout
    assert "--expected" in result.stdout
    assert "--report" in result.stdout
