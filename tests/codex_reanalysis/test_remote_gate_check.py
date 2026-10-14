from datetime import datetime, timedelta, timezone

import pytest

from scripts.codex_reanalysis.remote_gate_check import (
    GateError,
    check_pilot_gate,
    check_preflight,
)


NOW = datetime(2026, 8, 5, 12, 0, tzinfo=timezone.utc)


def preflight(**overrides):
    status = {
        "status": "ready",
        "checked_at_utc": (NOW - timedelta(minutes=30)).isoformat(),
        "delivery_root": "/mnt/mybook/Deliver_ADD",
        "recommended_jobs": 3,
    }
    status.update(overrides)
    return status


def gate(**overrides):
    status = {
        "accepted": True,
        "status": "passed",
        "excluded_from_acceptance": {},
        "checks": [
            {"name": "mapping", "status": "passed"},
            {"name": "secondary_saturation", "status": "passed"},
            {"name": "assembly_smoke", "status": "passed"},
            {"name": "sv_smoke", "status": "passed"},
        ],
    }
    status.update(overrides)
    return status


def test_a_fresh_ready_preflight_yields_its_verified_settings():
    assert check_preflight(preflight(), NOW) == ("/mnt/mybook/Deliver_ADD", 3)


def test_a_failed_preflight_is_refused():
    with pytest.raises(GateError, match="did not pass"):
        check_preflight(preflight(status="failed", error="mount is writable"), NOW)


def test_a_running_preflight_is_refused():
    """An interrupted preflight leaves 'running'; that is not a pass."""
    with pytest.raises(GateError, match="did not pass"):
        check_preflight({"status": "running", "checked_at_utc": NOW.isoformat()}, NOW)


def test_a_stale_preflight_is_refused():
    with pytest.raises(GateError, match="old"):
        check_preflight(preflight(checked_at_utc=(NOW - timedelta(hours=9)).isoformat()), NOW)


def test_an_unstamped_preflight_is_refused():
    status = preflight()
    del status["checked_at_utc"]
    with pytest.raises(GateError, match="timestamp"):
        check_preflight(status, NOW)


def test_a_future_stamp_is_refused():
    with pytest.raises(GateError, match="future"):
        check_preflight(preflight(checked_at_utc=(NOW + timedelta(hours=2)).isoformat()), NOW)


def test_a_fully_passed_pilot_is_accepted():
    check_pilot_gate(gate())


def test_a_scope_limited_pilot_cannot_start_the_cohort():
    """This is exactly the state the pilot is in today: accepted, but with the
    assembly smoke never run. Days of cohort compute must not follow it."""
    with pytest.raises(GateError, match="accepted_with_scope_limitation|excluded"):
        check_pilot_gate(
            gate(
                status="accepted_with_scope_limitation",
                excluded_from_acceptance={"assembly_smoke": "no stage yet"},
                checks=[
                    {"name": "mapping", "status": "passed"},
                    {"name": "assembly_smoke", "status": "not_run"},
                ],
            )
        )


def test_a_pilot_missing_the_sv_smoke_is_refused():
    with pytest.raises(GateError, match="sv_smoke"):
        check_pilot_gate(
            gate(checks=[{"name": "mapping", "status": "passed"},
                         {"name": "secondary_saturation", "status": "passed"},
                         {"name": "assembly_smoke", "status": "passed"}])
        )


def test_a_pilot_whose_multi_hit_evidence_was_truncated_is_refused():
    """A spike at the -N ceiling means alignments were still being discarded;
    the cohort would inherit censored multi-hit evidence."""
    with pytest.raises(GateError, match="secondary_saturation"):
        check_pilot_gate(
            gate(checks=[{"name": "secondary_saturation", "status": "failed"},
                         {"name": "assembly_smoke", "status": "passed"},
                         {"name": "sv_smoke", "status": "passed"}])
        )


def test_a_failed_smoke_is_refused_even_when_the_verdict_says_passed():
    with pytest.raises(GateError, match="assembly_smoke"):
        check_pilot_gate(
            gate(checks=[{"name": "secondary_saturation", "status": "passed"},
                         {"name": "assembly_smoke", "status": "failed"},
                         {"name": "sv_smoke", "status": "passed"}])
        )


def test_an_unaccepted_pilot_is_refused():
    with pytest.raises(GateError, match="did not accept"):
        check_pilot_gate(gate(accepted=False, blocking_reason="control depth too low"))
