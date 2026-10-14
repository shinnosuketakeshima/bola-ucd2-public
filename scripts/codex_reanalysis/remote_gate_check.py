"""Decide whether this host may start a long run, and emit the verified settings.

The shell runners must not re-read the status after checking it: everything a
runner needs is read once here, validated, and printed, so the values used are
the values that were verified.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.input_seal import InputSealError, validate_input_seal
# Imported rather than repeated: a check renamed in the gate must not silently
# stop being required here, which would reopen the cohort on weaker evidence.
from scripts.codex_reanalysis.pilot_gate import SMOKE_CHECK_NAMES


class GateError(RuntimeError):
    """Raised when the host or the pilot does not justify starting."""


MAX_PREFLIGHT_AGE_HOURS = 6.0

#: A cohort may only follow a pilot that passed every declared check. A verdict
#: of accepted_with_scope_limitation means conditions were deliberately left
#: unverified, which is not a basis for committing days of compute.
REQUIRED_GATE_STATUS = "passed"
REQUIRED_SMOKE_CHECKS = SMOKE_CHECK_NAMES


def check_preflight(
    status: Mapping[str, object],
    now: datetime,
    max_age_hours: float = MAX_PREFLIGHT_AGE_HOURS,
) -> tuple[str, int]:
    """Validate a preflight verdict and return the delivery root and job count."""
    state = status.get("status")
    if state != "ready":
        raise GateError(f"preflight did not pass: {status.get('error') or state}")
    stamped = status.get("checked_at_utc")
    if not isinstance(stamped, str) or not stamped:
        raise GateError("preflight status has no timestamp; re-run the preflight")
    age = (now - datetime.fromisoformat(stamped)).total_seconds() / 3600
    if age < 0:
        raise GateError("preflight status is stamped in the future; check the clock")
    if age > max_age_hours:
        raise GateError(
            f"preflight is {age:.1f} h old; re-run it against the current machine state"
        )
    delivery = status.get("delivery_root")
    jobs = status.get("recommended_jobs")
    if not isinstance(delivery, str) or not delivery:
        raise GateError("preflight status has no delivery root")
    if not isinstance(jobs, int) or jobs < 1:
        raise GateError("preflight status has no usable job count")
    return delivery, jobs


def check_pilot_gate(gate: Mapping[str, object]) -> None:
    """Refuse a cohort unless the pilot passed without scope limitations."""
    if not gate.get("accepted"):
        raise GateError(
            f"the pilot gate did not accept: {gate.get('blocking_reason') or gate.get('status')}"
        )
    if gate.get("status") != REQUIRED_GATE_STATUS:
        raise GateError(
            f"the pilot verdict is {gate.get('status')!r}, not {REQUIRED_GATE_STATUS!r}; "
            "conditions left out of scope must be verified before the cohort"
        )
    excluded = gate.get("excluded_from_acceptance") or {}
    if excluded:
        raise GateError(
            "the pilot was accepted with conditions excluded: " + ", ".join(sorted(excluded))
        )
    by_name = {
        str(check.get("name")): str(check.get("status"))
        for check in gate.get("checks", [])
        if isinstance(check, Mapping)
    }
    for required in REQUIRED_SMOKE_CHECKS:
        if by_name.get(required) != "passed":
            raise GateError(
                f"{required} is {by_name.get(required, 'absent')!r}; the cohort needs it passed"
            )


def check_mount_still_read_only(delivery: str, findmnt: str = "findmnt") -> str:
    """Confirm the delivery mount is read-only now, not merely when checked."""
    result = subprocess.run(
        [findmnt, "-no", "OPTIONS", "--target", delivery],
        capture_output=True,
        text=True,
        check=False,
    )
    options = result.stdout.strip().splitlines()
    observed = options[0] if options else ""
    if "ro" not in {token.strip() for token in observed.split(",")}:
        raise GateError(f"the delivery mount is no longer read-only ({observed!r}); remount it")
    return observed


def _read_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise GateError(f"{label} is missing or unreadable: {path}") from error
    if not isinstance(payload, dict):
        raise GateError(f"{label} must contain a JSON object: {path}")
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-status", required=True, type=Path)
    parser.add_argument("--pilot-status", type=Path, default=None)
    parser.add_argument("--require-pilot-gate", action="store_true")
    parser.add_argument("--input-seal", type=Path, default=None)
    parser.add_argument("--input-manifest", type=Path, default=None)
    parser.add_argument("--findmnt", default="findmnt")
    args = parser.parse_args(argv)

    try:
        status = _read_json(args.preflight_status, "preflight status")
        delivery, jobs = check_preflight(status, datetime.now(timezone.utc))
        options = check_mount_still_read_only(delivery, args.findmnt)
        if args.input_seal is not None:
            if args.input_manifest is None:
                raise GateError("an input manifest is required to validate the input seal")
            try:
                validate_input_seal(_read_json(args.input_seal, "input seal"), args.input_manifest)
            except InputSealError as error:
                raise GateError(f"input seal is not usable: {error}") from error
        if args.require_pilot_gate:
            if args.pilot_status is None:
                raise GateError("a pilot status path is required to start the cohort")
            check_pilot_gate(_read_json(args.pilot_status, "pilot status"))
    except GateError as error:
        sys.stderr.write(f"REFUSED: {error}\n")
        return 1

    sys.stderr.write(f"host verified; delivery mounted [{options}]\n")
    print(delivery)
    print(jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
