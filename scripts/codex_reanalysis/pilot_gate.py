"""Accept or reject the two-animal pilot before the cohort is expanded.

The gate never infers a passing state from absent evidence: a check whose
stage has not run is reported as ``not_run`` and keeps the pilot out of a
``passed`` verdict, so an unimplemented stage can never be mistaken for a
satisfied one.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from hashlib import sha256
from pathlib import Path
from typing import Iterable, Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance


class PilotGateError(RuntimeError):
    """Raised when the pilot cannot be accepted from the evidence on disk."""


FOCAL_ALLELE = "*016:01"
BOOLEAN_OBSERVATIONS = (
    "input_audit",
    "reference",
    "bam",
    "index",
    "secondary_saturation",
    "split_reads",
    "assembly",
)
REQUIRED_OBSERVATIONS = (*BOOLEAN_OBSERVATIONS, "coverage", "control_depth")
CONTROL_CATEGORY = "control_candidate"
DEFAULT_MINIMUM_DEPTH = 10

CHECK_NAMES = (
    "input_audit",
    "reference_registry",
    "mapping",
    "alignment_index",
    "target_coverage",
    "control_depth",
    "secondary_saturation",
    "sv_smoke",
    "assembly_smoke",
    "cohort_contrast",
)

#: Checks deliberately placed outside the pilot acceptance scope, with the
#: reason recorded here rather than in a run-time flag, so the decision is
#: reviewable in version control.
#:
#: This is empty on purpose. It previously excluded ``assembly_smoke`` because
#: no stage produced it; ``scripts/codex_reanalysis/pilot_smoke.py`` now does,
#: so the exclusion was retired rather than left standing as a permanent
#: allowance. The mechanism stays because a future exclusion must be declared
#: here, with its reason, and not passed as a run-time flag.
SCOPE_LIMITATIONS: dict[str, str] = {}

#: Observation key for each declared check, used to skip an excluded check.
_OBSERVATION_BY_CHECK = {
    "input_audit": "input_audit",
    "reference_registry": "reference",
    "mapping": "bam",
    "alignment_index": "index",
    "target_coverage": "coverage",
    "control_depth": "control_depth",
    "secondary_saturation": "secondary_saturation",
    "sv_smoke": "split_reads",
    "assembly_smoke": "assembly",
}


def _excluded_observations(excluded: Sequence[str]) -> set[str]:
    unknown = sorted(set(excluded) - set(CHECK_NAMES))
    if unknown:
        raise PilotGateError(f"unknown check excluded from acceptance: {', '.join(unknown)}")
    return {_OBSERVATION_BY_CHECK[name] for name in excluded if name in _OBSERVATION_BY_CHECK}


def evaluate_pilot(
    observations: Mapping[str, object],
    minimum_depth: int = DEFAULT_MINIMUM_DEPTH,
    excluded: Sequence[str] = (),
) -> bool:
    """Return True only when every in-scope pilot observation is satisfied."""
    skipped = _excluded_observations(excluded)
    for name in REQUIRED_OBSERVATIONS:
        if name in skipped:
            continue
        if name not in observations:
            raise PilotGateError(f"pilot observation is absent: {name}")
    for name in BOOLEAN_OBSERVATIONS:
        if name in skipped:
            continue
        if not observations[name]:
            raise PilotGateError(f"pilot check did not pass: {name}")
    if "coverage" not in skipped:
        coverage = float(observations["coverage"])  # type: ignore[arg-type]
        if coverage < 1.0:
            raise PilotGateError(
                f"target coverage is incomplete: {coverage:.3f} of declared targets were measured"
            )
    if "control_depth" not in skipped:
        control_depth = float(observations["control_depth"])  # type: ignore[arg-type]
        if control_depth < minimum_depth:
            raise PilotGateError(
                f"control depth {control_depth:.2f} is below the configured minimum {minimum_depth}"
            )
    return True


def parse_bedcov(lines: Iterable[str]) -> dict[str, float]:
    """Convert `samtools bedcov` rows into mean depth per named target."""
    depths: dict[str, float] = {}
    for line in lines:
        row = line.rstrip("\n\r")
        if not row:
            continue
        fields = row.split("\t")
        if len(fields) < 5:
            raise PilotGateError(f"malformed bedcov row: {row!r}")
        name = fields[3]
        try:
            start, end = int(fields[1]), int(fields[2])
            summed = float(fields[4])
        except ValueError as error:
            raise PilotGateError(f"malformed bedcov row: {row!r}") from error
        length = end - start
        if length <= 0:
            raise PilotGateError(f"bedcov region length is not positive for {name}")
        depths[name] = summed / length
    return depths


def carrier_contrast(
    registry: Sequence[Mapping[str, str]], pilot_samples: Sequence[str]
) -> dict[str, list[str]]:
    """Confirm the pilot contrasts a focal-allele carrier with a non-carrier."""
    by_sample = {row["sample_id"]: row for row in registry}
    carriers: list[str] = []
    non_carriers: list[str] = []
    for sample_id in pilot_samples:
        row = by_sample.get(sample_id)
        if row is None:
            raise PilotGateError(f"pilot sample is absent from the registry: {sample_id}")
        diplotype = {row.get("pcr_sbt_drb3_hap1", ""), row.get("pcr_sbt_drb3_hap2", "")}
        (carriers if FOCAL_ALLELE in diplotype else non_carriers).append(sample_id)
    if not carriers or not non_carriers:
        raise PilotGateError(
            "pilot must contrast one "
            f"{FOCAL_ALLELE} carrier with one non-carrier; "
            f"observed {len(carriers)} carriers and {len(non_carriers)} non-carriers"
        )
    return {"carriers": carriers, "non_carriers": non_carriers}


def summarize_verdict(
    checks: Sequence[tuple[str, str]] | Sequence[tuple[str, str, str]],
    excluded: Sequence[str] = (),
) -> str:
    """Reduce individual check statuses to one verdict without optimism.

    An excluded check is reported but cannot decide the verdict, and its
    exclusion is never hidden: the result is downgraded from ``passed`` to
    ``accepted_with_scope_limitation``.
    """
    _excluded_observations(excluded)
    in_scope = [check[1] for check in checks if check[0] not in set(excluded)]
    if not in_scope or "failed" in in_scope:
        return "failed"
    if "not_run" in in_scope:
        return "incomplete"
    return "accepted_with_scope_limitation" if excluded else "passed"


def build_report(
    *,
    verdict: str,
    checks: Sequence[tuple[str, str, str]],
    qc_rows: Sequence[Mapping[str, object]],
    contrast: Mapping[str, Sequence[str]],
    excluded: Mapping[str, str] | None = None,
) -> str:
    """Render a reviewable pilot report that names every unmet condition."""
    lines = [
        "# ARS-UCD2.0 two-animal pilot report",
        "",
        f"**Verdict: {verdict}**",
        "",
        "A verdict of `incomplete` means no check failed but at least one stage "
        "has not run; it is not an acceptance. A verdict of "
        "`accepted_with_scope_limitation` means the in-scope checks passed while "
        "the conditions listed below were deliberately excluded, and it must not "
        "be reported as an unqualified pass.",
        "",
    ]
    if excluded:
        lines.extend(["## Scope limitation", ""])
        lines.extend(
            [
                "The following declared conditions were excluded from this acceptance. "
                "They remain unverified.",
                "",
                "| Excluded check | Recorded reason |",
                "| --- | --- |",
                *(f"| {name} | {reason} |" for name, reason in sorted(excluded.items())),
                "",
            ]
        )
    lines += [
        "## Cohort contrast",
        "",
        f"- {FOCAL_ALLELE} carriers: {', '.join(contrast.get('carriers', [])) or 'none'}",
        f"- non-carriers: {', '.join(contrast.get('non_carriers', [])) or 'none'}",
        "",
        "## Checks",
        "",
        "| Check | Status | Detail |",
        "| --- | --- | --- |",
    ]
    lines.extend(f"| {name} | `{status}` | {detail} |" for name, status, detail in checks)

    zero_depth = [row for row in qc_rows if float(row["mean_depth"]) == 0.0]
    lines.extend(
        [
            "",
            "## Targets with zero measured depth",
            "",
            "A divergent haplotype and a broken pipeline both produce zero depth, "
            "so these are listed individually for scientific review rather than "
            "averaged away.",
            "",
        ]
    )
    if zero_depth:
        lines.extend(
            [
                "| Sample | Target | Category |",
                "| --- | --- | --- |",
                *(
                    f"| {row['sample_id']} | {row['target_name']} | {row['category']} |"
                    for row in zero_depth
                ),
            ]
        )
    else:
        lines.append("No target had zero measured depth.")

    lines.extend(["", "## Per-target mean depth", "", "| Sample | Target | Category | Mean depth |", "| --- | --- | --- | --- |"])
    lines.extend(
        f"| {row['sample_id']} | {row['target_name']} | {row['category']} | {float(row['mean_depth']):.2f} |"
        for row in qc_rows
    )
    return "\n".join(lines) + "\n"


def read_tsv(path: Path, label: str) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file():
        raise PilotGateError(f"{label} is absent: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise PilotGateError(f"{label} has no header")
        return [dict(row) for row in reader]


def read_json(path: Path, label: str) -> dict[str, object]:
    source = Path(path)
    if not source.is_file():
        raise PilotGateError(f"{label} is absent: {source}")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise PilotGateError(f"{label} is not valid JSON: {source}") from error
    if not isinstance(payload, dict):
        raise PilotGateError(f"{label} must contain a JSON object")
    return payload


def run_bedcov(samtools: str, bed: Path, bam: Path) -> list[str]:
    """Read summed depth per target region straight from the pilot BAM."""
    result = subprocess.run(
        [samtools, "bedcov", str(bed), str(bam)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise PilotGateError(f"samtools bedcov failed: {result.stderr.strip()}")
    return result.stdout.splitlines()


SMOKE_CHECK_NAMES = ("secondary_saturation", "sv_smoke", "assembly_smoke")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def smoke_verdicts(
    documents: Mapping[str, Mapping[str, object]],
    samples: Sequence[str],
    bam_digests: Mapping[str, str],
) -> list[tuple[str, str, str]]:
    """Fold per-sample smoke documents into one verdict per declared check.

    A check passes only when every mapped pilot sample passed it on the BAM
    this gate just hashed. A missing document is ``not_run``, never a pass, and
    a document describing some other BAM is a failure rather than evidence.
    """
    results: list[tuple[str, str, str]] = []
    for name in SMOKE_CHECK_NAMES:
        statuses: list[str] = []
        details: list[str] = []
        for sample in samples:
            document = documents.get(sample)
            if document is None:
                statuses.append("not_run")
                details.append(f"{sample}=absent")
                continue
            recorded = document.get("bam_sha256")
            expected_digest = bam_digests.get(sample)
            if expected_digest is None or recorded != expected_digest:
                statuses.append("failed")
                details.append(f"{sample}=stale (smoke describes a different BAM)")
                continue
            by_name = {
                str(check.get("name")): (str(check.get("status")), str(check.get("detail", "")))
                for check in document.get("checks", [])
                if isinstance(check, Mapping)
            }
            status, detail = by_name.get(name, ("not_run", "check absent from the smoke status"))
            statuses.append(status)
            details.append(f"{sample}={status}: {detail}")
        if not statuses:
            verdict = "not_run"
        elif "failed" in statuses:
            verdict = "failed"
        elif set(statuses) == {"passed"}:
            verdict = "passed"
        else:
            verdict = "not_run"
        results.append((name, verdict, "; ".join(details) or "no mapped pilot samples"))
    return results


def collect_smoke_checks(
    mapping_dir: Path, samples: Sequence[str], bam_digests: Mapping[str, str]
) -> list[tuple[str, str, str]]:
    documents: dict[str, Mapping[str, object]] = {}
    for sample in samples:
        path = mapping_dir / "pilot_smoke" / f"{sample}.smoke_status.json"
        if path.is_file():
            documents[sample] = read_json(path, "smoke status")
    return smoke_verdicts(documents, samples, bam_digests)


def collect_evidence(
    root: Path, expected: int, samtools: str, minimum_depth: int
) -> tuple[list[tuple[str, str, str]], list[dict[str, object]], dict[str, list[str]], dict[str, str]]:
    """Gather every observation the gate can make from the outputs on disk."""
    audit_dir = root / "00_input_audit"
    reference_dir = root / "01_reference"
    mapping_dir = root / "02_mapping"
    checks: list[tuple[str, str, str]] = []

    audit = read_json(audit_dir / "audit_status.json", "input audit status")
    audit_ok = audit.get("status") == "pass" and audit.get("analytical_sample_count") == 13
    checks.append(
        (
            "input_audit",
            "passed" if audit_ok else "failed",
            f"status={audit.get('status')}, analytical={audit.get('analytical_sample_count')}",
        )
    )

    reference = read_json(reference_dir / "reference_status.json", "reference status")
    reference_ok = reference.get("status") == "passed"
    checks.append(
        (
            "reference_registry",
            "passed" if reference_ok else "failed",
            f"status={reference.get('status')}, targets={reference.get('target_count')}",
        )
    )

    statuses = sorted(mapping_dir.glob("*.mapping_status.json"))
    samples = [path.name.removesuffix(".mapping_status.json") for path in statuses]
    mapping_ok = len(samples) == expected
    detail = f"{len(samples)} mapped sample(s), expected {expected}"
    for path in statuses:
        payload = read_json(path, "mapping status")
        if payload.get("status") != "passed":
            mapping_ok = False
            detail += f"; {payload.get('sample_id')} status={payload.get('status')}"
    checks.append(("mapping", "passed" if mapping_ok else "failed", detail))

    index_ok = bool(samples)
    index_detail = []
    bam_digests: dict[str, str] = {}
    for sample in samples:
        bam = mapping_dir / f"{sample}.targets.bam"
        bai = mapping_dir / f"{sample}.targets.bam.bai"
        present = bam.is_file() and bam.stat().st_size > 0 and bai.is_file() and bai.stat().st_size > 0
        if present:
            check = subprocess.run(
                [samtools, "quickcheck", "-v", str(bam)], capture_output=True, text=True, check=False
            )
            present = check.returncode == 0
        if present:
            # The verdict has to name the artefacts it was reached on, or a
            # later re-map silently inherits an acceptance it never earned.
            bam_digests[sample] = sha256_file(bam)
        index_ok = index_ok and present
        index_detail.append(f"{sample}={'ok' if present else 'missing'}")
    checks.append(
        ("alignment_index", "passed" if index_ok else "failed", ", ".join(index_detail) or "no samples")
    )

    registry_rows = read_tsv(reference_dir / "coordinate_registry.tsv", "coordinate registry")
    category_by_target = {row["target_name"]: row["category"] for row in registry_rows}
    declared = {
        row["target_name"] for row in registry_rows if row["target_type"] == "reference_feature"
    }

    qc_rows: list[dict[str, object]] = []
    measured_fractions: list[float] = []
    control_depths: list[float] = []
    bed = reference_dir / "targets.bed"
    for sample in samples:
        bam = mapping_dir / f"{sample}.targets.bam"
        if not bam.is_file():
            continue
        depths = parse_bedcov(run_bedcov(samtools, bed, bam))
        measured_fractions.append(len(declared & depths.keys()) / len(declared) if declared else 0.0)
        for target_name, mean_depth in sorted(depths.items()):
            category = category_by_target.get(target_name, "unknown")
            qc_rows.append(
                {
                    "sample_id": sample,
                    "target_name": target_name,
                    "category": category,
                    "mean_depth": round(mean_depth, 4),
                }
            )
            if category == CONTROL_CATEGORY:
                control_depths.append(mean_depth)

    coverage = min(measured_fractions) if measured_fractions else 0.0
    checks.append(
        (
            "target_coverage",
            "passed" if coverage >= 1.0 else "failed",
            f"{coverage:.3f} of {len(declared)} declared reference-feature targets measured",
        )
    )

    control_depth = min(control_depths) if control_depths else 0.0
    checks.append(
        (
            "control_depth",
            "passed" if control_depths and control_depth >= minimum_depth else "failed",
            f"lowest single-copy control mean depth {control_depth:.2f} (minimum {minimum_depth})",
        )
    )

    checks.extend(collect_smoke_checks(mapping_dir, samples, bam_digests))

    contrast: dict[str, list[str]] = {"carriers": [], "non_carriers": []}
    if samples:
        audit_registry = read_tsv(audit_dir / "sample_registry.tsv", "sample registry")
        try:
            contrast = carrier_contrast(audit_registry, samples)
            checks.append(
                (
                    "cohort_contrast",
                    "passed",
                    f"carriers={','.join(contrast['carriers'])}; non_carriers={','.join(contrast['non_carriers'])}",
                )
            )
        except PilotGateError as error:
            checks.append(("cohort_contrast", "failed", str(error)))
    else:
        checks.append(("cohort_contrast", "failed", "no mapped pilot samples"))

    return checks, qc_rows, contrast, bam_digests


#: Reference artefacts the pilot verdict is only meaningful against. If any of
#: these changes, the cohort is no longer running what the pilot accepted.
BOUND_REFERENCE_FILES = (
    "mapping_reference_seal.json",
    "reference_manifest.tsv",
    "targets.bed",
    "mapping_targets.bed",
)


def artefact_bindings(
    root: Path, bam_digests: Mapping[str, str], config: Path | None
) -> dict[str, object]:
    """Record the digests of everything the verdict was reached on."""
    mapping_dir = root / "02_mapping"
    reference_dir = root / "01_reference"
    bai: dict[str, str] = {}
    for sample in bam_digests:
        index = mapping_dir / f"{sample}.targets.bam.bai"
        if index.is_file():
            bai[sample] = sha256_file(index)
    reference: dict[str, str] = {}
    for name in BOUND_REFERENCE_FILES:
        path = reference_dir / name
        if path.is_file():
            reference[name] = sha256_file(path)
    return {
        "bam": dict(bam_digests),
        "bai": bai,
        "reference": reference,
        "config_sha256": sha256_file(config) if config and config.is_file() else "not_supplied",
    }


def write_qc_table(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.partial")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sample_id", "target_name", "category", "mean_depth"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f"{path.suffix}.partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--expected", required=True, type=int)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--minimum-depth", type=int, default=DEFAULT_MINIMUM_DEPTH)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="run configuration to bind into the verdict, so a cohort cannot inherit it",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="require every declared check, ignoring the recorded scope limitations",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    root = args.root.resolve()
    limitations: dict[str, str] = {} if args.strict else dict(SCOPE_LIMITATIONS)
    excluded = tuple(limitations)
    checks, qc_rows, contrast, bam_digests = collect_evidence(
        root, args.expected, args.samtools, args.minimum_depth
    )
    verdict = summarize_verdict(checks, excluded=excluded)

    by_name = {name: status for name, status, _ in checks}
    accepted = False
    reason = ""
    try:
        accepted = evaluate_pilot(
            {
                "input_audit": by_name.get("input_audit") == "passed",
                "reference": by_name.get("reference_registry") == "passed",
                "bam": by_name.get("mapping") == "passed",
                "index": by_name.get("alignment_index") == "passed",
                "coverage": 1.0 if by_name.get("target_coverage") == "passed" else 0.0,
                "control_depth": (
                    args.minimum_depth if by_name.get("control_depth") == "passed" else 0.0
                ),
                "secondary_saturation": by_name.get("secondary_saturation") == "passed",
                "split_reads": by_name.get("sv_smoke") == "passed",
                "assembly": by_name.get("assembly_smoke") == "passed",
            },
            minimum_depth=args.minimum_depth,
            excluded=excluded,
        )
    except PilotGateError as error:
        reason = str(error)

    bindings = artefact_bindings(root, bam_digests, args.config)
    write_qc_table(root / "02_mapping" / "pilot_qc.tsv", qc_rows)
    _write_json_atomic(
        root / "00_input_audit" / "pilot_status.json",
        {
            "status": verdict,
            "accepted": accepted,
            "acceptance_scope": "strict" if args.strict else "declared_scope_limitations",
            "excluded_from_acceptance": limitations,
            "blocking_reason": reason,
            "expected_samples": args.expected,
            "checks": [{"name": n, "status": s, "detail": d} for n, s, d in checks],
            "cohort_contrast": contrast,
            # A verdict that does not name the artefacts it was reached on can
            # be inherited by a later re-map that never earned it.
            "verified_artefacts": bindings,
        },
    )
    report = build_report(
        verdict=verdict, checks=checks, qc_rows=qc_rows, contrast=contrast, excluded=limitations
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report, encoding="utf-8")
    write_provenance(
        root / "02_mapping" / "pilot_provenance",
        {
            "command": f"pilot_gate.py --root [root] --expected {args.expected}",
            "versions": {"python": sys.version.split()[0]},
            "input_hashes": dict(bindings.get("bam", {})),
            "reference_hashes": dict(bindings.get("reference", {})),
            "config_hash": bindings.get("config_sha256", "not_applicable"),
        },
    )

    sys.stderr.write(f"pilot verdict: {verdict}\n")
    if reason:
        sys.stderr.write(f"blocking: {reason}\n")
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
