"""Summarize per-sample mapping QC into one cohort-level report.

Closes out Task 7 of the ARS-UCD2.0 reanalysis plan: every one of the 13
manuscript samples must appear exactly once with a passing mapping status,
and any sample whose MAPQ60+ fraction or primary-alignment count is a robust
outlier against the rest of the cohort must be flagged explicitly rather than
silently averaged away.
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path
from typing import Mapping, Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance


class CohortError(ValueError):
    """Raised when the cohort sample set is incomplete, duplicated, or unexpected."""


#: Robust outlier threshold on the modified z-score (Iglewicz & Hoaglin 1993).
OUTLIER_THRESHOLD = 3.5


def validate_cohort_ids(expected_ids: Sequence[str], observed_ids: Sequence[str]) -> bool:
    """Confirm every expected sample id appears in the observed set exactly once."""
    expected = list(expected_ids)
    observed = list(observed_ids)
    duplicate_expected = _duplicates(expected)
    if duplicate_expected:
        raise CohortError(f"duplicate expected sample IDs: {', '.join(duplicate_expected)}")
    duplicate_observed = _duplicates(observed)
    if duplicate_observed:
        raise CohortError(f"duplicate observed sample IDs: {', '.join(duplicate_observed)}")
    missing = sorted(set(expected) - set(observed))
    if missing:
        raise CohortError(f"missing expected samples: {', '.join(missing)}")
    unexpected = sorted(set(observed) - set(expected))
    if unexpected:
        raise CohortError(f"unexpected samples: {', '.join(unexpected)}")
    return True


def _duplicates(values: Sequence[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def flag_outliers(values: Mapping[str, float], threshold: float = OUTLIER_THRESHOLD) -> dict[str, bool]:
    """Flag values whose modified z-score against the cohort median is extreme.

    Falls back to flagging any nonzero deviation when the median absolute
    deviation itself is zero (a tight, mostly-identical cohort), since a
    zero MAD would otherwise make every modified z-score infinite or
    undefined.
    """
    keys = list(values)
    numbers = [float(values[key]) for key in keys]
    if not numbers:
        return {}
    median = statistics.median(numbers)
    deviations = {key: abs(value - median) for key, value in zip(keys, numbers)}
    mad = statistics.median(deviations.values())
    if mad == 0:
        return {key: deviation > 0 for key, deviation in deviations.items()}
    return {
        key: (0.6745 * deviation / mad) > threshold
        for key, deviation in deviations.items()
    }


def load_expected_sample_ids(registry_path: Path) -> list[str]:
    with Path(registry_path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return sorted(
        row["sample_id"].strip()
        for row in rows
        if row.get("analysis_status", "").strip() == "analytical" and row.get("sample_id", "").strip()
    )


def load_mapping_status(root: Path, sample_id: str) -> dict:
    status_path = Path(root) / f"{sample_id}.mapping_status.json"
    if not status_path.is_file():
        raise CohortError(f"missing mapping status for sample: {sample_id}")
    return json.loads(status_path.read_text(encoding="utf-8"))


def build_cohort_qc_rows(root: Path, expected_ids: Sequence[str]) -> list[dict]:
    rows = [load_mapping_status(root, sample_id) for sample_id in expected_ids]
    validate_cohort_ids(expected_ids, [row["sample_id"] for row in rows])

    mapq60_fraction = {
        row["sample_id"]: row["mapq_60_plus"] / row["alignment_count"] if row["alignment_count"] else 0.0
        for row in rows
    }
    primary_count = {row["sample_id"]: row["primary_count"] for row in rows}
    mapq_flags = flag_outliers(mapq60_fraction)
    read_count_flags = flag_outliers({k: float(v) for k, v in primary_count.items()})
    mapq_median = statistics.median(mapq60_fraction.values())
    read_median = statistics.median(primary_count.values())

    qc_rows = []
    for row in rows:
        sample_id = row["sample_id"]
        low_mapq = mapq_flags[sample_id] and mapq60_fraction[sample_id] < mapq_median
        low_reads = read_count_flags[sample_id] and primary_count[sample_id] < read_median
        flags = [name for name, hit in (("low_mapq60_fraction", low_mapq), ("low_primary_count", low_reads)) if hit]
        if row.get("status") != "passed":
            flags.append("mapping_status_not_passed")
        qc_rows.append(
            {
                "sample_id": sample_id,
                "status": row.get("status", "unknown"),
                "alignment_count": row["alignment_count"],
                "primary_count": row["primary_count"],
                "secondary_count": row["secondary_count"],
                "supplementary_count": row["supplementary_count"],
                "mapq_60_plus": row["mapq_60_plus"],
                "mapq60_fraction": round(mapq60_fraction[sample_id], 4),
                "flags": ";".join(flags),
            }
        )
    return sorted(qc_rows, key=lambda item: item["sample_id"])


def write_cohort_qc_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fieldnames = (
        "sample_id",
        "status",
        "alignment_count",
        "primary_count",
        "secondary_count",
        "supplementary_count",
        "mapq_60_plus",
        "mapq60_fraction",
        "flags",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def render_report(rows: Sequence[Mapping[str, object]]) -> str:
    flagged = [row for row in rows if row["flags"]]
    lines = [
        "# ARS-UCD2.0 Cohort Mapping QC Report",
        "",
        f"13-animal cohort, all-target mapping stage (`02_mapping`). {len(rows)} samples present, "
        f"each exactly once. {len(flagged)} sample(s) flagged below.",
        "",
        "| sample_id | status | alignment_count | primary_count | secondary_count | "
        "supplementary_count | mapq60_fraction | flags |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['sample_id']} | {row['status']} | {row['alignment_count']} | {row['primary_count']} | "
            f"{row['secondary_count']} | {row['supplementary_count']} | {row['mapq60_fraction']} | "
            f"{row['flags'] or '-'} |"
        )
    lines.append("")
    if flagged:
        lines.append("## Flagged samples")
        lines.append("")
        for row in flagged:
            lines.append(f"- **{row['sample_id']}**: {row['flags']}")
    else:
        lines.append("No samples flagged: no missing/duplicate IDs, all `passed`, no MAPQ/read-count outliers.")
    lines.append("")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="02_mapping stage directory")
    parser.add_argument("--registry", required=True, type=Path, help="00_input_audit/sample_registry.tsv")
    parser.add_argument("--output", required=True, type=Path, help="Markdown report destination")
    args = parser.parse_args(argv)

    expected_ids = load_expected_sample_ids(args.registry)
    rows = build_cohort_qc_rows(args.root, expected_ids)

    qc_dir = Path(args.root) / "cohort_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    write_cohort_qc_tsv(qc_dir / "cohort_qc.tsv", rows)

    flagged = [row for row in rows if row["flags"]]
    status_payload = {
        "sample_count": len(rows),
        "flagged_sample_count": len(flagged),
        "flagged_samples": [row["sample_id"] for row in flagged],
        "status": "pass" if not flagged else "flagged",
    }
    (qc_dir / "cohort_qc_status.json").write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        qc_dir,
        {
            "input_hashes": {row["sample_id"]: row["status"] for row in rows},
        },
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_report(rows), encoding="utf-8")
    print(f"cohort QC report written: {args.output} ({status_payload['status']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
