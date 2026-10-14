"""Build a privacy-aware, stable registry for delivered sequencing inputs."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from pathlib import Path, PureWindowsPath
from typing import Iterable, Mapping, Sequence


ABSENT_DELIVERABLES = ("DELIVERY_07", "DELIVERY_10")


class SampleRegistryError(ValueError):
    """Raised when manuscript samples do not have a one-to-one FASTQ mapping."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated metadata table without emitting its contents."""
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_sample_registry(
    rows: Iterable[Mapping[str, str]], expected_ids: Sequence[str]
) -> list[dict[str, str]]:
    """Validate and normalize one delivered FASTQ per manuscript sample.

    Input rows must contain ``sample_id`` and ``fastq``.  The function adds the
    two documented delivery-number gaps as non-analytical rows.  It deliberately
    does not infer a sample identity when metadata are incomplete.
    """
    expected = tuple(expected_ids)
    if not expected or len(set(expected)) != len(expected):
        raise SampleRegistryError("expected sample IDs must be unique and non-empty")

    grouped: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        fastq = str(row.get("fastq", "")).strip()
        if sample_id and fastq:
            grouped[sample_id].append(row)

    missing = sorted(set(expected) - set(grouped))
    extra = sorted(set(grouped) - set(expected))
    multiple = sorted(sample for sample, values in grouped.items() if len(values) != 1)
    if missing:
        raise SampleRegistryError(f"missing expected samples: {', '.join(missing)}")
    if extra:
        raise SampleRegistryError(f"unexpected samples: {', '.join(extra)}")
    if multiple:
        raise SampleRegistryError(f"multiple inputs for samples: {', '.join(multiple)}")

    registry: list[dict[str, str]] = []
    for sample_id in expected:
        source = grouped[sample_id][0]
        registry.append(
            {
                "deliverable_id": str(source.get("deliverable_id", "")).strip(),
                "sample_id": sample_id,
                "manuscript_id": str(source.get("manuscript_id", sample_id)).strip(),
                "input_relative_path": str(source["fastq"]).replace("\\", "/"),
                "analysis_status": "analytical",
                "pcr_sbt_status": "external_metadata",
                "pcr_sbt_drb3_hap1": str(source.get("drb3_hap1", "")).strip(),
                "pcr_sbt_drb3_hap2": str(source.get("drb3_hap2", "")).strip(),
            }
        )
    present_deliverables = {item["deliverable_id"] for item in registry}
    for deliverable_id in ABSENT_DELIVERABLES:
        if deliverable_id not in present_deliverables:
            registry.append(
                {
                    "deliverable_id": deliverable_id,
                    "sample_id": "",
                    "manuscript_id": "",
                    "input_relative_path": "",
                    "analysis_status": "absent_cohort",
                    "pcr_sbt_status": "not_applicable",
                    "pcr_sbt_drb3_hap1": "",
                    "pcr_sbt_drb3_hap2": "",
                }
            )
    return sorted(registry, key=lambda row: row["deliverable_id"])


def build_rows_from_metadata(
    samplesheet_path: Path,
    deliverable_mapping_path: Path,
    drb3_path: Path,
    deliver_root: Path,
) -> tuple[list[dict[str, str]], list[str]]:
    """Join tracked metadata and resolve only relative paths below deliver_root."""
    samples = read_tsv(samplesheet_path)
    mappings = read_tsv(deliverable_mapping_path)
    drb3_by_sample = {row["sample_id"]: row for row in read_tsv(drb3_path)}
    mapping_by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
    for mapping in mappings:
        if mapping.get("fastq_path", "").strip():
            mapping_by_sample[mapping.get("sample_id", "").strip()].append(mapping)

    expected_ids = [row["sample_id"].strip() for row in samples]
    duplicate_expected = [name for name, count in Counter(expected_ids).items() if count > 1]
    if duplicate_expected:
        raise SampleRegistryError(
            f"duplicate sample IDs in samplesheet: {', '.join(sorted(duplicate_expected))}"
        )

    joined: list[dict[str, str]] = []
    root = Path(deliver_root).resolve()
    for sample in samples:
        sample_id = sample["sample_id"].strip()
        candidates = mapping_by_sample.get(sample_id, [])
        if len(candidates) != 1:
            # Let the shared one-to-one validation give a stable error message.
            continue
        mapping = candidates[0]
        # The source metadata records Windows drive paths.  PureWindowsPath
        # accepts both separators, so the delivered file name is recovered
        # identically on Windows and on the POSIX host that runs mapping; a
        # POSIX parser would read the whole backslash string as one name and
        # silently drop every sample.
        supplied = PureWindowsPath(mapping["fastq_path"])
        # The source metadata may contain an old absolute drive path.  Bind the
        # audit strictly to the explicit current delivery root by filename.
        candidate = (root / "Filtered" / supplied.name).resolve()
        if root not in candidate.parents or not candidate.is_file():
            continue
        relative = candidate.relative_to(root).as_posix()
        typed = drb3_by_sample.get(sample_id, {})
        joined.append(
            {
                "deliverable_id": mapping["deliverable_id"].strip(),
                "sample_id": sample_id,
                "manuscript_id": sample_id,
                "fastq": relative,
                "drb3_hap1": typed.get("drb3_hap1", ""),
                "drb3_hap2": typed.get("drb3_hap2", ""),
            }
        )
    return joined, expected_ids


def write_tsv(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    """Write a deterministic TSV table, including headers for an empty result."""
    fieldnames = (
        "deliverable_id",
        "sample_id",
        "manuscript_id",
        "input_relative_path",
        "analysis_status",
        "pcr_sbt_status",
        "pcr_sbt_drb3_hap1",
        "pcr_sbt_drb3_hap2",
    )
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
