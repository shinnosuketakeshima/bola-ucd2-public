"""Write public DDBJ BioSample / DRA metadata worksheets for the 13-animal cohort.

Public identifiers are manuscript IDs ID01–ID13 only. Geographic origin is
reported as Japan. Real deliverable names, prefectures, and the anonymization
key must never appear in these files.

This script does not upload data, contact DDBJ, or mint accession numbers.
Fill accession_register.tsv after D-way issues PRJDB / SAMD / DRX / DRR
accessions, then replace the placeholders in Declarations.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

if __package__ in {None, ""}:
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# Table 1 inferred_sex. ID13 remains uncertain in the manuscript; BioSample uses
# INSDC missing-value language rather than inventing a sex call.
COHORT: tuple[dict[str, str], ...] = (
    {"sample_id": "ID01", "sex": "female", "drb3": "*005:02/*005:03", "batch": "2021"},
    {"sample_id": "ID02", "sex": "male", "drb3": "*015:01/*016:01", "batch": "2021"},
    {"sample_id": "ID03", "sex": "male", "drb3": "*002:01/*015:01", "batch": "2021"},
    {"sample_id": "ID04", "sex": "male", "drb3": "*015:01/*016:01", "batch": "2021"},
    {"sample_id": "ID05", "sex": "male", "drb3": "*011:01/*016:01", "batch": "2021"},
    {"sample_id": "ID06", "sex": "male", "drb3": "*015:01/*016:01", "batch": "2021"},
    {"sample_id": "ID07", "sex": "female", "drb3": "*005:03/*016:01", "batch": "later"},
    {"sample_id": "ID08", "sex": "male", "drb3": "*010:01/*010:01", "batch": "2021"},
    {"sample_id": "ID09", "sex": "male", "drb3": "*012:01/*007:01", "batch": "2021"},
    {"sample_id": "ID10", "sex": "female", "drb3": "*010:01/*016:01", "batch": "later"},
    {"sample_id": "ID11", "sex": "female", "drb3": "*002:01/*011:01", "batch": "2021"},
    {"sample_id": "ID12", "sex": "female", "drb3": "*015:01/*012:01", "batch": "2021"},
    {"sample_id": "ID13", "sex": "female", "drb3": "*011:01/*005:03", "batch": "2021"},
)

ORGANISM = "Bos taurus"
BREED = "Japanese Black"
TISSUE = "blood"
GEO_LOC_NAME = "Japan"
# INSDC missing-value vocabulary. D-way rejects a bare "not collected" for
# collection_date and requires "missing: <reporting-level term>"; the authors
# chose "third party data" (D-way validation, 2026-10-09) rather than inventing
# a year or site.
COLLECTION_DATE = "missing: third party data"
DEV_STAGE = "not collected"
INSTRUMENT_MODEL = "PromethION"
PLATFORM = "OXFORD_NANOPORE"
LIBRARY_STRATEGY = "WGS"
LIBRARY_SOURCE = "GENOMIC"
LIBRARY_SELECTION = "RANDOM"
LIBRARY_LAYOUT = "SINGLE"

BIOSAMPLE_FIELDS = (
    "sample_name",
    "organism",
    "isolate",
    "breed",
    "sex",
    "tissue",
    "geo_loc_name",
    "collection_date",
    "dev_stage",
    "genotype",
    "description",
)

EXPERIMENT_FIELDS = (
    "sample_name",
    "title",
    "library_name",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "platform",
    "instrument_model",
    "design_description",
)

RUN_FIELDS = (
    "sample_name",
    "file_name",
    "file_type",
    "md5_checksum",
)

REGISTER_FIELDS = (
    "sample_id",
    "planned_fastq",
    "biosample_accession",
    "experiment_accession",
    "run_accession",
    "bioproject_accession",
    "md5_checksum",
)

BIOPROJECT_FIELDS = ("field", "value")

BIOPROJECT_ROWS: tuple[tuple[str, str], ...] = (
    ("title", "Long-read whole-genome sequencing of the BoLA region in Japanese Black cattle"),
    ("organism", ORGANISM),
    ("project_data_type", "Genome Sequencing"),
    ("relevance", "Agricultural"),
    (
        "description",
        "Untargeted long-read whole-genome sequencing of 13 Japanese Black (Wagyu) cattle "
        "using Oxford Nanopore PromethION sequencers, generated to resolve the structure, "
        "copy-number variation, and haplotype organisation of the bovine major "
        "histocompatibility complex (BoLA) region on chromosome 23, including the classical "
        "class I and class II (DRB3/DQ) loci and their linkage. The cohort was ascertained "
        "to include animals carrying the BoLA-DRB3*016:01 allele. Reads support a manuscript "
        "reporting corrected BoLA gene coordinates and structural variation on the "
        "ARS-UCD2.0 reference assembly. Public sample identifiers are manuscript IDs "
        "ID01–ID13; geographic origin is reported at the country level (Japan) only, "
        "without farm, prefecture, or owner identifying information.",
    ),
    ("release", "Hold until publication of the linked DRA records"),
)


def planned_fastq_name(sample_id: str) -> str:
    """Public FASTQ filename for DRA. Alphanumeric / underscore / hyphen / dot only."""
    return f"{sample_id}.ont.fastq.gz"


def _library_protocol(batch: str) -> str:
    if batch == "later":
        return (
            "Oxford Nanopore PromethION pass reads; reads shorter than 500 bp discarded "
            "to match the 2021 cohort length floor."
        )
    return (
        "1D Ligation Sequencing Kit, R9 flow cells, Guppy 5.0.12; NanoFilt retained "
        "reads >=500 bp after a 75 bp head trim."
    )


def biosample_rows() -> list[dict[str, str]]:
    rows = []
    for animal in COHORT:
        sample_id = animal["sample_id"]
        rows.append(
            {
                "sample_name": sample_id,
                "organism": ORGANISM,
                "isolate": sample_id,
                "breed": BREED,
                "sex": animal["sex"],
                "tissue": TISSUE,
                "geo_loc_name": GEO_LOC_NAME,
                "collection_date": COLLECTION_DATE,
                "dev_stage": DEV_STAGE,
                "genotype": f"BoLA-DRB3 PCR-SBT {animal['drb3']}",
                "description": (
                    "Japanese Black cattle blood DNA; untargeted Oxford Nanopore "
                    "whole-genome sequencing. Cohort preferentially ascertained for "
                    "BoLA-DRB3*016:01. Public identifier matches the manuscript sample ID."
                ),
            }
        )
    return rows


def experiment_rows() -> list[dict[str, str]]:
    rows = []
    for animal in COHORT:
        sample_id = animal["sample_id"]
        rows.append(
            {
                "sample_name": sample_id,
                "title": f"PromethION WGS of Japanese Black cattle {sample_id}",
                "library_name": f"{sample_id}_ont_wgs",
                "library_strategy": LIBRARY_STRATEGY,
                "library_source": LIBRARY_SOURCE,
                "library_selection": LIBRARY_SELECTION,
                "library_layout": LIBRARY_LAYOUT,
                "platform": PLATFORM,
                "instrument_model": INSTRUMENT_MODEL,
                "design_description": _library_protocol(animal["batch"]),
            }
        )
    return rows


def run_rows() -> list[dict[str, str]]:
    rows = []
    for animal in COHORT:
        sample_id = animal["sample_id"]
        rows.append(
            {
                "sample_name": sample_id,
                "file_name": planned_fastq_name(sample_id),
                "file_type": "fastq",
                "md5_checksum": "",
            }
        )
    return rows


def register_rows() -> list[dict[str, str]]:
    rows = []
    for animal in COHORT:
        sample_id = animal["sample_id"]
        rows.append(
            {
                "sample_id": sample_id,
                "planned_fastq": planned_fastq_name(sample_id),
                "biosample_accession": "",
                "experiment_accession": "",
                "run_accession": "",
                "bioproject_accession": "",
                "md5_checksum": "",
            }
        )
    return rows


def bioproject_rows() -> list[dict[str, str]]:
    return [{"field": field, "value": value} for field, value in BIOPROJECT_ROWS]


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fieldnames})


def _register_has_accessions(path: Path) -> bool:
    """True when a register already records issued accessions or checksums.

    The register is the only tracked record of INSDC identifiers, so a re-run
    must never blank it back to the empty template.
    """
    if not path.is_file():
        return False
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return any(cell.strip() for row in rows[1:] for cell in row[2:])


def write_ddbj_metadata(output_dir: Path) -> dict[str, Path]:
    """Write the four public worksheets. Returns paths keyed by basename."""
    output_dir = Path(output_dir)
    written = {
        "biosample_model_organism.tsv": biosample_rows(),
        "dra_experiment.tsv": experiment_rows(),
        "dra_run.tsv": run_rows(),
        "accession_register.tsv": register_rows(),
        "bioproject.tsv": bioproject_rows(),
    }
    fields = {
        "biosample_model_organism.tsv": BIOSAMPLE_FIELDS,
        "dra_experiment.tsv": EXPERIMENT_FIELDS,
        "dra_run.tsv": RUN_FIELDS,
        "accession_register.tsv": REGISTER_FIELDS,
        "bioproject.tsv": BIOPROJECT_FIELDS,
    }
    paths = {}
    for name, rows in written.items():
        path = output_dir / name
        if name == "accession_register.tsv" and _register_has_accessions(path):
            paths[name] = path
            continue
        _write_tsv(path, fields[name], rows)
        paths[name] = path
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/codex_reanalysis/ddbj_metadata"),
        help="directory for public TSV worksheets",
    )
    args = parser.parse_args()
    paths = write_ddbj_metadata(args.output_dir)
    for name, path in paths.items():
        print(f"wrote {path} ({name})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
