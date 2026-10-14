from __future__ import annotations

from pathlib import Path

from scripts.codex_reanalysis.build_evidence_matrix import contains_real_sample_id
from scripts.codex_reanalysis.write_ddbj_metadata import (
    COHORT,
    COLLECTION_DATE,
    GEO_LOC_NAME,
    INSTRUMENT_MODEL,
    write_ddbj_metadata,
)

REPO = Path(__file__).resolve().parents[2]
TRACKED_DIR = REPO / "docs" / "codex_reanalysis" / "ddbj_metadata"


def test_cohort_has_thirteen_manuscript_ids() -> None:
    ids = [row["sample_id"] for row in COHORT]
    assert ids == [f"ID{i:02d}" for i in range(1, 14)]


def test_id13_sex_is_author_confirmed() -> None:
    assert next(row["sex"] for row in COHORT if row["sample_id"] == "ID13") == "female"


def test_writer_emits_japan_only_geo_and_promethion(tmp_path: Path) -> None:
    paths = write_ddbj_metadata(tmp_path)
    biosample = paths["biosample_model_organism.tsv"].read_text(encoding="utf-8")
    experiment = paths["dra_experiment.tsv"].read_text(encoding="utf-8")
    assert biosample.count(f"\t{GEO_LOC_NAME}\t") == 13
    assert "Site A" not in biosample
    assert INSTRUMENT_MODEL in experiment
    assert "ID01.ont.fastq.gz" in paths["dra_run.tsv"].read_text(encoding="utf-8")


def test_collection_date_uses_dway_missing_term_not_invented_year(tmp_path: Path) -> None:
    assert COLLECTION_DATE == "missing: third party data"
    biosample = write_ddbj_metadata(tmp_path)["biosample_model_organism.tsv"].read_text(
        encoding="utf-8"
    )
    assert biosample.count("\tmissing: third party data\t") == 13
    assert "missing: not collected" not in biosample


def test_writer_does_not_mint_accessions(tmp_path: Path) -> None:
    paths = write_ddbj_metadata(tmp_path)
    register = paths["accession_register.tsv"].read_text(encoding="utf-8")
    bioproject = paths["bioproject.tsv"].read_text(encoding="utf-8")
    assert "PRJDB" not in register
    assert "SAMD" not in register
    assert "DRR" not in register
    assert "PRJDB" not in bioproject


def test_example_fastq_map_is_placeholders_only() -> None:
    text = (TRACKED_DIR / "fastq_map.example.tsv").read_text(encoding="utf-8")
    assert not contains_real_sample_id(text)
    assert "Deliver_ADD" not in text
    rows = [line for line in text.splitlines() if line.strip()]
    assert len(rows) == 14
    assert all(f"ID{i:02d}\t/path/on/analysis/host/" in text for i in range(1, 14))


def test_public_worksheets_have_no_identity_leak() -> None:
    for path in TRACKED_DIR.glob("*.tsv"):
        text = path.read_text(encoding="utf-8")
        assert not contains_real_sample_id(text), path
        assert "Deliver_ADD" not in text
        assert "ANONYMIZATION" not in text


def test_tracked_worksheets_match_writer(tmp_path: Path) -> None:
    generated = write_ddbj_metadata(tmp_path)
    for name, path in generated.items():
        tracked = TRACKED_DIR / name
        if name == "accession_register.tsv":
            # hand-maintained once D-way issues accessions; only the shape is fixed
            assert tracked.read_text(encoding="utf-8").splitlines()[0] == path.read_text(
                encoding="utf-8"
            ).splitlines()[0]
            continue
        assert tracked.read_text(encoding="utf-8") == path.read_text(encoding="utf-8")


def test_rewrite_preserves_issued_accessions(tmp_path: Path) -> None:
    register = write_ddbj_metadata(tmp_path)["accession_register.tsv"]
    filled = register.read_text(encoding="utf-8").replace("ID01.ont.fastq.gz\t", "ID01.ont.fastq.gz\tSAMD00000001", 1)
    register.write_text(filled, encoding="utf-8")
    write_ddbj_metadata(tmp_path)
    assert "SAMD00000001" in register.read_text(encoding="utf-8")


