from __future__ import annotations

from pathlib import Path

import pytest

from scripts.codex_reanalysis.build_evidence_matrix import (
    ACCEPTED,
    LIMITED,
    REJECTED,
    EvidenceMatrixError,
    assert_no_identity_leak,
    build_matrix,
    classify_claim,
    contains_real_sample_id,
    validate_evidence_paths,
    write_matrix_tsv,
)
from scripts.codex_reanalysis.compare_legacy_claims import (
    LegacyComparisonError,
    excerpt_found,
    extract_docx_text,
    normalize_whitespace,
    render_comparison,
    verify_excerpts,
)


def _minimal_docx(tmp_path: Path, paragraphs: list[str]) -> Path:
    """Build a syntactically minimal .docx whose body text is `paragraphs`."""
    import zipfile

    body = "".join(f"<w:p><w:r><w:t>{text}</w:t></w:r></w:p>" for text in paragraphs)
    document_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        f"<w:body>{body}</w:body></w:document>"
    )
    path = tmp_path / "legacy.docx"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("word/document.xml", document_xml)
    return path


def _claim(**overrides) -> dict:
    base = {
        "claim_id": "example_claim",
        "section": "Results",
        "manuscript_files": ("docs/RESULTS_EN.md",),
        "legacy_excerpt": None,
        "legacy_value": "old value",
        "new_evidence_paths": (),
        "new_value": "new value",
        "source_exists": True,
        "reproducible": True,
        "independent_support": False,
        "limitation": False,
        "note": "note",
        "action": "keep as-is",
    }
    base.update(overrides)
    return base


def test_claim_without_new_source_is_rejected():
    assert classify_claim(source_exists=False, reproducible=False) == REJECTED


def test_claim_with_two_methods_and_scope_limit_is_limited():
    assert classify_claim(source_exists=True, reproducible=True, independent_support=True, limitation=True) == LIMITED


def test_claim_with_source_and_no_limitation_is_accepted():
    assert classify_claim(source_exists=True, reproducible=True) == ACCEPTED


def test_claim_not_reproduced_is_rejected_even_with_source():
    assert classify_claim(source_exists=True, reproducible=False) == REJECTED


@pytest.mark.parametrize("token", ["DELIVERY_02", "DELIVERY_10", "DELIVERY_08", "DELIVERY_01"])
def test_contains_real_sample_id_detects_known_prefixes(token):
    assert contains_real_sample_id(f"see {token} for detail")


def test_contains_real_sample_id_false_for_anonymized_id():
    assert not contains_real_sample_id("ID08 shows the highest ratio")


def test_assert_no_identity_leak_raises_on_real_sample_id():
    rows = [_claim(new_value="DELIVERY_08 reached 72x depth")]
    with pytest.raises(EvidenceMatrixError, match="real sample_id"):
        assert_no_identity_leak(rows)


def test_assert_no_identity_leak_passes_for_anonymized_only_rows():
    rows = [_claim(new_value="ID08 reached 72x depth")]
    assert_no_identity_leak(rows) is None


def test_validate_evidence_paths_raises_on_missing_file(tmp_path):
    with pytest.raises(EvidenceMatrixError, match="do not exist"):
        validate_evidence_paths(["does/not/exist.tsv"], tmp_path)


def test_validate_evidence_paths_passes_for_existing_file(tmp_path):
    (tmp_path / "evidence.tsv").write_text("data", encoding="utf-8")
    validate_evidence_paths(["evidence.tsv"], tmp_path)


def test_build_matrix_rejects_duplicate_claim_id(tmp_path):
    claims = [_claim(claim_id="dup"), _claim(claim_id="dup")]
    with pytest.raises(EvidenceMatrixError, match="duplicate claim_id"):
        build_matrix(claims, tmp_path)


def test_build_matrix_rejects_identity_leak_even_when_paths_are_valid(tmp_path):
    (tmp_path / "evidence.tsv").write_text("data", encoding="utf-8")
    claims = [_claim(new_evidence_paths=("evidence.tsv",), note="from DELIVERY_02")]
    with pytest.raises(EvidenceMatrixError, match="real sample_id"):
        build_matrix(claims, tmp_path)


def test_build_matrix_attaches_classification(tmp_path):
    (tmp_path / "evidence.tsv").write_text("data", encoding="utf-8")
    claims = [_claim(new_evidence_paths=("evidence.tsv",))]
    rows = build_matrix(claims, tmp_path)
    assert rows[0]["classification"] == ACCEPTED


def test_build_matrix_skips_path_check_when_source_does_not_exist(tmp_path):
    claims = [_claim(source_exists=False, reproducible=False, new_evidence_paths=("missing.tsv",))]
    rows = build_matrix(claims, tmp_path)
    assert rows[0]["classification"] == REJECTED


def test_write_matrix_tsv_round_trips_fields(tmp_path):
    rows = build_matrix([_claim()], tmp_path)
    output = tmp_path / "matrix.tsv"
    write_matrix_tsv(output, rows)
    text = output.read_text(encoding="utf-8")
    assert "example_claim" in text
    assert "accepted" in text


def test_normalize_whitespace_collapses_newlines_and_tabs():
    assert normalize_whitespace("a\n\tb   c") == "a b c"


def test_extract_docx_text_reads_minimal_docx(tmp_path):
    docx = _minimal_docx(tmp_path, ["Hello world.", "Second paragraph."])
    text = extract_docx_text(docx)
    assert "Hello world." in text
    assert "Second paragraph." in text


def test_extract_docx_text_raises_on_non_docx(tmp_path):
    fake = tmp_path / "not_a_docx.docx"
    fake.write_text("plain text, not a zip", encoding="utf-8")
    with pytest.raises(LegacyComparisonError):
        extract_docx_text(fake)


def test_excerpt_found_normalizes_whitespace():
    docx_text = normalize_whitespace("The   value   was  0.74   in this cohort.")
    assert excerpt_found(docx_text, "the value\nwas 0.74 in this cohort.".replace("the", "The"))


def test_excerpt_found_false_when_absent():
    docx_text = normalize_whitespace("Nothing relevant here.")
    assert not excerpt_found(docx_text, "a phrase that was never written")


def test_verify_excerpts_flags_missing_excerpt():
    docx_text = normalize_whitespace("Only this sentence exists.")
    rows = [{"claim_id": "a", "legacy_excerpt": "a phrase not in the docx"}]
    assert verify_excerpts(rows, docx_text) == ["a"]


def test_verify_excerpts_skips_rows_without_excerpt():
    docx_text = normalize_whitespace("Only this sentence exists.")
    rows = [{"claim_id": "a", "legacy_excerpt": ""}]
    assert verify_excerpts(rows, docx_text) == []


def test_verify_excerpts_passes_when_excerpt_present():
    docx_text = normalize_whitespace("The squared binary correlation was r squared 0.74.")
    rows = [{"claim_id": "a", "legacy_excerpt": "The squared binary correlation was r squared 0.74."}]
    assert verify_excerpts(rows, docx_text) == []


def test_render_comparison_includes_every_claim():
    rows = [
        {
            "claim_id": "a",
            "classification": "accepted",
            "legacy_value": "old",
            "new_value": "new",
            "action": "keep as-is",
        }
    ]
    rendered = render_comparison(rows)
    assert "a" in rendered
    assert "accepted" in rendered
