import json
from hashlib import sha256

import pytest

from scripts.codex_reanalysis.input_seal import (
    InputSealError,
    build_input_seal,
    validate_input_seal,
    verified_sample,
)


def make_delivery(tmp_path, payloads=(("ID01", b"one"), ("ID02", b"two"))):
    root = tmp_path / "delivery"
    (root / "Filtered").mkdir(parents=True)
    rows = []
    for sample_id, payload in payloads:
        path = root / "Filtered" / f"{sample_id}.fq.gz"
        path.write_bytes(payload)
        rows.append(
            {
                "sample_id": sample_id,
                "input_relative_path": f"Filtered/{sample_id}.fq.gz",
                "byte_size": str(len(payload)),
                "sha256": sha256(payload).hexdigest(),
            }
        )
    manifest = tmp_path / "input_manifest.tsv"
    manifest.write_text(
        "sample_id\tinput_relative_path\tbyte_size\tsha256\n"
        + "".join(
            f"{r['sample_id']}\t{r['input_relative_path']}\t{r['byte_size']}\t{r['sha256']}\n"
            for r in rows
        ),
        encoding="utf-8",
    )
    return root, manifest


def test_seal_records_every_sample_it_verified(tmp_path):
    root, manifest = make_delivery(tmp_path)

    seal = build_input_seal(root, manifest)

    assert seal["status"] == "verified"
    assert seal["sample_count"] == 2
    assert set(seal["samples"]) == {"ID01", "ID02"}
    assert seal["samples"]["ID01"]["sha256"] == sha256(b"one").hexdigest()


def test_changed_content_is_rejected_not_recorded(tmp_path):
    root, manifest = make_delivery(tmp_path)
    (root / "Filtered" / "ID01.fq.gz").write_bytes(b"ONE")  # same length, different bytes

    with pytest.raises(InputSealError, match="content changed"):
        build_input_seal(root, manifest)


def test_changed_size_is_rejected_before_hashing(tmp_path):
    root, manifest = make_delivery(tmp_path)
    (root / "Filtered" / "ID01.fq.gz").write_bytes(b"one-longer")

    with pytest.raises(InputSealError, match="size changed"):
        build_input_seal(root, manifest)


def test_an_absent_delivered_file_is_rejected(tmp_path):
    root, manifest = make_delivery(tmp_path)
    (root / "Filtered" / "ID02.fq.gz").unlink()

    with pytest.raises(InputSealError, match="absent"):
        build_input_seal(root, manifest)


def test_a_valid_seal_passes_the_cheap_check(tmp_path):
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)

    validate_input_seal(seal, manifest)


def test_a_seal_built_against_a_different_manifest_is_rejected(tmp_path):
    """Re-running the audit changes the manifest, so the old seal no longer
    proves anything about what mapping is about to read."""
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)
    manifest.write_text(manifest.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(InputSealError, match="different audit manifest"):
        validate_input_seal(seal, manifest)


def test_a_seal_missing_a_sample_is_rejected(tmp_path):
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)
    del seal["samples"]["ID02"]

    with pytest.raises(InputSealError, match="missing samples"):
        validate_input_seal(seal, manifest)


def test_the_cheap_check_does_not_reread_the_delivered_files(tmp_path):
    """The host runs this before every long job; it must stay cheap."""
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)
    for sample in ("ID01", "ID02"):
        (root / "Filtered" / f"{sample}.fq.gz").unlink()

    validate_input_seal(seal, manifest)


def test_a_sample_outside_the_seal_cannot_be_claimed_as_verified(tmp_path):
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)

    with pytest.raises(InputSealError, match="does not cover"):
        verified_sample(seal, "ID99")


def test_a_verified_sample_returns_its_recorded_digest(tmp_path):
    root, manifest = make_delivery(tmp_path)
    seal = build_input_seal(root, manifest)

    assert verified_sample(seal, "ID01")["sha256"] == sha256(b"one").hexdigest()
