from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

from scripts.codex_reanalysis.stage_dra_fastq import (
    EXPECTED_IDS,
    StageError,
    is_delivery_path,
    public_read_id,
    stage_dra_fastq,
)
from scripts.codex_reanalysis.write_ddbj_metadata import planned_fastq_name


def _write_gzip_fastq(path: Path, header: str, sequence: str = "ACGT") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    quality = "I" * len(sequence)
    with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
        handle.write(f"{header}\n{sequence}\n+\n{quality}\n")


def _write_map(path: Path, sources: dict[str, Path]) -> None:
    lines = ["manuscript_id\tfastq"]
    for sample_id in EXPECTED_IDS:
        lines.append(f"{sample_id}\t{sources[sample_id]}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_public_read_id_keeps_uuid_and_replaces_vendor_token() -> None:
    assert public_read_id("@abcd-uuid runid=x sampleid=DELIVERY_02", "ID01", 1) == "abcd-uuid"
    assert public_read_id("@DELIVERY_02 runid=x", "ID01", 7) == "ID01.read7"


def test_is_delivery_path_detects_deliver_add(tmp_path: Path) -> None:
    nested = tmp_path / "Deliver_ADD" / "Filtered" / "x.fastq.gz"
    nested.parent.mkdir(parents=True)
    nested.write_bytes(b"x")
    assert is_delivery_path(nested)


def test_stage_refuses_delivery_root(tmp_path: Path) -> None:
    delivery = tmp_path / "Deliver_ADD" / "x.fastq.gz"
    _write_gzip_fastq(delivery, "@r1")
    sources = {sample_id: tmp_path / f"{sample_id}.fastq.gz" for sample_id in EXPECTED_IDS}
    for sample_id, path in sources.items():
        _write_gzip_fastq(path, "@ok")
    sources["ID01"] = delivery
    fastq_map = tmp_path / "map.tsv"
    _write_map(fastq_map, sources)
    with pytest.raises(StageError, match="delivery root"):
        stage_dra_fastq(fastq_map, tmp_path / "stage", sanitize_headers=False, header_scan_reads=10)


def test_stage_allows_delivery_root_when_sanitizing(tmp_path: Path) -> None:
    delivery = tmp_path / "Deliver_ADD" / "x.fastq.gz"
    _write_gzip_fastq(delivery, "@uuid-1 runid=abc sampleid=DELIVERY_02")
    sources = {sample_id: tmp_path / f"{sample_id}.fastq.gz" for sample_id in EXPECTED_IDS}
    for sample_id, path in sources.items():
        _write_gzip_fastq(path, "@ok")
    sources["ID01"] = delivery
    fastq_map = tmp_path / "map.tsv"
    _write_map(fastq_map, sources)
    before = delivery.read_bytes()
    rows = stage_dra_fastq(
        fastq_map, tmp_path / "stage", sanitize_headers=True, header_scan_reads=10
    )
    assert len(rows) == 13
    assert delivery.read_bytes() == before


def test_stage_requires_sanitize_when_headers_leak(tmp_path: Path) -> None:
    sources = {}
    for sample_id in EXPECTED_IDS:
        path = tmp_path / "in" / f"{sample_id}.fastq.gz"
        _write_gzip_fastq(path, "@r1 sampleid=DELIVERY_02")
        sources[sample_id] = path
    fastq_map = tmp_path / "map.tsv"
    _write_map(fastq_map, sources)
    with pytest.raises(StageError, match="sanitize-headers"):
        stage_dra_fastq(fastq_map, tmp_path / "stage", sanitize_headers=False, header_scan_reads=10)


def test_sanitize_with_external_compressor_matches_python_gzip(tmp_path: Path) -> None:
    pigz = shutil.which("pigz")
    if pigz is None:
        pytest.skip("pigz is not installed")
    sources = {}
    for sample_id in EXPECTED_IDS:
        path = tmp_path / "in" / f"{sample_id}.fastq.gz"
        _write_gzip_fastq(path, "@uuid-1 runid=abc sampleid=DELIVERY_10")
        sources[sample_id] = path
    fastq_map = tmp_path / "map.tsv"
    _write_map(fastq_map, sources)
    stage_dir = tmp_path / "stage"
    rows = stage_dra_fastq(
        fastq_map,
        stage_dir,
        sanitize_headers=True,
        header_scan_reads=10,
        compressor=pigz,
    )
    assert len(rows) == 13
    staged = stage_dir / planned_fastq_name("ID01")
    with gzip.open(staged, "rt", encoding="ascii") as handle:
        header = handle.readline().rstrip("\n")
        sequence = handle.readline().rstrip("\n")
    assert header == "@uuid-1"
    assert sequence == "ACGT"
    assert rows[0]["md5_checksum"]


def test_sanitize_strips_vendor_tags_and_fills_md5(tmp_path: Path) -> None:
    sources = {}
    for sample_id in EXPECTED_IDS:
        path = tmp_path / "in" / f"{sample_id}.fastq.gz"
        _write_gzip_fastq(path, "@uuid-1 runid=abc sampleid=DELIVERY_10")
        sources[sample_id] = path
    fastq_map = tmp_path / "map.tsv"
    _write_map(fastq_map, sources)
    stage_dir = tmp_path / "stage"
    rows = stage_dra_fastq(fastq_map, stage_dir, sanitize_headers=True, header_scan_reads=10)
    assert len(rows) == 13
    staged = stage_dir / planned_fastq_name("ID01")
    with gzip.open(staged, "rt", encoding="ascii") as handle:
        header = handle.readline().rstrip("\n")
        sequence = handle.readline().rstrip("\n")
    assert header == "@uuid-1"
    assert sequence == "ACGT"
    filled = (stage_dir / "dra_run.filled.tsv").read_text(encoding="utf-8")
    assert "DELIVERY_" not in filled
    assert "DELIVERY_" not in filled
    assert rows[0]["md5_checksum"]
    assert "ID01.ont.fastq.gz" in (stage_dir / "sftp_commands.txt").read_text(encoding="utf-8")
    assert "ftp-private.ddbj.nig.ac.jp" in (stage_dir / "sftp_commands.txt").read_text(
        encoding="utf-8"
    )
