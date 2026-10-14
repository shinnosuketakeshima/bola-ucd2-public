import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.codex_reanalysis.stream_target_mapping import (
    MappingError,
    build_mapping_argv,
    build_mapping_command,
    build_target_filter_argv,
    build_target_filter_command,
    resolve_audited_sample,
    build_redacted_mapping_command,
    build_redacted_filter_commands,
    create_reference_seal,
    publish_staged_artifacts,
    parse_flagstat_total,
    invalidate_status_before_cleanup,
    reference_resource_sha256,
    summarize_sam_alignments,
    validate_mapping_summary,
    validate_mapping_config,
    validate_reference_seal,
    validate_stage_gates,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_mapping_command_uses_map_ont_and_keeps_secondary():
    command = build_mapping_command(
        "sample.fastq.gz",
        "ARS-UCD2.0.fa",
        threads=8,
        sample_id="ID01",
        deliverable_id="DELIVERY_01",
    )
    assert "-ax map-ont" in command
    assert " -L " in command
    assert "--secondary=yes" in command
    assert "-N 20" in command
    assert "-t 8" in command
    assert "ID:ID01" in command
    assert "PU:DELIVERY_01" in command


def test_mapping_argv_is_a_shell_free_argument_vector():
    argv = build_mapping_argv(
        "reads with spaces.fastq.gz",
        "reference with spaces.fa",
        threads=2,
        sample_id="ID01",
        deliverable_id="DELIVERY_01",
    )
    assert argv[0] == "minimap2"
    assert "-L" in argv
    assert argv[-2:] == ["reference with spaces.fa", "reads with spaces.fastq.gz"]
    assert "|" not in argv


def test_target_filter_command_never_writes_whole_genome_sam_or_bam():
    command = build_target_filter_command(
        "targets.bed",
        "out.targets.bam.partial",
        temp_prefix="tmp/ID01/sort",
        threads=4,
    )
    assert "samtools view" in command
    assert "-L targets.bed" in command
    assert "samtools sort" in command
    assert ".sam" not in command
    assert "whole" not in command


def test_target_filter_argv_retains_secondary_and_supplementary_alignments():
    view_argv, sort_argv = build_target_filter_argv(
        "targets.bed",
        "out.targets.bam.partial",
        temp_prefix="tmp/ID01/sort",
        threads=4,
    )
    assert view_argv == ["samtools", "view", "-h", "-b", "-L", "targets.bed", "-"]
    assert "-F" not in view_argv
    assert sort_argv == [
        "samtools",
        "sort",
        "-@",
        "4",
        "-T",
        "tmp/ID01/sort",
        "-o",
        "out.targets.bam.partial",
        "-",
    ]


def test_stage_gate_requires_passing_audit_and_reference_manifest(tmp_path):
    audit = tmp_path / "audit_status.json"
    input_manifest = tmp_path / "input_manifest.tsv"
    sample_registry = tmp_path / "sample_registry.tsv"
    reference = tmp_path / "reference_status.json"
    manifest = tmp_path / "reference_manifest.tsv"
    input_manifest.write_text("sample_id\tsha256\nID01\t" + "a" * 64 + "\n", encoding="utf-8")
    sample_registry.write_text("sample_id\nID01\n", encoding="utf-8")
    audit.write_text(
        json.dumps(
            {
                "status": "pass",
                "analytical_sample_count": 13,
                "input_manifest_sha256": __import__("hashlib").sha256(
                    input_manifest.read_bytes()
                ).hexdigest(),
                "sample_registry_sha256": __import__("hashlib").sha256(
                    sample_registry.read_bytes()
                ).hexdigest(),
            }
        ),
        encoding="utf-8",
    )
    manifest.write_text("resource\tsha256\nreference.fa\tabc\n", encoding="utf-8")
    reference.write_text(
        json.dumps(
            {
                "status": "passed",
                "target_count": 57,
                "reference_manifest_sha256": __import__("hashlib").sha256(
                    manifest.read_bytes()
                ).hexdigest(),
            }
        ),
        encoding="utf-8",
    )

    validate_stage_gates(audit, input_manifest, sample_registry, reference, manifest)

    input_manifest.write_text(input_manifest.read_text() + "ID02\t" + "b" * 64 + "\n")
    with pytest.raises(MappingError, match="audit status is not bound"):
        validate_stage_gates(audit, input_manifest, sample_registry, reference, manifest)

    audit.write_text(json.dumps({"status": "failed"}), encoding="utf-8")
    with pytest.raises(MappingError, match="input audit is not passing"):
        validate_stage_gates(audit, input_manifest, sample_registry, reference, manifest)


def test_reference_fasta_hash_is_reused_from_the_verified_manifest(tmp_path):
    manifest = tmp_path / "reference_manifest.tsv"
    manifest.write_text(
        "resource\tsha256\nARS-UCD2.0.fa\tverified-digest\n",
        encoding="utf-8",
    )
    assert reference_resource_sha256(manifest, "ARS-UCD2.0.fa") == "verified-digest"


def test_mapq_summary_counts_primary_secondary_and_supplementary_records():
    rows = summarize_sam_alignments(
        [
            "read1\t0\tNC_1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
            "read1\t256\tNC_2\t1\t10\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
            "read2\t2048\tNC_1\t20\t0\t5S5M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF",
        ]
    )
    assert rows == {
        "alignment_count": 3,
        "primary_count": 1,
        "secondary_count": 1,
        "supplementary_count": 1,
        "mapq_0": 1,
        "mapq_1_9": 0,
        "mapq_10_19": 1,
        "mapq_20_59": 0,
        "mapq_60_plus": 1,
    }


def test_zero_target_alignments_are_a_hard_failure():
    with pytest.raises(MappingError, match="zero target alignments"):
        validate_mapping_summary({"alignment_count": 0, "primary_count": 0})


def test_mapping_script_supports_direct_cli():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "stream_target_mapping.py"),
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_audited_sample_resolution_stays_inside_external_delivery_root(tmp_path):
    deliver_root = tmp_path / "delivery"
    fastq = deliver_root / "pass" / "ID01.fastq.gz"
    fastq.parent.mkdir(parents=True)
    fastq.write_bytes(b"synthetic")
    registry = tmp_path / "sample_registry.tsv"
    registry.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tanalysis_status\n"
        "DELIVERY_01\tID01\tpass/ID01.fastq.gz\tanalytical\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "input_manifest.tsv"
    manifest.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tbyte_size\tsha256\n"
        "DELIVERY_01\tID01\tpass/ID01.fastq.gz\t9\t"
        "2b68c8fb159f5cceba9f666b3f39a4d5c53f5ad6f2ef92d47b1f91e20d2a3b3c\n",
        encoding="utf-8",
    )

    resolved = resolve_audited_sample("ID01", registry, manifest, deliver_root)

    digest = "2b68c8fb159f5cceba9f666b3f39a4d5c53f5ad6f2ef92d47b1f91e20d2a3b3c"
    assert resolved == {
        "deliverable_id": "DELIVERY_01",
        "fastq": fastq.resolve(),
        "input_sha256": digest,
        # Without a normalisation register the delivered file is both the source
        # and the mapped input, and that must be stated rather than implied.
        "source_sha256": digest,
        "input_normalization": "delivered",
        # Mapping re-checks the size, not the hash, so the record must say so
        # rather than let provenance imply a fresh verification.
        "input_hash_verified_at_mapping": False,
        "input_hash_verification": "not_verified",
    }


def test_an_input_seal_makes_the_delivered_hash_a_verified_fact(tmp_path):
    """The seal was produced by hashing the delivery on this host, so mapping
    can record the hash as verified here rather than carried over."""
    deliver_root, registry, manifest, digest = _delivery_fixture(tmp_path)
    seal = {"status": "verified", "samples": {"ID01": {"sha256": digest, "byte_size": 9}}}

    resolved = resolve_audited_sample("ID01", registry, manifest, deliver_root, input_seal=seal)

    assert resolved["input_hash_verified_at_mapping"] is True
    assert resolved["input_hash_verification"] == "input_seal"


def test_a_seal_disagreeing_with_the_audit_stops_mapping(tmp_path):
    deliver_root, registry, manifest, _ = _delivery_fixture(tmp_path)
    seal = {"status": "verified", "samples": {"ID01": {"sha256": "0" * 64, "byte_size": 9}}}

    with pytest.raises(MappingError, match="seal"):
        resolve_audited_sample("ID01", registry, manifest, deliver_root, input_seal=seal)


def test_a_seal_that_does_not_cover_the_sample_stops_mapping(tmp_path):
    """Silently falling back to unverified would make the record dishonest."""
    deliver_root, registry, manifest, _ = _delivery_fixture(tmp_path)
    seal = {"status": "verified", "samples": {"ID99": {"sha256": "0" * 64, "byte_size": 9}}}

    with pytest.raises(MappingError, match="seal"):
        resolve_audited_sample("ID01", registry, manifest, deliver_root, input_seal=seal)


def _delivery_fixture(tmp_path):
    deliver_root = tmp_path / "delivery"
    fastq = deliver_root / "pass" / "ID01.fastq.gz"
    fastq.parent.mkdir(parents=True)
    fastq.write_bytes(b"delivered")
    digest = __import__("hashlib").sha256(b"delivered").hexdigest()
    registry = tmp_path / "sample_registry.tsv"
    registry.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tanalysis_status\n"
        "DELIVERY_14\tID01\tpass/ID01.fastq.gz\tanalytical\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "input_manifest.tsv"
    manifest.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tbyte_size\tsha256\n"
        f"DELIVERY_14\tID01\tpass/ID01.fastq.gz\t9\t{digest}\n",
        encoding="utf-8",
    )
    return deliver_root, registry, manifest, digest


def _write_normalized(tmp_path, source_digest, payload=b"filtered", sha_override=None):
    normalized_dir = tmp_path / "00b_normalized"
    normalized_dir.mkdir()
    derived = normalized_dir / "ID01.min500.fq.gz"
    derived.write_bytes(payload)
    digest = sha_override or __import__("hashlib").sha256(payload).hexdigest()
    path = normalized_dir / "normalized_manifest.tsv"
    path.write_text(
        "deliverable_id\tsample_id\taction\tsource_sha256\tnormalized_relative_path\t"
        "normalized_sha256\tnormalized_byte_size\n"
        f"DELIVERY_14\tID01\tlength_filter\t{source_digest}\tID01.min500.fq.gz\t"
        f"{digest}\t{len(payload)}\n",
        encoding="utf-8",
    )
    return path, derived


def test_audited_sample_uses_a_registered_normalized_derivative(tmp_path):
    deliver_root, registry, manifest, digest = _delivery_fixture(tmp_path)
    normalized_manifest, derived = _write_normalized(tmp_path, digest)

    resolved = resolve_audited_sample(
        "ID01", registry, manifest, deliver_root, normalized_manifest_path=normalized_manifest
    )

    assert resolved["fastq"] == derived.resolve()
    assert resolved["input_normalization"] == "length_filtered"
    assert resolved["source_sha256"] == digest
    assert resolved["input_sha256"] != digest


def test_audited_sample_keeps_the_delivered_file_when_normalization_was_not_required(tmp_path):
    deliver_root, registry, manifest, digest = _delivery_fixture(tmp_path)
    normalized_dir = tmp_path / "00b_normalized"
    normalized_dir.mkdir()
    normalized_manifest = normalized_dir / "normalized_manifest.tsv"
    normalized_manifest.write_text(
        "deliverable_id\tsample_id\taction\tsource_sha256\tnormalized_relative_path\t"
        "normalized_sha256\tnormalized_byte_size\n"
        f"DELIVERY_14\tID01\tnot_required\t{digest}\t\t\t\n",
        encoding="utf-8",
    )

    resolved = resolve_audited_sample(
        "ID01", registry, manifest, deliver_root, normalized_manifest_path=normalized_manifest
    )

    assert resolved["fastq"] == (deliver_root / "pass" / "ID01.fastq.gz").resolve()
    assert resolved["input_normalization"] == "delivered"
    assert resolved["input_sha256"] == digest


def test_audited_sample_resolves_a_normalized_register_reached_through_a_symlink(tmp_path):
    """The remote host links the 240 GB register from a read-only mount.

    Resolution must follow the link rather than reject it, or the documented
    remote layout would fail only once the drive is on the other machine.
    """
    deliver_root, registry, manifest, digest = _delivery_fixture(tmp_path)
    real_manifest, derived = _write_normalized(tmp_path, digest)
    linked_project = tmp_path / "project"
    linked_project.mkdir()
    link = linked_project / "00b_normalized"
    link.symlink_to(real_manifest.parent, target_is_directory=True)

    resolved = resolve_audited_sample(
        "ID01",
        registry,
        manifest,
        deliver_root,
        normalized_manifest_path=link / "normalized_manifest.tsv",
    )

    assert resolved["fastq"] == derived.resolve()
    assert resolved["input_normalization"] == "length_filtered"


def test_audited_sample_rejects_a_normalized_derivative_whose_hash_changed(tmp_path):
    deliver_root, registry, manifest, digest = _delivery_fixture(tmp_path)
    normalized_manifest, _ = _write_normalized(tmp_path, digest, sha_override="0" * 64)

    with pytest.raises(MappingError, match="normalized"):
        resolve_audited_sample(
            "ID01", registry, manifest, deliver_root, normalized_manifest_path=normalized_manifest
        )


def test_audited_sample_rejects_a_derivative_built_from_a_different_source(tmp_path):
    """The derivative must be bound to the FASTQ the audit actually hashed."""
    deliver_root, registry, manifest, _ = _delivery_fixture(tmp_path)
    normalized_manifest, _ = _write_normalized(tmp_path, "f" * 64)

    with pytest.raises(MappingError, match="source"):
        resolve_audited_sample(
            "ID01", registry, manifest, deliver_root, normalized_manifest_path=normalized_manifest
        )


def test_audited_sample_rejects_forged_hash_and_size(tmp_path):
    deliver_root = tmp_path / "delivery"
    fastq = deliver_root / "ID01.fastq.gz"
    deliver_root.mkdir()
    fastq.write_bytes(b"synthetic")
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tanalysis_status\n"
        "D01\tID01\tID01.fastq.gz\tanalytical\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "deliverable_id\tsample_id\tinput_relative_path\tbyte_size\tsha256\n"
        "D01\tID01\tID01.fastq.gz\t8\tnot-a-sha\n",
        encoding="utf-8",
    )
    with pytest.raises(MappingError, match="byte size|SHA-256"):
        resolve_audited_sample("ID01", registry, manifest, deliver_root)


def test_production_cli_has_no_unaudited_direct_fastq_mode():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "stream_target_mapping.py"),
            "map",
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "--sample-registry" in result.stdout
    assert "--input-manifest" in result.stdout
    assert "--deliver-root" in result.stdout
    assert "--normalized-manifest" in result.stdout
    assert "--fastq" not in result.stdout
    assert "--deliverable-id" not in result.stdout


def test_reference_seal_detects_changed_reference(tmp_path):
    reference = tmp_path / "reference.fa"
    targets = tmp_path / "targets.bed"
    manifest = tmp_path / "reference_manifest.tsv"
    config = tmp_path / "config.yaml"
    seal = tmp_path / "seal.json"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")
    targets.write_text("chr1\t0\t4\ttarget\n", encoding="utf-8")
    config.write_text("mapping: {}\n", encoding="utf-8")
    manifest.write_text(
        "resource\tsha256\n"
        f"reference.fa\t{__import__('hashlib').sha256(reference.read_bytes()).hexdigest()}\n"
        f"targets.bed\t{__import__('hashlib').sha256(targets.read_bytes()).hexdigest()}\n",
        encoding="utf-8",
    )
    create_reference_seal(reference, targets, manifest, config, seal)
    validate_reference_seal(seal, reference, targets, manifest, config)
    reference.write_text(">chr1\nTTTT\n", encoding="utf-8")
    with pytest.raises(MappingError, match="reference seal"):
        validate_reference_seal(seal, reference, targets, manifest, config)


def test_reference_seal_detects_same_size_edit_even_with_unchanged_stat(tmp_path):
    """A same-size content edit must be caught even if size/mtime happen to match.

    exFAT mtime resolution is coarse (and even finer filesystems are not
    immune to two writes landing in the same tick), so the seal must not
    trust a size/mtime match alone -- it stores a SHA-256 specifically so a
    same-tick edit is still detected.
    """
    reference = tmp_path / "reference.fa"
    targets = tmp_path / "targets.bed"
    manifest = tmp_path / "reference_manifest.tsv"
    config = tmp_path / "config.yaml"
    seal = tmp_path / "seal.json"
    reference.write_text(">chr1\nAAAA\n", encoding="utf-8")
    targets.write_text("chr1\t0\t4\ttarget\n", encoding="utf-8")
    config.write_text("mapping: {}\n", encoding="utf-8")
    manifest.write_text(
        "resource\tsha256\n"
        f"reference.fa\t{__import__('hashlib').sha256(reference.read_bytes()).hexdigest()}\n"
        f"targets.bed\t{__import__('hashlib').sha256(targets.read_bytes()).hexdigest()}\n",
        encoding="utf-8",
    )
    create_reference_seal(reference, targets, manifest, config, seal)
    stat_before = reference.stat()

    # Same byte count as "AAAA", and the stat is pinned to the sealed values
    # via os.utime so this reproduces a same-tick collision deterministically
    # instead of depending on filesystem timestamp resolution.
    reference.write_text(">chr1\nTTTT\n", encoding="utf-8")
    os.utime(reference, ns=(stat_before.st_atime_ns, stat_before.st_mtime_ns))
    assert reference.stat().st_size == stat_before.st_size
    assert reference.stat().st_mtime_ns == stat_before.st_mtime_ns

    with pytest.raises(MappingError, match="reference seal"):
        validate_reference_seal(seal, reference, targets, manifest, config)


def test_logged_mapping_command_redacts_external_fastq_path():
    command = build_redacted_mapping_command(
        "private/pass/ID01.fastq.gz", "reference.fa", threads=8,
        sample_id="ID01", deliverable_id="D01",
    )
    assert "private/pass" not in command
    assert "<audited-fastq:ID01>" in command

    view, sort = build_redacted_filter_commands(threads=8)
    assert "<validated-target-bed>" in view
    assert "<sample-temp-prefix>" in sort
    assert "/mnt/" not in view + sort


def test_mapping_config_enforces_preset_retention_and_versions(tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text(
        "tools:\n  minimap2_version: 2.31-r1302\n  samtools_version: 1.23.1\n"
        "mapping:\n  threads: 8\n  preset: map-ont\n"
        "  retain_secondary: true\n  retain_supplementary: true\n  secondary_max: 20\n",
        encoding="utf-8",
    )
    validate_mapping_config(config, 8, "2.31-r1302", "samtools 1.23.1")
    config.write_text(config.read_text().replace("map-ont", "map-pb"), encoding="utf-8")
    with pytest.raises(MappingError, match="map-ont"):
        validate_mapping_config(config, 8, "2.31-r1302", "samtools 1.23.1")
    config.write_text(config.read_text().replace("map-pb", "map-ont"), encoding="utf-8")
    with pytest.raises(MappingError, match="samtools version"):
        validate_mapping_config(config, 8, "2.31-r1302", "samtools 1.23.10")
    config.write_text(config.read_text().replace("secondary_max: 20", "secondary_max: 5"), encoding="utf-8")
    with pytest.raises(MappingError, match="secondary"):
        validate_mapping_config(config, 8, "2.31-r1302", "samtools 1.23.1", secondary_max=20)


def test_flagstat_total_is_parsed_for_count_reconciliation(tmp_path):
    flagstat = tmp_path / "flagstat.txt"
    flagstat.write_text(
        "3 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "3 + 0 primary\n",
        encoding="utf-8",
    )
    assert parse_flagstat_total(flagstat) == 3


def test_artifacts_are_not_published_until_every_staged_file_exists(tmp_path):
    staged_bam = tmp_path / "stage.bam"
    staged_bai = tmp_path / "stage.bam.bai"
    final_bam = tmp_path / "final.bam"
    final_bai = tmp_path / "final.bam.bai"
    staged_bam.write_bytes(b"bam")
    with pytest.raises(MappingError, match="staged mapping artifact"):
        publish_staged_artifacts({staged_bam: final_bam, staged_bai: final_bai})
    assert not final_bam.exists()
    assert not final_bai.exists()


def test_the_staging_directory_is_emptied_not_selectively_cleaned(tmp_path):
    """An artefact left by an earlier run can satisfy a later existence check.
    Enumerating files to delete leaves behind anything not on the list, which is
    how a stale idxstats.txt survived a re-map; clearing the directory cannot."""
    status = tmp_path / "sample.mapping_status.json"
    staging = tmp_path / "publish"
    staging.mkdir()
    (staging / "sample.idxstats.txt").write_text("stale\n", encoding="utf-8")
    (staging / "provenance").mkdir()

    invalidate_status_before_cleanup(status, "sample", [], staging_dir=staging)

    assert staging.is_dir()
    assert list(staging.iterdir()) == []


def test_the_staging_directory_is_recreated_even_when_absent(tmp_path):
    status = tmp_path / "sample.mapping_status.json"
    staging = tmp_path / "publish"

    invalidate_status_before_cleanup(status, "sample", [], staging_dir=staging)

    assert staging.is_dir()


def test_cleanup_failure_cannot_leave_a_passed_status(tmp_path):
    status = tmp_path / "sample.mapping_status.json"
    status.write_text('{"status":"passed"}\n', encoding="utf-8")
    invalid_cleanup_target = tmp_path / "directory"
    invalid_cleanup_target.mkdir()
    with pytest.raises(OSError):
        invalidate_status_before_cleanup(status, "sample", [invalid_cleanup_target])
    assert json.loads(status.read_text(encoding="utf-8"))["status"] == "running"
