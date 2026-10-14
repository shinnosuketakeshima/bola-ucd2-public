import gzip
from hashlib import md5, sha256
from pathlib import Path
import subprocess
import sys

import yaml

import pytest
import scripts.codex_reanalysis.reference_registry as reference_registry_module

from scripts.codex_reanalysis.reference_registry import (
    ReferenceRegistryError,
    _decompress_gzip,
    _safe_resource_name,
    build_reference_registry,
    build_index_commands,
    compare_sequence_records,
    compare_fasta_contigs,
    parse_assembly_report_aliases,
    parse_official_checksums,
    sha256_file,
    summarize_ipd_alignments,
    validate_tool_version,
    verify_official_checksums,
    write_reference_manifest,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_sha256_file_streams_an_exact_digest(tmp_path):
    path = tmp_path / "reference.fa"
    path.write_bytes(b">chr23\nACGT\n")
    assert sha256_file(path) == sha256(b">chr23\nACGT\n").hexdigest()


def test_official_checksum_manifest_is_verified(tmp_path):
    fasta = tmp_path / "genomic.fna.gz"
    fasta.write_bytes(b"fixed-reference")
    manifest = tmp_path / "md5checksums.txt"
    manifest.write_text(
        f"{md5(b'fixed-reference').hexdigest()}  ./genomic.fna.gz\n",
        encoding="ascii",
    )

    expected = parse_official_checksums(manifest)
    verified = verify_official_checksums(tmp_path, expected, ["genomic.fna.gz"])

    assert verified["genomic.fna.gz"] == md5(b"fixed-reference").hexdigest()


def test_checksum_mismatch_is_a_hard_failure(tmp_path):
    (tmp_path / "genomic.gff.gz").write_bytes(b"changed")
    with pytest.raises(ReferenceRegistryError, match="checksum mismatch"):
        verify_official_checksums(
            tmp_path,
            {"genomic.gff.gz": md5(b"expected").hexdigest()},
            ["genomic.gff.gz"],
        )


def test_checksum_manifest_rejects_paths_outside_package(tmp_path):
    manifest = tmp_path / "md5checksums.txt"
    manifest.write_text("d41d8cd98f00b204e9800998ecf8427e  ../escape\n", encoding="ascii")
    with pytest.raises(ReferenceRegistryError, match="unsafe"):
        parse_official_checksums(manifest)


def test_index_commands_pin_fasta_and_never_use_shell_strings():
    commands = build_index_commands("samtools", "minimap2", "reference.fa")
    assert commands == [
        ["samtools", "faidx", "-o", "reference.fa.fai.partial", "reference.fa"],
        ["minimap2", "-d", "reference.fa.mmi.partial", "reference.fa"],
    ]


def test_existing_decompressed_resource_is_regenerated_from_verified_gzip(tmp_path):
    source = tmp_path / "reference.fa.gz"
    source.write_bytes(gzip.compress(b">ref\nACGT\n"))
    destination = tmp_path / "reference.fa"
    destination.write_bytes(b"stale-or-corrupt")

    _decompress_gzip(source, destination)

    assert destination.read_bytes() == b">ref\nACGT\n"


@pytest.mark.parametrize("name", [r"..\escape.gz", "../escape.gz", "C:escape.gz", "x?y.gz"])
def test_resource_name_rejects_cross_platform_path_and_url_syntax(name):
    with pytest.raises(ReferenceRegistryError, match="unsafe"):
        _safe_resource_name(name)


def test_different_autosomal_records_require_scientific_review():
    comparison = compare_sequence_records(
        contig="BTA23",
        old_accession="NC_037350.1",
        old_sequence="ACGT",
        new_accession="NC_037350.1",
        new_sequence="ACGA",
    )
    assert comparison.status == "decision_required"
    assert comparison.old_length == comparison.new_length == 4
    assert comparison.old_sha256 != comparison.new_sha256


def test_identical_autosomal_records_are_explicitly_accepted():
    comparison = compare_sequence_records(
        contig="BTA22",
        old_accession="NC_037349.1",
        old_sequence="AACCGG",
        new_accession="NC_037349.1",
        new_sequence="AACCGG",
    )
    assert comparison.status == "identical"


def test_reference_manifest_records_both_official_and_sha256_checksums(tmp_path):
    resource = tmp_path / "genomic.gff.gz"
    resource.write_bytes(b"annotation")
    path = write_reference_manifest(
        tmp_path / "reference_manifest.tsv",
        resources=[
            {
                "name": "genomic.gff.gz",
                "path": resource,
                "url": "https://ftp.ncbi.nlm.nih.gov/genomic.gff.gz",
                "official_md5": md5(b"annotation").hexdigest(),
            }
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="GCF_002263795.3-RS_2024_12",
        retrieval_date="2026-07-31",
    )
    rows = path.read_text(encoding="utf-8").splitlines()
    assert rows[0].split("\t") == [
        "assembly_accession",
        "annotation_release",
        "retrieval_date",
        "resource",
        "url",
        "size_bytes",
        "official_md5",
        "sha256",
    ]
    assert sha256(b"annotation").hexdigest() in rows[1]


def test_assembly_report_defines_one_refseq_accession_per_chromosome(tmp_path):
    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "1\tassembled-molecule\t1\tChromosome\tCM000000.1\t=\tNC_000001.1\tPrimary Assembly\t100\tchr1\n"
        "23\tassembled-molecule\t23\tChromosome\tCM000023.1\t=\tNC_037350.1\tPrimary Assembly\t200\tchr23\n",
        encoding="utf-8",
    )
    aliases = parse_assembly_report_aliases(report)
    assert aliases["23"] == ["NC_037350.1"]
    assert aliases["chr23"] == ["NC_037350.1"]


def test_assembly_report_rejects_ambiguous_chromosome_alias(tmp_path):
    report = tmp_path / "assembly_report.txt"
    header = "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    report.write_text(
        header
        + "23\tassembled-molecule\t23\tChromosome\tCM1\t=\tNC_A\tPrimary\t200\tchr23\n"
        + "23_alt\tassembled-molecule\t23\tChromosome\tCM2\t=\tNC_B\tPrimary\t200\tchr23\n",
        encoding="utf-8",
    )
    with pytest.raises(ReferenceRegistryError, match="ambiguous"):
        parse_assembly_report_aliases(report)


def test_assembly_report_skips_molecule_without_refseq_when_canonical_row_exists(tmp_path):
    report = tmp_path / "assembly_report.txt"
    header = "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
    report.write_text(
        header
        + "M\tassembled-molecule\tMT\tMitochondrion\tCM1\t<>\tna\tnon-nuclear\t100\tna\n"
        + "MT\tassembled-molecule\tMT\tMitochondrion\tna\t<>\tNC_MT\tnon-nuclear\t100\tchrM\n"
        + "23\tassembled-molecule\t23\tChromosome\tCM23\t=\tNC_23\tPrimary\t200\tchr23\n",
        encoding="utf-8",
    )
    aliases = parse_assembly_report_aliases(report)
    assert aliases["MT"] == ["NC_MT"]
    assert aliases["23"] == ["NC_23"]


def test_build_reference_registry_runs_from_verified_cached_package(tmp_path):
    project = tmp_path / "project"
    (project / ".git").mkdir(parents=True)
    config_dir = project / "configs"
    config_dir.mkdir()
    output = project / "results" / "codex_reanalysis_v1" / "01_reference"
    downloads = output / "downloads"
    downloads.mkdir(parents=True)

    fasta_name = "GCF_test_genomic.fna.gz"
    gff_name = "GCF_test_genomic.gff.gz"
    assembly_name = "GCF_test_assembly_report.txt"
    annotation_name = "GCF_test_annotation_report.xml"
    fasta_bytes = gzip.compress(b">NC_037350.1 chromosome 23\n" + b"A" * 9 + b"ACGTACGTACG\n")
    gff_bytes = gzip.compress(
        b"##gff-version 3\n"
        b"NC_037350.1\tRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene-DRB3;gene=DRB3\n"
    )
    assembly_bytes = (
        "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
        "23\tassembled-molecule\t23\tChromosome\tCM23\t=\tNC_037350.1\tPrimary Assembly\t20\tchr23\n"
    ).encode()
    resources = {
        fasta_name: fasta_bytes,
        gff_name: gff_bytes,
        assembly_name: assembly_bytes,
        annotation_name: b"<annotation />\n",
    }
    for name, content in resources.items():
        (downloads / name).write_bytes(content)
    (downloads / "md5checksums.txt").write_text(
        "".join(f"{md5(content).hexdigest()}  ./{name}\n" for name, content in resources.items()),
        encoding="ascii",
    )
    target_spec = config_dir / "targets.tsv"
    target_spec.write_text(
        "name\tgene_symbols\tcategory\tabsent_policy\tprobe_id\n"
        "DRB3\tDRB3\tclass_II\terror\t\n"
        "DQB2\tDQB2\tclass_II\tprobe_only\tIPD-DQB2\n",
        encoding="utf-8",
    )
    config = config_dir / "example.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "tools": {"samtools": "samtools", "minimap2": "minimap2"},
                "reference": {
                    "assembly": "ARS-UCD2.0",
                    "accession": "GCF_002263795.3",
                    "annotation_release": "GCF_002263795.3-RS_2024_12",
                    "base_url": "https://ftp.ncbi.nlm.nih.gov/test-release",
                    "checksum_file": "md5checksums.txt",
                    "fasta_file": fasta_name,
                    "gff_file": gff_name,
                    "assembly_report_file": assembly_name,
                    "annotation_report_file": annotation_name,
                    "target_spec_relative_path": "configs/targets.tsv",
                    "comparison_contigs": ["23"],
                },
                "outputs": {"root": "results/codex_reanalysis_v1"},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    commands = []

    def fake_runner(command):
        commands.append(command)
        if command[1] == "faidx":
            Path(command[command.index("-o") + 1]).write_text("index\n")
        else:
            Path(command[command.index("-d") + 1]).write_text("index\n")

    result = build_reference_registry(
        config,
        output,
        command_runner=fake_runner,
        retrieval_date="2026-07-31",
    )

    assert result["status"] == "incomplete"
    assert len(commands) == 2
    registry = (output / "coordinate_registry.tsv").read_text(encoding="utf-8")
    assert "DRB3" in registry
    assert "haplotype_probe" in registry
    manifest_text = (output / "reference_manifest.tsv").read_text(encoding="utf-8")
    assert "\ttargets.bed\tgenerated-coordinate-registry\t" in manifest_text
    assert (output / "provenance.json").is_file()


def test_reference_registry_script_supports_the_planned_direct_cli():
    result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "scripts" / "codex_reanalysis" / "reference_registry.py"),
            "--help",
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def test_fasta_contig_comparison_writes_alignment_and_decision_table(tmp_path):
    old_fasta = tmp_path / "old.fa"
    new_fasta = tmp_path / "new.fa"
    old_fasta.write_text(">chr23 old\nACGT\n", encoding="ascii")
    new_fasta.write_text(">NC_23 new\nACGT\n", encoding="ascii")
    commands = []

    def fake_alignment(command, output_path):
        commands.append(command)
        output_path.write_text("old\t4\t0\t4\t+\tnew\t4\t0\t4\t4\t4\t60\n", encoding="ascii")

    comparisons = compare_fasta_contigs(
        old_fasta=old_fasta,
        new_fasta=new_fasta,
        old_aliases={"23": ["NC_23"], "chr23": ["NC_23"]},
        aliases={"23": ["NC_23"], "chr23": ["NC_23"]},
        contigs=["23"],
        output_dir=tmp_path / "comparison",
        minimap2="minimap2",
        alignment_runner=fake_alignment,
    )

    assert comparisons[0]["status"] == "identical"
    assert comparisons[0]["old_accession"] == "NC_23"
    assert comparisons[0]["old_fasta_name"] == "chr23"
    assert comparisons[0]["alignment_status"] == "completed"
    assert commands[0][:4] == ["minimap2", "-x", "asm5", "-c"]
    table = (tmp_path / "comparison" / "sequence_comparison.tsv").read_text(encoding="utf-8")
    assert "old_sha256\tnew_sha256\tstatus\talignment_status" in table


def test_ipd_alignment_summary_keeps_loci_with_no_alignment(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "ipd_id\tallele_name\tlocus\tseq_length\tsource\n"
        "ID1\tBoLA-DQB*001\tBoLA-DQB\t100\tpanel.fa\n"
        "ID2\tBoLA-NC10*001\tBoLA-NC10\t100\tpanel.fa\n",
        encoding="utf-8",
    )
    paf = tmp_path / "panel.paf"
    paf.write_text(
        "BoLA-DQB*001\t100\t0\t100\t+\tNC_23\t1000\t10\t110\t98\t100\t60\n",
        encoding="ascii",
    )
    summary = summarize_ipd_alignments(manifest, paf, tmp_path / "summary.tsv")
    assert summary["BoLA-DQB"]["status"] == "aligned"
    assert summary["BoLA-DQB"]["best_identity"] == pytest.approx(0.98)
    assert summary["BoLA-NC10"]["status"] == "no_alignment"


def test_ipd_alignment_summary_rejects_partial_hits_as_target_evidence(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "allele_name\tlocus\nID1\tBoLA-NC1\n",
        encoding="utf-8",
    )
    paf = tmp_path / "panel.paf"
    paf.write_text(
        "ID1\t100\t0\t28\t+\tNC_23\t1000\t10\t38\t28\t28\t60\n",
        encoding="ascii",
    )

    summary = summarize_ipd_alignments(manifest, paf, tmp_path / "summary.tsv")

    assert summary["BoLA-NC1"]["status"] == "partial_alignment"
    assert summary["BoLA-NC1"]["best_identity"] == pytest.approx(1.0)
    assert summary["BoLA-NC1"]["best_coverage"] == pytest.approx(0.28)


def test_tool_version_validation_rejects_a_non_pinned_binary():
    assert validate_tool_version("samtools 1.23.1\n", "1.23.1", "samtools") == "samtools 1.23.1"
    with pytest.raises(ReferenceRegistryError, match="version mismatch"):
        validate_tool_version("samtools 1.19.2\n", "1.23.1", "samtools")
    with pytest.raises(ReferenceRegistryError, match="version mismatch"):
        validate_tool_version("2.310-r9999\n", "2.31-r1302", "minimap2")


def test_final_passed_status_is_not_published_when_provenance_fails(tmp_path, monkeypatch):
    project = tmp_path / "project"
    (project / ".git").mkdir(parents=True)
    output = project / "results" / "codex_reanalysis_v1" / "01_reference"
    output.mkdir(parents=True)
    status_path = output / "reference_status.json"
    status_path.write_text('{"status": "running"}\n', encoding="utf-8")

    def fail_provenance(*args, **kwargs):
        raise RuntimeError("provenance failed")

    monkeypatch.setattr(reference_registry_module, "write_provenance", fail_provenance)
    with pytest.raises(RuntimeError, match="provenance failed"):
        reference_registry_module._publish_reference_completion(
            output,
            {"status": "passed"},
            {"config_hash": "test"},
        )

    assert '"running"' in status_path.read_text(encoding="utf-8")
