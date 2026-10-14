from hashlib import sha256
from pathlib import Path

import pytest

from scripts.codex_reanalysis.coordinate_registry import (
    CoordinateError,
    TargetSpec,
    apply_probe_alignment_statuses,
    generate_target_records,
    gff_to_bed,
    interval_length,
    load_fasta_sequences,
    load_target_specs,
    parse_gene_features,
    resolve_contig,
    write_coordinate_outputs,
)


def test_gff_one_based_closed_to_bed_zero_based_half_open():
    assert gff_to_bed(start=10, end=20) == (9, 20)


def test_interval_end_is_exclusive_in_bed():
    assert interval_length(9, 20) == 11


@pytest.mark.parametrize("start,end", [(0, 1), (2, 1), (-1, 5)])
def test_invalid_gff_intervals_are_rejected(start, end):
    with pytest.raises(CoordinateError, match="GFF interval"):
        gff_to_bed(start=start, end=end)


def test_alias_resolution_rejects_ambiguous_aliases():
    with pytest.raises(CoordinateError, match="ambiguous"):
        resolve_contig("23", {"23": ["NC_037350.1", "chr23"]})


def test_alias_resolution_requires_a_known_alias():
    with pytest.raises(CoordinateError, match="unknown"):
        resolve_contig("BTA99", {"23": ["NC_037350.1"]})


def test_gff_features_generate_provenanced_bed_coordinates(tmp_path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "NC_037350.1\tRefSeq\tgene\t10\t20\t.\t-\t.\tID=gene-DRB3;gene=DRB3\n",
        encoding="utf-8",
    )
    features = parse_gene_features(gff)
    records = generate_target_records(
        features,
        [TargetSpec(name="DRB3", gene_symbols=("DRB3",), category="class_II")],
        assembly_accession="GCF_002263795.3",
        annotation_release="GCF_002263795.3-RS_2024_12",
        aliases={"23": ["NC_037350.1"]},
        sequences={"NC_037350.1": "A" * 9 + "ACGTACGTACG" + "T" * 10},
    )

    assert len(records) == 1
    record = records[0]
    assert (record.gff_start, record.gff_end) == (10, 20)
    assert (record.bed_start, record.bed_end) == (9, 20)
    assert record.strand == "-"
    assert record.feature_id == "gene-DRB3"
    assert record.sequence_sha256 == sha256(b"ACGTACGTACG").hexdigest()


def test_duplicate_gene_features_are_not_chosen_silently(tmp_path):
    gff = tmp_path / "duplicate.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "NC_037350.1\tRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene-a;gene=DQA1\n"
        "NC_037350.1\tRefSeq\tgene\t30\t40\t.\t+\t.\tID=gene-b;gene=DQA1\n",
        encoding="utf-8",
    )

    with pytest.raises(CoordinateError, match="non-unique"):
        generate_target_records(
            parse_gene_features(gff),
            [TargetSpec(name="DQA1", gene_symbols=("DQA1",), category="class_II")],
            assembly_accession="GCF_002263795.3",
            annotation_release="release",
            aliases={"23": ["NC_037350.1"]},
            sequences={"NC_037350.1": "A" * 100},
        )


def test_absent_dqb2_is_probe_only_not_an_invented_reference_gene():
    records = generate_target_records(
        [],
        [
            TargetSpec(
                name="DQB2",
                gene_symbols=("DQB2",),
                category="class_II",
                absent_policy="probe_only",
                probe_id="IPD-MHC-DQB2",
            )
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_037350.1"]},
        sequences={"NC_037350.1": "A" * 100},
    )

    assert len(records) == 1
    assert records[0].target_type == "haplotype_probe"
    assert records[0].feature_id == ""
    assert records[0].probe_alignment_status == "pending"


def test_target_specification_rejects_legacy_coordinate_columns(tmp_path):
    spec = tmp_path / "targets.tsv"
    spec.write_text(
        "name\tgene_symbols\tcategory\tabsent_policy\tprobe_id\tstart\n"
        "DRB3\tDRB3\tclass_II\terror\t\t25723716\n",
        encoding="utf-8",
    )
    with pytest.raises(CoordinateError, match="coordinate columns"):
        load_target_specs(spec)


def test_fasta_loader_selects_named_records_without_loading_unrequested_contigs(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(
        ">NC_1 description\nAAAA\n>NC_2\nCC\nGG\n>NC_3\nTTTT\n",
        encoding="ascii",
    )
    assert load_fasta_sequences(fasta, {"NC_2"}) == {"NC_2": "CCGG"}


def test_coordinate_outputs_exclude_probe_only_targets_from_bed_and_fasta(tmp_path):
    reference = generate_target_records(
        [
            parse_gene_features(_write_one_gene_gff(tmp_path))[0],
        ],
        [TargetSpec(name="DRB3", gene_symbols=("DRB3",), category="class_II")],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_037350.1"]},
        sequences={"NC_037350.1": "A" * 9 + "ACGTACGTACG"},
    )[0]
    probe = generate_target_records(
        [],
        [
            TargetSpec(
                name="DQB2",
                gene_symbols=("DQB2",),
                category="class_II",
                absent_policy="probe_only",
                probe_id="IPD-DQB2",
            )
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_037350.1"]},
        sequences={"NC_037350.1": "A" * 20},
    )[0]

    paths = write_coordinate_outputs(
        [reference, probe],
        sequences={"NC_037350.1": "A" * 9 + "ACGTACGTACG"},
        output_dir=tmp_path / "out",
    )

    registry_rows = paths["registry"].read_text(encoding="utf-8").splitlines()
    assert len(registry_rows) == 3
    assert paths["bed"].read_text(encoding="utf-8") == "NC_037350.1\t9\t20\tDRB3\n"
    assert paths["fasta"].read_text(encoding="ascii") == ">DRB3|NC_037350.1:10-20\nACGTACGTACG\n"


def test_probe_statuses_use_available_ipd_loci_without_inventing_dqb2():
    records = generate_target_records(
        [],
        [
            TargetSpec(
                name="DQB2",
                gene_symbols=("DQB2",),
                category="class_II",
                absent_policy="probe_only",
                probe_id="IPD-MHC-DQB2",
                probe_loci=(),
            )
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_23"]},
        sequences={"NC_23": "A"},
    )
    updated = apply_probe_alignment_statuses(
        records,
        {"BoLA-DQB": "aligned", "BoLA-NC10": "no_alignment"},
    )
    assert updated[0].target_type == "haplotype_probe"
    assert updated[0].feature_id == ""
    assert updated[0].probe_alignment_status == "probe_unavailable"


def test_probe_statuses_require_an_explicit_target_to_panel_locus_mapping():
    records = generate_target_records(
        [],
        [
            TargetSpec(
                name="DQB_ref",
                gene_symbols=("DQB_ref",),
                category="class_II",
                absent_policy="probe_only",
                probe_id="IPD-MHC-DQB",
                probe_loci=("BoLA-DQB",),
            )
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_23"]},
        sequences={"NC_23": "A"},
    )
    updated = apply_probe_alignment_statuses(records, {"BoLA-DQB": "aligned"})
    assert updated[0].probe_alignment_status == "panel_aligned"


def test_shared_panel_locus_is_labeled_as_family_not_target_specific_evidence():
    records = generate_target_records(
        [],
        [
            TargetSpec(
                name=name,
                gene_symbols=(name,),
                category="class_II",
                absent_policy="probe_only",
                probe_id=f"IPD-{name}",
                probe_loci=("BoLA-DQA",),
            )
            for name in ("DQA1", "DQA2")
        ],
        assembly_accession="GCF_002263795.3",
        annotation_release="release",
        aliases={"23": ["NC_23"]},
        sequences={"NC_23": "A"},
    )
    updated = apply_probe_alignment_statuses(records, {"BoLA-DQA": "aligned"})
    assert {row.probe_alignment_status for row in updated} == {"panel_family_aligned"}


def test_target_specification_rejects_duplicate_target_names(tmp_path):
    spec = tmp_path / "targets.tsv"
    spec.write_text(
        "name\tgene_symbols\tcategory\tabsent_policy\tprobe_id\tprobe_loci\n"
        "DRB3\tDRB3\tclass_II\terror\t\t\n"
        "DRB3\tDRB3\tclass_II\terror\t\t\n",
        encoding="utf-8",
    )
    with pytest.raises(CoordinateError, match="duplicate target"):
        load_target_specs(spec)


def test_tracked_variant_targets_are_coordinate_free_and_have_separate_host_genes():
    specs = load_target_specs(
        Path(__file__).resolve().parents[2] / "configs" / "codex_reanalysis.targets.tsv"
    )
    by_name = {spec.name: spec for spec in specs}
    for rs_id, host_name in (
        ("rs110616206", "host_CNTN3"),
        ("rs29026690", "host_ABHD16A"),
        ("rs17872126", "host_PRR3"),
    ):
        assert by_name[rs_id].absent_policy == "probe_only"
        assert by_name[rs_id].gene_symbols == (rs_id,)
        assert by_name[host_name].absent_policy == "error"


def _write_one_gene_gff(tmp_path):
    gff = tmp_path / "one.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "NC_037350.1\tRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene-DRB3;gene=DRB3\n",
        encoding="utf-8",
    )
    return gff
