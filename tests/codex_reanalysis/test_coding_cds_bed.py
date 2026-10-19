from __future__ import annotations

import json
from pathlib import Path

from scripts.codex_reanalysis.build_coding_cds_bed import build_models, render_bed


def _write_fixtures(tmp_path: Path) -> tuple[Path, Path, Path]:
    btn_json = tmp_path / "btn_cds_models.json"
    btn_json.write_text(
        json.dumps({"BTNL2": {"chrom": "NC_037350.1", "cds": [[100, 150], [200, 250]]}}),
        encoding="utf-8",
    )
    cntn3_bed = tmp_path / "cntn3_cds.bed"
    cntn3_bed.write_text("NC_037349.1\t99\t120\tCNTN3_CDS\n", encoding="utf-8")
    gwas_bed = tmp_path / "bta23_gwas_cds.bed"
    gwas_bed.write_text(
        "NC_037350.1\t9999\t10050\tABHD16A\t+\n"
        "NC_037350.1\t20000\t20060\tPRR3\t-\n"
        "NC_037350.1\t30000\t30060\tBAG6\t+\n",  # out of scope, must be dropped
        encoding="utf-8",
    )
    return btn_json, cntn3_bed, gwas_bed


def test_build_models_renames_cntn3_and_drops_out_of_scope_genes(tmp_path):
    btn_json, cntn3_bed, gwas_bed = _write_fixtures(tmp_path)
    models = build_models(btn_json, cntn3_bed, gwas_bed)
    assert set(models) == {"BTNL2", "CNTN3", "ABHD16A", "PRR3"}
    assert "BAG6" not in models


def test_render_bed_is_zero_based_half_open_and_sorted(tmp_path):
    btn_json, cntn3_bed, gwas_bed = _write_fixtures(tmp_path)
    models = build_models(btn_json, cntn3_bed, gwas_bed)
    bed_text = render_bed(models)
    lines = bed_text.splitlines()
    # CNTN3 is on NC_037349.1, everything else on NC_037350.1; contig-sorted first.
    assert lines[0] == "NC_037349.1\t99\t120\tCNTN3"
    assert "NC_037350.1\t99\t150\tBTNL2" in lines
    assert "NC_037350.1\t199\t250\tBTNL2" in lines
    assert "NC_037350.1\t9999\t10050\tABHD16A" in lines
    assert "NC_037350.1\t20000\t20060\tPRR3" in lines
    assert len(lines) == 5
