from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from scripts.codex_reanalysis.package_manuscript import (
    ManuscriptPackageError,
    package_manuscript,
)


def test_package_manuscript_builds_reviewed_bundle(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()

    evidence_dir = repo / "results/codex_reanalysis_v1/07_evidence_audit"
    evidence_dir.mkdir(parents=True)
    (evidence_dir / "evidence_matrix_status.json").write_text(
        json.dumps({"status": "passed", "claim_count": 1, "accepted_count": 1, "limited_count": 0, "rejected_count": 0})
        + "\n",
        encoding="utf-8",
    )

    figure_dir = repo / "docs/chatgpt_figures_submission"
    figure_dir.mkdir(parents=True)
    for name in (
        "Figure1_BoLA_framework",
        "Figure2_DQ_reference_and_assembly",
        "Figure3_BTN_depth_heterogeneity",
        "Figure4_classI_depth_variation",
        "Figure5_DQ_by_DRB3_016",
        "Figure6_GWAS_marker_context",
        "Figure7_correlation_phase_unresolved",
        "Figure8_WhatsHap_phase_blocks",
        "FigureS1_classIII_exploratory",
    ):
        (figure_dir / f"{name}.svg").write_text("<svg viewBox='0 0 10 10'></svg>\n", encoding="utf-8")

    manifest_dir = repo / "results/chatgpt_figure_reanalysis"
    manifest_dir.mkdir(parents=True)
    (manifest_dir / "export_manifest.json").write_text("{}\n", encoding="utf-8")
    (manifest_dir / "figure_qa.json").write_text('{"status":"pass","checks":[]}\n', encoding="utf-8")
    (manifest_dir / "plot_data.json").write_text("{}\n", encoding="utf-8")

    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True)

    output = repo / "results/codex_reanalysis_v1/08_manuscript"
    payload = package_manuscript(repo, output)
    assert payload["status"] == "passed"
    assert (output / "figures/Figure1_BoLA_framework.svg").is_file()
    assert (output / "manuscript_status.json").is_file()
    assert (output / "provenance.json").is_file()


def test_package_manuscript_refuses_failed_figure_qa(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    evidence_dir = repo / "results/codex_reanalysis_v1/07_evidence_audit"
    evidence_dir.mkdir(parents=True)
    (evidence_dir / "evidence_matrix_status.json").write_text('{"status":"passed"}\n', encoding="utf-8")
    manifest_dir = repo / "results/chatgpt_figure_reanalysis"
    manifest_dir.mkdir(parents=True)
    (manifest_dir / "figure_qa.json").write_text('{"status":"fail"}\n', encoding="utf-8")
    with pytest.raises(ManuscriptPackageError, match="figure QA"):
        package_manuscript(repo, repo / "results/codex_reanalysis_v1/08_manuscript")
