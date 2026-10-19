from __future__ import annotations

from pathlib import Path

import pytest

from scripts.codex_reanalysis.build_haploblock_figdata import build_payload


REPO_ROOT = Path(__file__).resolve().parents[2]
PHASE_BLOCKS = REPO_ROOT / "results/codex_reanalysis_v1/04b_mhc_phase_blocks/phase_blocks.tsv"


@pytest.mark.skipif(not PHASE_BLOCKS.exists(), reason="needs the untracked 04b phase-block stage output")
def test_build_payload_matches_ucd2_phase_block_peak() -> None:
    repo_root = REPO_ROOT
    payload = build_payload(repo_root)
    assert len(payload["samples"]) == 13
    peak = max(payload["boundary_hist"], key=lambda item: item["n"])
    assert peak["mb"] == 26.25
    assert peak["n"] == 21
    assert payload["blocks"]["ID01"][0]["s"] < payload["blocks"]["ID01"][0]["e"]
