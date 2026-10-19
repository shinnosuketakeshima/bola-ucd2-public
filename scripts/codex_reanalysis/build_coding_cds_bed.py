"""Build the combined CDS restriction BED for `call_coding_variants.py`'s call mode.

`call_coding_variants.py --mode call` restricts `bcftools mpileup` to one `--cds-bed`
per sample, but the BTN + CNTN3 + ABHD16A/PRR3 exon spans live in three separate legacy
files loaded by `--mode merge` (`btn_cds_models.json`, `cntn3_cds.bed`,
`bta23_gwas_cds.bed`). This script reuses those same loaders and the same
`_GENES_IN_SCOPE` gene set so the region a sample is called against and the region its
calls are later scored against can never drift apart -- there is exactly one place gene
scope is defined.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.call_coding_variants import (
    _BTN_STRAND,
    _GENES_IN_SCOPE,
    load_cds_models_from_bed,
    load_cds_models_from_json,
    load_cntn3_cds_models,
)


def build_models(btn_cds_json: Path, cntn3_cds_bed: Path, bta23_gwas_cds_bed: Path) -> dict[str, dict]:
    """Load and scope-filter exon models exactly as call_coding_variants.py's merge mode does."""
    models: dict[str, dict] = {}
    models.update(load_cds_models_from_json(btn_cds_json, _BTN_STRAND))
    models.update(load_cntn3_cds_models(cntn3_cds_bed))
    models.update(load_cds_models_from_bed(bta23_gwas_cds_bed))
    return {gene: model for gene, model in models.items() if gene in _GENES_IN_SCOPE}


def render_bed(models: dict[str, dict]) -> str:
    """Emit one 0-based half-open BED row per exon, sorted by (contig, start, end)."""
    rows = [
        (model["contig"], start - 1, end, gene)
        for gene, model in models.items()
        for start, end in model["exons"]
    ]
    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    return "".join(f"{contig}\t{start}\t{end}\t{gene}\n" for contig, start, end, gene in rows)


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--btn-cds-json", type=Path, required=True)
    parser.add_argument("--cntn3-cds-bed", type=Path, required=True)
    parser.add_argument("--bta23-gwas-cds-bed", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    models = build_models(args.btn_cds_json, args.cntn3_cds_bed, args.bta23_gwas_cds_bed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_bed(models), encoding="utf-8")
    exon_count = sum(len(model["exons"]) for model in models.values())
    print(f"Coding CDS restriction BED written: {args.output} ({len(models)} genes, {exon_count} exons)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
