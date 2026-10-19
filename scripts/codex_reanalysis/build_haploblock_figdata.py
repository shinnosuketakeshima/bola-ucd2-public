"""Build haploblock figure data from UCD2.0 MHC phase-block outputs."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

# Corrected BTA23 landmarks (Mb) for the haploblock map overlay.
LOCI = [
    {"name": "DQ (class II)", "s": 25.58303, "e": 25.843091},
    {"name": "DRB3", "s": 25.723716, "e": 25.734818},
    {"name": "BTN", "s": 25.86, "e": 27.15},
    {"name": "rs29026690 (class III)", "s": 27.603651, "e": 27.603651},
    {"name": "class I", "s": 28.493421, "e": 28.724293},
    {"name": "rs17872126 (class I prox)", "s": 28.419662, "e": 28.419662},
]

DISPLAY_RANGE = [25.4, 28.9]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(lines, delimiter="\t"))


def is_016_carrier(hap1: str, hap2: str) -> int:
    return int("*016:01" in {hap1.strip(), hap2.strip()})


def load_blocks(phase_blocks_path: Path) -> tuple[list[str], dict[str, list[dict[str, float]]]]:
    rows = read_tsv(phase_blocks_path)
    order: list[str] = []
    blocks: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        sample = row["sample_id"].strip()
        if sample not in order:
            order.append(sample)
        blocks[sample].append(
            {
                "s": round(int(row["start"]) / 1_000_000, 6),
                "e": round(int(row["end"]) / 1_000_000, 6),
            }
        )
    return order, dict(blocks)


def load_boundary_hist(hotspots_path: Path) -> list[dict[str, float | int]]:
    rows = read_tsv(hotspots_path)
    hist = []
    for row in rows:
        start = int(row["window_start"])
        end = int(row["window_end"])
        hist.append(
            {
                "mb": round((start + end) / 2 / 1_000_000, 2),
                "n": int(row["n_breaks"]),
            }
        )
    hist.sort(key=lambda item: item["mb"])
    return hist


def load_drb3_from_alleles(drb3_path: Path) -> dict[str, str]:
    rows = read_tsv(drb3_path)
    return {
        row["sample_id"].strip(): f"{row['drb3_hap1'].strip()}/{row['drb3_hap2'].strip()}"
        for row in rows
    }


def load_gwas(marker_path: Path, drb3_path: Path) -> dict[str, dict[str, object]]:
    markers = read_tsv(marker_path)
    drb3_rows = {row["sample_id"].strip(): row for row in read_tsv(drb3_path)}
    by_sample: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in markers:
        sample = row["sample_id"].strip()
        marker = row["marker_id"].strip()
        by_sample[sample][marker] = row["primary_genotype"].strip()

    gwas: dict[str, dict[str, object]] = {}
    for sample, genotypes in by_sample.items():
        drb3 = drb3_rows[sample]
        gwas[sample] = {
            "rs29026690": genotypes.get("rs29026690", ""),
            "rs17872126": genotypes.get("rs17872126", ""),
            "c016": is_016_carrier(drb3["drb3_hap1"], drb3["drb3_hap2"]),
        }
    return gwas


def build_payload(repo_root: Path) -> dict[str, object]:
    phase_root = repo_root / "results/codex_reanalysis_v1/04b_mhc_phase_blocks"
    samples, blocks = load_blocks(phase_root / "phase_blocks.tsv")
    return {
        "loci": LOCI,
        "samples": samples,
        "blocks": blocks,
        "boundary_hist": load_boundary_hist(phase_root / "all_boundary_hotspots_100kb.tsv"),
        "drb3": load_drb3_from_alleles(repo_root / "data/meta/drb3_alleles.tsv"),
        "gwas": load_gwas(
            repo_root / "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv",
            repo_root / "data/meta/drb3_alleles.tsv",
        ),
        "range": DISPLAY_RANGE,
        "source": {
            "phase_blocks": str(phase_root / "phase_blocks.tsv"),
            "boundary_hotspots": str(phase_root / "all_boundary_hotspots_100kb.tsv"),
            "markers": str(repo_root / "results/codex_reanalysis_v1/04_small_variants/marker_genotypes.tsv"),
        },
    }


def write_new_json(path: Path, payload: dict[str, object]) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ": ")) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/figures/haploblock_figdata.json"),
    )
    parser.add_argument("--force", action="store_true", help="overwrite an existing output file")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    output = args.output if args.output.is_absolute() else repo_root / args.output
    if output.exists() and not args.force:
        raise FileExistsError(f"refusing to overwrite {output}; pass --force to replace")

    payload = build_payload(repo_root)
    if args.force and output.exists():
        output.unlink()
    write_new_json(output, payload)
    peak = max(payload["boundary_hist"], key=lambda item: item["n"])
    print(
        f"Wrote {output} ({len(payload['samples'])} samples; "
        f"peak boundary {peak['mb']} Mb, n={peak['n']})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
