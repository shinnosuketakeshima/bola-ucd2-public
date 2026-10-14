"""Regenerate docs/tables/table2_dq_copynumber.tsv and table3_btn_gd_classI.tsv.

Both tables previously carried ARS-UCD1.2-generation values (three-fixed-control
depth normalization) even after the manuscript prose and evidence matrix moved to
the ARS-UCD2.0 codex_reanalysis pipeline: table3's BTNL2 cohort median (1.361) is
the legacy value named in docs/codex_reanalysis/LEGACY_COMPARISON.md, not the
UCD2.0 value (1.3239) already quoted in docs/RESULTS_EN.md's own prose. This
script rebuilds both tables directly from the pipeline's own status-gated outputs
(03_depth_cnv/normalized_depth_matrix.tsv at MAPQ>=20, 06_local_assembly/dq_dr's
recovery TSV) so table and text cite the same numbers again. See
docs/codex_reanalysis/UCD2_MIGRATION_STATUS.md for the fuller narrative.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

SAMPLE_ORDER = [
    "ID01",
    "ID02",
    "ID03",
    "ID04",
    "ID05",
    "ID06",
    "ID07",
    "ID08",
    "ID09",
    "ID10",
    "ID11",
    "ID12",
    "ID13",
]

# table3 column name -> normalized_depth_matrix.tsv locus name
BTN_LOCI = {
    "BTNL2": "BTNL2",
    "BTN2A1L_LOC104969893": "BTN_LOC104969893",
    "BTN3A2": "BTN_LOC615223",
    "BTN2A2": "BTN_LOC107131749",
    "BTN3A3": "BTN_LOC786760",
}
CLASS_I_LOCI = {
    "classI_NC1": "NC1",
    "classI_JSP.1": "JSP1",
    # Target "BOLA" (not "classI_BoLA_A") is classical class I BoLA at
    # 28,720,499-28,724,294 (coordinate_registry.tsv, Table S1) with real IPD
    # BoLA-1..6 probes; confirmed against RESULTS_EN.md's own prose ("classical
    # BoLA ... roughly 2.3-fold elevation in ID08 ... ratio 2.30"), which only
    # "BOLA" reproduces (DEPTH_CNV_REPORT.md max 2.3). "classI_BoLA_A" is a
    # separate, unprobed alias target and is not classical BoLA.
    "classI_BOLA": "BOLA",
}
DQ_LOCI = {
    "depth_DQA2": "DQA2",
    "depth_DQB1": "DQB1",
    "depth_DRB3": "DRB3",
}

CALL_MAP = {
    "elevated_relative_to_cohort": "dup",
    "reduced_relative_to_cohort": "del",
    "consistent_with_cohort": "normal",
    "present_at_otherwise_unrepresented_locus": "present",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def index_depth(rows: list[dict[str, str]], mapq: str) -> dict[tuple[str, str], dict[str, str]]:
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if row["mapq"] != mapq:
            continue
        index[(row["locus"], row["sample_id"])] = row
    return index


def dq_interpretation(record: dict[str, str]) -> str:
    raw = record["interpretation"]
    dqa = int(record["n_dqa_clusters"])
    dqb = int(record["n_dqb_clusters"])
    if raw == "under_assembled":
        # Exact string "under-assembly" is matched by scripts/chatgpt_figures/lib/data.mjs's
        # underAssembly flag (dqRow.interpretation === 'under-assembly'); keep it unsuffixed
        # here and add "(excluded)" only in the human-readable manuscript table.
        return "under-assembly"
    if raw == "diploid_baseline":
        return "diploid-baseline"
    if raw == "expanded_dq":
        if dqa > 2 and dqb > 2:
            return "DQ expanded (DQA+DQB)"
        if dqa > 2:
            return "DQ expanded (DQA)"
        if dqb > 2:
            return "DQ expanded (DQB)"
        return "DQ expanded"
    raise ValueError(f"unrecognized dq_dr_recovery interpretation: {raw!r}")


def build_table2(depth_mapq20: dict[tuple[str, str], dict[str, str]], recovery: dict[str, dict[str, str]]) -> tuple[list[str], list[list[str]]]:
    header = [
        "sample_id",
        "DRB3",
        "depth_DQA2",
        "depth_DQB1",
        "depth_DRB3",
        "asm_DQA_clusters",
        "asm_DQB_clusters",
        "interpretation",
    ]
    rows = []
    for sample_id in SAMPLE_ORDER:
        recovery_row = recovery[sample_id]
        row = [sample_id, recovery_row["pcr_sbt"]]
        for column, locus in DQ_LOCI.items():
            row.append(depth_mapq20[(locus, sample_id)]["normalized_ratio"])
        row.append(recovery_row["n_dqa_clusters"])
        row.append(recovery_row["n_dqb_clusters"])
        row.append(dq_interpretation(recovery_row))
        rows.append(row)
    return header, rows


def build_table3(depth_mapq20: dict[tuple[str, str], dict[str, str]]) -> tuple[list[str], list[list[str]]]:
    btn_columns = list(BTN_LOCI.items())
    header = (
        ["sample_id"]
        + [column for column, _ in btn_columns]
        + [f"{column}_call" for column, _ in btn_columns]
        + list(CLASS_I_LOCI.keys())
    )
    rows = []
    for sample_id in SAMPLE_ORDER:
        ratios = [depth_mapq20[(locus, sample_id)]["normalized_ratio"] for _, locus in btn_columns]
        calls = [
            CALL_MAP[depth_mapq20[(locus, sample_id)]["interpretation"]]
            for _, locus in btn_columns
        ]
        class_i = [depth_mapq20[(locus, sample_id)]["normalized_ratio"] for locus in CLASS_I_LOCI.values()]
        rows.append([sample_id, *ratios, *calls, *class_i])
    return header, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.project_root.resolve()

    depth_rows = read_tsv(root / "results/codex_reanalysis_v1/03_depth_cnv/normalized_depth_matrix.tsv")
    depth_mapq20 = index_depth(depth_rows, "20")

    recovery_rows = read_tsv(root / "results/codex_reanalysis_v1/06_local_assembly/dq_dr/dq_dr_recovery.tsv")
    recovery = {row["sample_id"]: row for row in recovery_rows}
    missing = [sample for sample in SAMPLE_ORDER if sample not in recovery]
    if missing:
        raise SystemExit(f"dq_dr_recovery.tsv missing samples: {missing}")

    header2, rows2 = build_table2(depth_mapq20, recovery)
    write_tsv(root / "docs/tables/table2_dq_copynumber.tsv", header2, rows2)

    header3, rows3 = build_table3(depth_mapq20)
    write_tsv(root / "docs/tables/table3_btn_gd_classI.tsv", header3, rows3)

    print("Rebuilt docs/tables/table2_dq_copynumber.tsv and table3_btn_gd_classI.tsv from ARS-UCD2.0 outputs (MAPQ>=20).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
