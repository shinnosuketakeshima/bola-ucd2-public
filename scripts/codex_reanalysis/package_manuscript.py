"""Assemble the stageable manuscript artifact bundle under 08_manuscript/.

Consumes the evidence-audit gate, regenerated figure SVGs, haploblock Fig. 7
data, and the manuscript number audit. Writes a machine-readable inventory plus
status/provenance for the codex_reanalysis pipeline's final stage.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file

FIGURE_NAMES = (
    "Figure1_BoLA_framework",
    "Figure2_DQ_reference_and_assembly",
    "Figure3_BTN_depth_heterogeneity",
    "Figure4_classI_depth_variation",
    "Figure5_DQ_by_DRB3_016",
    "Figure6_GWAS_marker_context",
    "Figure7_correlation_phase_unresolved",
    "Figure8_WhatsHap_phase_blocks",
    "FigureS1_classIII_exploratory",
)


class ManuscriptPackageError(ValueError):
    """Raised when a prerequisite artifact is missing or fails QA."""


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise ManuscriptPackageError(f"missing required {label}: {path}")
    return path


def _validate_evidence_gate(repo_root: Path) -> dict:
    status_path = repo_root / "results/codex_reanalysis_v1/07_evidence_audit/evidence_matrix_status.json"
    payload = _read_json(_require_file(status_path, "evidence audit status"))
    if payload.get("status") not in {"passed", "passed_with_scope_limitation"}:
        raise ManuscriptPackageError(
            f"evidence audit status is not passed: {payload.get('status')!r}"
        )
    return payload


def _validate_figure_qa(repo_root: Path) -> dict:
    qa_path = repo_root / "results/chatgpt_figure_reanalysis/figure_qa.json"
    payload = _read_json(_require_file(qa_path, "figure QA report"))
    if payload.get("status") != "pass":
        raise ManuscriptPackageError(f"figure QA status is not pass: {payload.get('status')!r}")
    return payload


def _copy_file(source: Path, destination: Path, force: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not force:
            raise FileExistsError(f"refusing to overwrite {destination}")
        destination.unlink()
    shutil.copy2(source, destination)


def _sanitize_export_manifest(source_path: Path) -> str:
    payload = _read_json(source_path)
    figures = []
    for entry in payload.get("figures", []):
        name = Path(entry["path"]).name
        figures.append(
            {
                "format": entry["format"],
                "path": f"results/codex_reanalysis_v1/08_manuscript/figures/{name}",
                "bytes": entry["bytes"],
                "sha256": entry["sha256"],
            }
        )
    return json.dumps({"figures": figures}, indent=2, ensure_ascii=False) + "\n"


def _write_new_text(path: Path, text: str, *, force: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if not force:
            raise FileExistsError(f"refusing to overwrite {path}")
        path.unlink()
    path.write_text(text, encoding="utf-8")


def build_inventory(
    *,
    evidence_status: dict,
    figure_qa: dict,
    copied_files: Sequence[Path],
) -> str:
    lines = [
        "# Manuscript artifact package (ARS-UCD2.0 reanalysis)",
        "",
        "Reviewed, anonymized figure exports and QA manifests for the manuscript chain",
        "(`docs/RESULTS_EN.md`, `docs/FULL_MANUSCRIPT_EN.md`, ...). Source sequencing",
        "data and internal sample keys are not included.",
        "",
        "## Evidence audit gate",
        "",
        f"- status: `{evidence_status.get('status')}`",
        f"- claims: {evidence_status.get('claim_count')} "
        f"({evidence_status.get('accepted_count')} accepted, "
        f"{evidence_status.get('limited_count')} limited, "
        f"{evidence_status.get('rejected_count')} rejected)",
        "",
        "## Figure QA",
        "",
        f"- chatgpt figure QA: `{figure_qa.get('status')}`",
        f"- checked composites: {len([row for row in figure_qa.get('checks', []) if row.get('figure', '').startswith('Figure')])}",
        "",
        "## Bundled files",
        "",
    ]
    for path in copied_files:
        lines.append(f"- `{path.as_posix()}`")
    lines.extend(
        [
            "",
        "## Related manuscript sources (not copied here)",
        "",
        "- `docs/codex_reanalysis/EVIDENCE_MATRIX.tsv`",
        "- `docs/codex_reanalysis/MANUSCRIPT_VALIDATION.md`",
        "- `results/codex_reanalysis_v1/04b_mhc_phase_blocks/` (WhatsHap tables behind Figure 8;",
        "  keyed by real sample_id and kept untracked)",
        "",
        ]
    )
    return "\n".join(lines)


def package_manuscript(repo_root: Path, output_root: Path, *, force: bool = False) -> dict:
    repo_root = repo_root.resolve()
    output_root = output_root.resolve()
    if output_root.name != "08_manuscript":
        raise ManuscriptPackageError("output root must be results/codex_reanalysis_v1/08_manuscript")

    evidence_status = _validate_evidence_gate(repo_root)
    figure_qa = _validate_figure_qa(repo_root)

    figure_source = repo_root / "docs/chatgpt_figures_submission"
    manifests_source = repo_root / "results/chatgpt_figure_reanalysis"

    copied: list[Path] = []
    input_hashes: dict[str, str] = {}

    figures_dir = output_root / "figures"
    for name in FIGURE_NAMES:
        source = _require_file(figure_source / f"{name}.svg", f"figure SVG {name}")
        destination = figures_dir / f"{name}.svg"
        _copy_file(source, destination, force)
        copied.append(destination.relative_to(repo_root))
        input_hashes[str(source.relative_to(repo_root))] = sha256_file(source)

    for relative_name in ("figure_qa.json", "plot_data.json"):
        source = _require_file(manifests_source / relative_name, relative_name)
        destination = output_root / relative_name
        _copy_file(source, destination, force)
        copied.append(destination.relative_to(repo_root))
        input_hashes[str(source.relative_to(repo_root))] = sha256_file(source)

    export_source = _require_file(manifests_source / "export_manifest.json", "export_manifest.json")
    export_destination = output_root / "export_manifest.json"
    _write_new_text(
        export_destination,
        _sanitize_export_manifest(export_source),
        force=force,
    )
    copied.append(export_destination.relative_to(repo_root))
    input_hashes[str(export_source.relative_to(repo_root))] = sha256_file(export_source)

    inventory_path = output_root / "MANUSCRIPT_PACKAGE.md"
    if inventory_path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite {inventory_path}")
    inventory_path.write_text(
        build_inventory(
            evidence_status=evidence_status,
            figure_qa=figure_qa,
            copied_files=copied,
        ),
        encoding="utf-8",
    )
    copied.append(inventory_path.relative_to(repo_root))

    status_payload = {
        "status": "passed",
        "figure_count": len(FIGURE_NAMES),
        "evidence_audit_status": evidence_status.get("status"),
        "figure_qa_status": figure_qa.get("status"),
        "artifacts": [path.as_posix() for path in copied],
    }
    status_path = output_root / "manuscript_status.json"
    if status_path.exists() and not force:
        raise FileExistsError(f"refusing to overwrite {status_path}")
    status_path.write_text(
        json.dumps(status_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    write_provenance(
        output_root,
        {
            "input_hashes": input_hashes,
            "reference_hashes": {
                "evidence_matrix_status": sha256_file(
                    repo_root / "results/codex_reanalysis_v1/07_evidence_audit/evidence_matrix_status.json"
                ),
            },
        },
    )
    return status_payload


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/codex_reanalysis_v1/08_manuscript"),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing 08_manuscript bundle",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    repo_root = args.repo_root.resolve()
    output_root = args.output if args.output.is_absolute() else repo_root / args.output
    if args.force and output_root.exists():
        shutil.rmtree(output_root)
    payload = package_manuscript(repo_root, output_root, force=args.force)
    print(
        f"Manuscript package written: {output_root} "
        f"({payload['figure_count']} SVG composites, status={payload['status']})."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
