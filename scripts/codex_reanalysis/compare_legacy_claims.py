"""Compare the evidence matrix against the legacy manuscript baseline.

Baseline is `docs/INTEGRATED_MANUSCRIPT_EN.docx` per the plan's Global Constraints
("The original docs/INTEGRATED_MANUSCRIPT_EN.docx is immutable and remains the
comparison baseline"), not `docs/CHATGPT_EVIDENCE_AUDIT.md` (a separate manuscript
lineage -- see `build_evidence_matrix.py`'s module docstring for why).

For every matrix row with a `legacy_excerpt`, this verifies the excerpt is actually
present in the docx text -- so `LEGACY_COMPARISON.md` reports a checked fact, not an
assumed one. A row without an excerpt is reported as having no legacy counterpart.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.build_evidence_matrix import SCOPE_LIMITATIONS

_TAG_RE = re.compile(r"<[^>]+>")
_WHITESPACE_RE = re.compile(r"\s+")


class LegacyComparisonError(ValueError):
    """Raised when a docx cannot be read or an excerpt cannot be located."""


def extract_docx_text(path: Path) -> str:
    """Return the flattened body text of a .docx file's main document part."""
    try:
        with zipfile.ZipFile(path) as archive:
            xml = archive.read("word/document.xml").decode("utf-8", errors="ignore")
    except (zipfile.BadZipFile, KeyError) as exc:
        raise LegacyComparisonError(f"could not read word/document.xml from {path}: {exc}") from exc
    return normalize_whitespace(_TAG_RE.sub("", xml))


def normalize_whitespace(text: str) -> str:
    return _WHITESPACE_RE.sub(" ", text).strip()


def excerpt_found(docx_text: str, excerpt: str) -> bool:
    return normalize_whitespace(excerpt) in docx_text


def read_matrix_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def render_comparison(rows: Sequence[dict[str, str]]) -> str:
    lines = [
        "# Legacy Baseline Comparison",
        "",
        "Baseline: `docs/INTEGRATED_MANUSCRIPT_EN.docx` (immutable, per the reanalysis "
        "plan's Global Constraints). Every row below is drawn from "
        "`docs/codex_reanalysis/EVIDENCE_MATRIX.tsv`; rows citing a legacy excerpt were "
        "checked against the docx text directly, not assumed.",
        "",
        "| claim_id | classification | legacy value | new value | action |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['claim_id']} | {row['classification']} | {row['legacy_value']} | "
            f"{row['new_value']} | {row['action']} |"
        )
    lines.append("")
    lines.append(f"{len(rows)} claims compared; every declared legacy excerpt was verified present in the docx text.")
    if SCOPE_LIMITATIONS:
        lines.append("")
        lines.append(
            "**Not yet covered by this seed** (absent here means \"not yet audited,\" never "
            "\"no finding\" -- see `build_evidence_matrix.py`'s `SCOPE_LIMITATIONS`):"
        )
        for key, reason in SCOPE_LIMITATIONS.items():
            lines.append(f"- `{key}`: {reason}")
    return "\n".join(lines) + "\n"


def verify_excerpts(rows: Sequence[dict[str, str]], docx_text: str) -> list[str]:
    """Return claim_ids whose declared legacy excerpt is not found in the docx text.

    A row without an excerpt (empty string, written from a claim whose legacy_excerpt
    was None) is not checked -- it declares no legacy counterpart to verify.
    """
    unverified = []
    for row in rows:
        excerpt = row.get("legacy_excerpt", "")
        if excerpt and not excerpt_found(docx_text, excerpt):
            unverified.append(row["claim_id"])
    return unverified


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy", type=Path, required=True, help="docs/INTEGRATED_MANUSCRIPT_EN.docx")
    parser.add_argument("--matrix", type=Path, required=True, help="docs/codex_reanalysis/EVIDENCE_MATRIX.tsv")
    parser.add_argument("--output", type=Path, required=True, help="docs/codex_reanalysis/LEGACY_COMPARISON.md")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    docx_text = extract_docx_text(args.legacy)
    rows = read_matrix_tsv(args.matrix)

    unverified = verify_excerpts(rows, docx_text)
    if unverified:
        raise LegacyComparisonError(
            f"legacy excerpt not found in {args.legacy} for claim_id(s): {unverified}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_comparison(rows), encoding="utf-8")
    print(f"Legacy comparison written: {args.output} ({len(rows)} claims, 0 unverified excerpts)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
