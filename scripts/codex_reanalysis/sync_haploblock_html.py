"""Embed haploblock_figdata.json into self-contained HTML figure files."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

TARGETS = (
    ("docs/figures/haploblock_ld.html", "const D="),
    ("docs/figures/paper_figures.html", "const HAP="),
)


def sync_file(html_path: Path, marker: str, payload_text: str) -> None:
    text = html_path.read_text(encoding="utf-8")
    pattern = re.compile(re.escape(marker) + r"\{.*?\};", re.DOTALL)
    replacement = f"{marker}{payload_text};"
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise ValueError(f"expected one {marker} assignment in {html_path}, found {count}")
    html_path.write_text(updated, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument(
        "--data",
        type=Path,
        default=Path("docs/figures/haploblock_figdata.json"),
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    data_path = args.data if args.data.is_absolute() else repo_root / args.data
    payload_text = data_path.read_text(encoding="utf-8").strip()
    json.loads(payload_text)

    for relative, marker in TARGETS:
        sync_file(repo_root / relative, marker, payload_text)
        print(f"Updated {relative}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
