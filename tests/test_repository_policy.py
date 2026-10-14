from __future__ import annotations

import subprocess
from pathlib import Path


MAX_TRACKED_BYTES = 20 * 1024 * 1024
FORBIDDEN_SUFFIXES = {
    ".bam", ".bai", ".bcf", ".cram", ".crai", ".fast5", ".fasta",
    ".fastq", ".fa", ".fna", ".fq", ".gvcf", ".mmi", ".pod5", ".sam",
    ".vcf",
}
SOURCE_DATA_ROOT = Path("data/raw").resolve()
DELIVERY_ROOT = Path(r"H:\Deliver_ADD").resolve()
PROTECTED_ROOTS = (SOURCE_DATA_ROOT, DELIVERY_ROOT)


def is_forbidden(path: Path) -> bool:
    """Return whether a possible Git candidate violates repository policy."""
    resolved = path.resolve()
    if any(resolved.is_relative_to(root) for root in PROTECTED_ROOTS):
        return True
    if any(suffix.lower() in FORBIDDEN_SUFFIXES for suffix in path.suffixes):
        return True
    return path.exists() and path.is_file() and path.stat().st_size > MAX_TRACKED_BYTES


def tracked_candidates(root: Path) -> list[Path]:
    """List files actually staged in Git's index for a commit."""
    output = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--cached", "-z"],
        check=True,
        capture_output=True,
    ).stdout
    return [root / Path(name) for name in output.decode().split("\0") if name]


def test_repository_policy_rejects_genomic_extensions(tmp_path):
    bad = tmp_path / "sample.bam"
    bad.write_bytes(b"BAM")
    assert is_forbidden(bad)


def test_repository_policy_allows_source_and_small_evidence_files(tmp_path):
    assert not is_forbidden(tmp_path / "script.py")
    assert not is_forbidden(tmp_path / "evidence.tsv")


def test_repository_policy_rejects_large_files(tmp_path):
    large = tmp_path / "large.txt"
    large.write_bytes(b"0" * (20 * 1024 * 1024 + 1))
    assert is_forbidden(large)


def test_repository_policy_rejects_a_link_resolving_into_delivery_root(tmp_path):
    link = tmp_path / "delivery_root_link"
    link.symlink_to(DELIVERY_ROOT, target_is_directory=True)
    assert is_forbidden(link / "sample.txt")


def test_repository_policy_rejects_a_delivery_root_path():
    assert is_forbidden(DELIVERY_ROOT / "sample.txt")


def test_workspace_has_no_source_data_symlink_in_tracked_candidates():
    candidates = tracked_candidates(Path("."))
    assert all(not p.is_symlink() and not is_forbidden(p) for p in candidates)
