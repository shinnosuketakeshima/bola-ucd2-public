"""Extract target-BAM reads and assemble one BoLA interval with Flye.

Class I interval NC_037350.1:28,450,000-28,780,000 (NC10 through classical
BOLA). Reads come from the UCD2.0 target BAM, never from
``H:\\Deliver_ADD``. Flye stays out of the mapping environment: version and
prefix are recorded in the per-sample status.

DQ/DR and wide-interval assembly are reported separately in
`docs/codex_reanalysis/ASSEMBLY_DQ_DR_REPORT.md` and
`docs/codex_reanalysis/ASSEMBLY_WIDE_REPORT.md`.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.codex_reanalysis.provenance import write_provenance
from scripts.codex_reanalysis.reference_registry import sha256_file


CLASS_I_REGION = "NC_037350.1:28450000-28780000"
CLASS_I_GENOME_SIZE = "0.35m"
DQ_DR_REGION = "NC_037350.1:25400000-26000000"
DQ_DR_GENOME_SIZE = "0.6m"
WIDE_REGION = "NC_037350.1:25400000-27700000"
WIDE_GENOME_SIZE = "2.3m"

EXTRACT_FLAGS = "0x900"
MIN_OVERLAP = "3000"
FLYE_VERSION = "2.9.6-b1802"

INTERVALS = {
    "class_i": {
        "region": CLASS_I_REGION,
        "genome_size": CLASS_I_GENOME_SIZE,
        "panel_kind": "class_i",
        "paf_name": "classI.paf",
    },
    "dq_dr": {
        "region": DQ_DR_REGION,
        "genome_size": DQ_DR_GENOME_SIZE,
        "panel_kind": "class_ii",
        "paf_name": "dq_dr.paf",
    },
    "wide": {
        "region": WIDE_REGION,
        "genome_size": WIDE_GENOME_SIZE,
        "panel_kind": "wide",
        "paf_name": "wide.paf",
    },
}


def paf_filename(interval_name: str) -> str:
    return INTERVALS[interval_name]["paf_name"]


def build_extract_view_argv(
    bam: str, region: str, samtools: str = "samtools"
) -> list[str]:
    """Primary alignments only; secondaries in the target BAM would duplicate reads."""
    return [samtools, "view", "-b", "-F", EXTRACT_FLAGS, bam, region]


def build_extract_fastq_argv(samtools: str = "samtools") -> list[str]:
    return [samtools, "fastq", "-"]


def build_flye_argv(
    reads: str,
    out_dir: str,
    genome_size: str = CLASS_I_GENOME_SIZE,
    threads: int = 4,
    flye: str = "flye",
) -> list[str]:
    return [
        flye,
        "--nano-raw",
        reads,
        "--keep-haplotypes",
        "--min-overlap",
        MIN_OVERLAP,
        "--genome-size",
        genome_size,
        "--threads",
        str(threads),
        "--out-dir",
        out_dir,
    ]


def build_ipd_map_argv(
    assembly: str,
    panel: str,
    minimap2: str = "minimap2",
    threads: int = 4,
) -> list[str]:
    return [
        minimap2,
        "-x",
        "splice",
        "-p",
        "0.2",
        "-N",
        "30",
        "--secondary=yes",
        "-t",
        str(threads),
        assembly,
        panel,
    ]


def parse_assembly_info(text: str) -> list[dict[str, object]]:
    contigs: list[dict[str, object]] = []
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(f"malformed assembly_info row: {line!r}")
        contigs.append({"name": fields[0], "length": int(fields[1])})
    return contigs


def n50(lengths: Sequence[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def count_fastq_reads_and_bases(path: Path) -> tuple[int, int]:
    n_reads = 0
    n_bases = 0
    with path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                n_reads += 1
                n_bases += len(line.strip())
    return n_reads, n_bases


def _run(argv: Sequence[str], *, stdout=None, stderr=None) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(list(argv), stdout=stdout, stderr=stderr, check=False)


def _tool_version(executable: str) -> str:
    result = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=False
    )
    text = (result.stdout or result.stderr or "").strip().splitlines()
    return text[0] if text else "unknown"


def assemble_one(
    *,
    sample_id: str,
    bam: Path,
    panel: Path,
    output_dir: Path,
    interval_name: str = "class_i",
    samtools: str = "samtools",
    flye: str = "flye",
    minimap2: str = "minimap2",
    threads: int = 4,
) -> dict:
    if interval_name not in INTERVALS:
        raise ValueError(f"unknown interval: {interval_name}")
    spec = INTERVALS[interval_name]
    region = spec["region"]
    genome_size = spec["genome_size"]
    output_dir.mkdir(parents=True, exist_ok=True)

    reads = output_dir / "reads.fastq"
    flye_dir = output_dir / "flye"
    flye_log = output_dir / "flye.log"
    if flye_dir.exists():
        shutil.rmtree(flye_dir)

    view = subprocess.Popen(
        build_extract_view_argv(str(bam), region, samtools=samtools),
        stdout=subprocess.PIPE,
    )
    with reads.open("wb") as handle:
        fastq = subprocess.Popen(
            build_extract_fastq_argv(samtools=samtools),
            stdin=view.stdout,
            stdout=handle,
        )
        if view.stdout is not None:
            view.stdout.close()
        if fastq.wait() != 0 or view.wait() != 0:
            raise RuntimeError(
                f"extract failed for {sample_id}: view={view.returncode} fastq={fastq.returncode}"
            )

    n_reads, n_bases = count_fastq_reads_and_bases(reads)
    contig, start_s, end_s = region.split(":")[0], *region.split(":")[1].split("-")
    region_span = int(end_s) - int(start_s) + 1
    estimated_coverage = (n_bases / region_span) if region_span else 0.0

    flye_argv = build_flye_argv(
        str(reads), str(flye_dir), genome_size=genome_size, threads=threads, flye=flye
    )
    with flye_log.open("wb") as log_handle:
        flye_proc = _run(flye_argv, stdout=log_handle, stderr=subprocess.STDOUT)

    assembly_src = flye_dir / "assembly.fasta"
    info_src = flye_dir / "assembly_info.txt"
    assembly_dst = output_dir / "assembly.fasta"
    info_dst = output_dir / "assembly_info.txt"
    paf_dst = output_dir / paf_filename(interval_name)

    passed = flye_proc.returncode == 0 and assembly_src.is_file() and info_src.is_file()
    contig_records: list[dict[str, object]] = []
    if passed:
        shutil.copy2(assembly_src, assembly_dst)
        shutil.copy2(info_src, info_dst)
        contig_records = parse_assembly_info(info_dst.read_text(encoding="utf-8"))
        map_argv = build_ipd_map_argv(
            str(assembly_dst), str(panel), minimap2=minimap2, threads=threads
        )
        with paf_dst.open("wb") as paf_handle, (output_dir / "minimap2.log").open("wb") as map_log:
            map_proc = _run(map_argv, stdout=paf_handle, stderr=map_log)
        if map_proc.returncode != 0:
            passed = False
            paf_dst.write_text("", encoding="utf-8")
    else:
        if not paf_dst.exists():
            paf_dst.write_text("", encoding="utf-8")

    lengths = [int(row["length"]) for row in contig_records]
    flye_version = _tool_version(flye)
    status = {
        "sample_id": sample_id,
        "interval_name": interval_name,
        "region": region,
        "contig": contig,
        "n_reads": n_reads,
        "n_bases": n_bases,
        "estimated_coverage": round(estimated_coverage, 2),
        "flye_version": flye_version,
        "flye_returncode": flye_proc.returncode,
        "n_contigs": len(lengths),
        "total_assembly_length": sum(lengths),
        "n50": n50(lengths),
        "longest_contig": max(lengths) if lengths else 0,
        "bam_sha256": sha256_file(bam),
        "panel_sha256": sha256_file(panel),
        "assembly_sha256": sha256_file(assembly_dst) if assembly_dst.is_file() else "",
        "status": "passed" if passed else "failed",
    }
    (output_dir / f"{sample_id}.assembly_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        output_dir,
        {
            "versions": {
                "python": sys.version.split()[0],
                "flye": flye_version,
                "samtools": _tool_version(samtools),
                "minimap2": _tool_version(minimap2),
            },
            "input_hashes": {
                "target_bam": status["bam_sha256"],
                "ipd_panel": status["panel_sha256"],
            },
            "reference_hashes": {
                "assembly": status["assembly_sha256"],
            },
        },
    )
    reads.unlink(missing_ok=True)
    if flye_dir.exists():
        shutil.rmtree(flye_dir)
    return status


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--interval", default="class_i", choices=sorted(INTERVALS))
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--flye", default="flye")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    status = assemble_one(
        sample_id=args.sample_id,
        bam=args.bam,
        panel=args.panel,
        output_dir=args.output_dir,
        interval_name=args.interval,
        samtools=args.samtools,
        flye=args.flye,
        minimap2=args.minimap2,
        threads=args.threads,
    )
    print(
        f"{args.sample_id} {status['status']} reads={status['n_reads']} "
        f"contigs={status['n_contigs']} n50={status['n50']}"
    )
    return 0 if status["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
