#!/usr/bin/env python3
"""Run MOB-suite mob_recon per genome with bounded parallelism."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


STATUS_FIELDS = [
    "sample_id",
    "fasta_path",
    "output_dir",
    "status",
    "exit_code",
    "runtime_seconds",
    "message",
]


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def sample_id(path: Path) -> str:
    return path.stem


def has_reusable_output(out_dir: Path) -> bool:
    if not out_dir.exists():
        return False
    for path in out_dir.rglob("*"):
        if not path.is_file() or path.stat().st_size == 0:
            continue
        if path.name in {"mob_recon.stdout", "mob_recon.stderr"}:
            continue
        if path.suffix.lower() in {".tsv", ".csv", ".txt", ".tab"}:
            return True
    return False


def run_one(fasta: Path, raw_dir: Path, database_dir: str, threads: int, reuse_existing: bool) -> dict[str, str]:
    sid = sample_id(fasta)
    out_dir = raw_dir / sid
    out_dir.mkdir(parents=True, exist_ok=True)
    if reuse_existing and has_reusable_output(out_dir):
        return {
            "sample_id": sid,
            "fasta_path": str(fasta),
            "output_dir": str(out_dir),
            "status": "REUSED",
            "exit_code": "0",
            "runtime_seconds": "0.00",
            "message": "Existing MOB-suite output reused.",
        }

    cmd = [
        "mob_recon",
        "--infile",
        str(fasta),
        "--outdir",
        str(out_dir),
        "--num_threads",
        str(threads),
        "--force",
    ]
    if database_dir:
        cmd.extend(["--database_directory", database_dir])

    started = time.monotonic()
    stdout = out_dir / "mob_recon.stdout"
    stderr = out_dir / "mob_recon.stderr"
    with stdout.open("w", encoding="utf-8") as out, stderr.open("w", encoding="utf-8") as err:
        proc = subprocess.run(cmd, stdout=out, stderr=err, check=False)
    runtime = time.monotonic() - started
    status = "PASS" if proc.returncode == 0 else "FAIL"
    return {
        "sample_id": sid,
        "fasta_path": str(fasta),
        "output_dir": str(out_dir),
        "status": status,
        "exit_code": str(proc.returncode),
        "runtime_seconds": f"{runtime:.2f}",
        "message": "mob_recon completed." if status == "PASS" else f"mob_recon exited with {proc.returncode}.",
    }


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-dir", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--status-file", required=True)
    parser.add_argument("--database-directory", default="")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads-per-sample", type=int, default=1)
    parser.add_argument("--reuse-existing", type=as_bool, default=True)
    args = parser.parse_args()

    sequence_dir = Path(args.sequence_dir)
    raw_dir = Path(args.raw_dir)
    status_file = Path(args.status_file)
    raw_dir.mkdir(parents=True, exist_ok=True)
    fasta_files = sorted(sequence_dir.glob("*.fna")) if sequence_dir.exists() else []

    rows: list[dict[str, str]] = []
    if not fasta_files:
        write_status(status_file, rows)
        return
    if subprocess.run(["bash", "-lc", "command -v mob_recon"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
        rows = [
            {
                "sample_id": sample_id(fasta),
                "fasta_path": str(fasta),
                "output_dir": str(raw_dir / sample_id(fasta)),
                "status": "FAIL",
                "exit_code": "127",
                "runtime_seconds": "0.00",
                "message": "mob_recon executable was not found.",
            }
            for fasta in fasta_files
        ]
        write_status(status_file, rows)
        raise SystemExit(1)

    jobs = max(1, int(args.jobs))
    threads = max(1, int(args.threads_per_sample))
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(run_one, fasta, raw_dir, args.database_directory, threads, args.reuse_existing)
            for fasta in fasta_files
        ]
        for future in as_completed(futures):
            rows.append(future.result())
    rows.sort(key=lambda row: row["sample_id"])
    write_status(status_file, rows)
    if any(row["status"] == "FAIL" for row in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
