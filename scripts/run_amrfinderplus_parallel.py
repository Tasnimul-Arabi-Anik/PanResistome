#!/usr/bin/env python3
"""Run AMRFinderPlus per assembly with bounded parallelism."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def run_one(
    fasta: Path,
    raw_dir: Path,
    threads_per_sample: int,
    organism: str | None,
    reuse_existing: bool,
) -> dict[str, str]:
    prefix = fasta.stem
    out = raw_dir / f"{prefix}.tsv"
    log = raw_dir / f"{prefix}.log"
    if reuse_existing and out.exists() and out.stat().st_size > 0:
        return {
            "sample_id": prefix,
            "fasta_path": str(fasta),
            "output_path": str(out),
            "status": "PASS",
            "exit_code": "0",
            "duration_seconds": "0.000",
            "message": "existing_output_reused",
        }

    cmd = [
        "amrfinder",
        "-n",
        str(fasta),
        "--threads",
        str(threads_per_sample),
        "-o",
        str(out),
    ]
    if organism:
        cmd.extend(["--organism", organism])

    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    duration = time.monotonic() - start

    status = "PASS" if completed.returncode == 0 else "FAIL"
    message = "completed" if status == "PASS" else f"see {log}"
    return {
        "sample_id": prefix,
        "fasta_path": str(fasta),
        "output_path": str(out),
        "status": status,
        "exit_code": str(completed.returncode),
        "duration_seconds": f"{duration:.3f}",
        "message": message,
    }


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "fasta_path", "output_path", "status", "exit_code", "duration_seconds", "message"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run AMRFinderPlus per FASTA with bounded parallel workers.")
    parser.add_argument("--sequence-dir", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--status-file", required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads-per-sample", type=int, default=1)
    parser.add_argument("--organism")
    parser.add_argument("--update-db", default="true")
    parser.add_argument(
        "--reuse-existing",
        default="true",
        help="Reuse non-empty per-sample TSVs instead of rerunning AMRFinderPlus.",
    )
    parser.add_argument("--progress-every", type=int, default=10, help="Print progress every N completed samples.")
    args = parser.parse_args()

    sequence_dir = Path(args.sequence_dir)
    raw_dir = Path(args.raw_dir)
    status_file = Path(args.status_file)
    raw_dir.mkdir(parents=True, exist_ok=True)

    amrfinder = shutil.which("amrfinder")
    if not amrfinder:
        write_status(status_file, [{
            "sample_id": "all",
            "fasta_path": str(sequence_dir),
            "output_path": "",
            "status": "FAIL",
            "exit_code": "127",
            "duration_seconds": "0.000",
            "message": "AMRFinderPlus executable missing",
        }])
        return 127

    fasta_files = sorted(sequence_dir.glob("*.fna"))
    if not fasta_files:
        write_status(status_file, [{
            "sample_id": "all",
            "fasta_path": str(sequence_dir),
            "output_path": "",
            "status": "FAIL",
            "exit_code": "2",
            "duration_seconds": "0.000",
            "message": "sequence directory empty",
        }])
        return 2

    if as_bool(args.update_db):
        update_log = raw_dir / "amrfinder_update.log"
        with update_log.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(["amrfinder", "-u"], stdout=handle, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            write_status(status_file, [{
                "sample_id": "all",
                "fasta_path": str(sequence_dir),
                "output_path": str(update_log),
                "status": "FAIL",
                "exit_code": str(completed.returncode),
                "duration_seconds": "0.000",
                "message": f"AMRFinderPlus database update failed; see {update_log}",
            }])
            return completed.returncode

    jobs = max(1, min(args.jobs, len(fasta_files)))
    threads_per_sample = max(1, args.threads_per_sample)
    reuse_existing = as_bool(args.reuse_existing)
    progress_every = max(1, args.progress_every)
    print(
        f"AMRFinderPlus sample run: {len(fasta_files)} FASTA(s), jobs={jobs}, "
        f"threads_per_sample={threads_per_sample}, reuse_existing={str(reuse_existing).lower()}",
        flush=True,
    )
    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(run_one, fasta, raw_dir, threads_per_sample, args.organism, reuse_existing): fasta
            for fasta in fasta_files
        }
        for future in as_completed(futures):
            rows.append(future.result())
            completed_count = len(rows)
            if completed_count == len(fasta_files) or completed_count % progress_every == 0:
                passed_so_far = sum(1 for row in rows if row["status"] == "PASS")
                reused_so_far = sum(1 for row in rows if row.get("message") == "existing_output_reused")
                failed_so_far = sum(1 for row in rows if row["status"] == "FAIL")
                print(
                    f"AMRFinderPlus progress: {completed_count}/{len(fasta_files)} completed "
                    f"({passed_so_far} PASS, {failed_so_far} FAIL, {reused_so_far} reused)",
                    flush=True,
                )

    rows.sort(key=lambda row: row["sample_id"])
    write_status(status_file, rows)

    passed = sum(1 for row in rows if row["status"] == "PASS")
    failed = sum(1 for row in rows if row["status"] == "FAIL")
    if passed == 0 and failed > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
