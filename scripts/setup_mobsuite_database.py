#!/usr/bin/env python3
"""Prepare and audit a MOB-suite database directory."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "database_dir",
    "auto_init_requested",
    "auto_init_taxa_requested",
    "mob_init_status",
    "taxa_init_status",
    "core_status",
    "taxa_status",
    "status",
    "message",
]

REQUIRED_FILES = [
    "clusters.txt",
    "host_range_literature_plasmidDB.txt",
    "mob.proteins.faa",
    "mpf.proteins.faa",
    "ncbi_plasmid_full_seqs.fas",
    "ncbi_plasmid_full_seqs.fas.msh",
    "orit.fas",
    "rep.dna.fas",
    "repetitive.dna.fas",
]

REQUIRED_PREFIXES = [
    "ncbi_plasmid_full_seqs.fas.n",
    "repetitive.dna.fas.n",
]


def as_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "on"}


def run_command(command: list[str], timeout: int = 7200) -> tuple[int, str]:
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return 127, str(exc)
    return completed.returncode, completed.stdout.strip()


def first_line(text: str) -> str:
    for line in str(text or "").splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def core_status(path: Path) -> tuple[str, str]:
    if not path.exists():
        return "FAIL", "database directory does not exist"
    if not path.is_dir():
        return "FAIL", "database path is not a directory"
    missing_files = [name for name in REQUIRED_FILES if not (path / name).exists()]
    missing_prefixes = [
        prefix for prefix in REQUIRED_PREFIXES if not any(child.name.startswith(prefix) for child in path.iterdir())
    ]
    if missing_files or missing_prefixes:
        missing = missing_files + [f"{prefix}*" for prefix in missing_prefixes]
        return "FAIL", "missing " + ", ".join(missing)
    return "PASS", "core MOB-suite database files are present"


def taxa_status(path: Path) -> tuple[str, str]:
    taxa = path / "taxa.sqlite"
    if taxa.exists() and taxa.stat().st_size > 0:
        return "PASS", "taxa.sqlite is present"
    return "WARNING_TAXA_MISSING", "taxa.sqlite is missing; mob_recon may initialize ETE taxonomy at runtime"


def run_mob_init(path: Path) -> tuple[str, str]:
    mob_init = shutil.which("mob_init")
    if not mob_init:
        return "FAIL", "mob_init executable was not found"
    path.mkdir(parents=True, exist_ok=True)
    attempts = [
        [mob_init, "-d", str(path)],
        [mob_init, "--database_directory", str(path)],
    ]
    messages: list[str] = []
    for command in attempts:
        code, output = run_command(command)
        messages.append(first_line(output) or f"{' '.join(command)} exited {code}")
        if code == 0:
            return "PASS", "; ".join(messages)
    return "FAIL", "; ".join(messages)


def run_taxa_init(path: Path) -> tuple[str, str]:
    taxa_path = path / "taxa.sqlite"
    code, output = run_command(
        [
            sys.executable,
            "-c",
            (
                "from ete3 import NCBITaxa; "
                f"NCBITaxa(dbfile={str(taxa_path)!r}); "
                f"print({str(taxa_path)!r})"
            ),
        ],
        timeout=7200,
    )
    status, message = taxa_status(path)
    if code == 0 and status == "PASS":
        return "PASS", first_line(output) or "ETE taxonomy database initialized"
    return "FAIL", first_line(output) or message


def write_row(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--auto-init", type=as_bool, default=True)
    parser.add_argument("--auto-init-taxa", type=as_bool, default=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    db_dir = Path(args.db_dir)
    before_core, before_core_message = core_status(db_dir)

    mob_init_status = "SKIPPED"
    mob_init_message = "core database files were already present"
    if before_core != "PASS" and args.auto_init:
        mob_init_status, mob_init_message = run_mob_init(db_dir)

    after_core, after_core_message = core_status(db_dir)
    before_taxa, before_taxa_message = taxa_status(db_dir)
    taxa_init_status = "SKIPPED"
    taxa_init_message = before_taxa_message
    if after_core == "PASS" and before_taxa != "PASS" and args.auto_init_taxa:
        taxa_init_status, taxa_init_message = run_taxa_init(db_dir)

    after_taxa, after_taxa_message = taxa_status(db_dir)
    if after_core != "PASS":
        status = "FAIL"
    elif after_taxa != "PASS":
        status = "WARNING_TAXA_MISSING"
    else:
        status = "PASS"

    message = (
        f"core before: {before_core_message}; "
        f"mob_init: {mob_init_message}; "
        f"core after: {after_core_message}; "
        f"taxa: {taxa_init_message or after_taxa_message}"
    )
    row = {
        "database_dir": str(db_dir),
        "auto_init_requested": str(bool(args.auto_init)).lower(),
        "auto_init_taxa_requested": str(bool(args.auto_init_taxa)).lower(),
        "mob_init_status": mob_init_status,
        "taxa_init_status": taxa_init_status,
        "core_status": after_core,
        "taxa_status": after_taxa,
        "status": status,
        "message": message,
    }
    write_row(Path(args.out), row)
    if args.strict and status == "FAIL":
        print(f"MOB-suite database setup failed; see {args.out}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
