#!/usr/bin/env python3
"""Fast readiness check for optional geNomad validation.

This script intentionally does not download the geNomad database. It gives a
quick answer before a user launches a long Nextflow run: is the geNomad command
available, does the database path contain files, and what command should be run
next?
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


FIELDS = [
    "genomad_available",
    "genomad_executable",
    "genomad_version",
    "database_requested",
    "database_resolved",
    "database_has_files",
    "status",
    "message",
]


def has_files(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    return any(child.is_file() for child in path.rglob("*"))


def resolve_database_dir(path: Path) -> Path:
    nested = path / "genomad_db"
    if nested.exists() and has_files(nested):
        return nested
    return path


def capture_version(executable: str | None) -> str:
    if not executable:
        return ""
    for command in ([executable, "--version"], [executable, "-h"]):
        try:
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
                timeout=20,
            )
        except Exception as exc:  # pragma: no cover - defensive CLI helper
            return f"unavailable: {exc}"
        output = (completed.stdout or "").strip()
        if completed.returncode == 0 and output:
            return output.splitlines()[0]
    return "unavailable"


def write_report(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-dir", default="", help="Existing or planned geNomad database directory.")
    parser.add_argument("--out", default="", help="Optional TSV report path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    executable = shutil.which("genomad")
    requested = Path(args.db_dir).expanduser() if args.db_dir else Path("")
    resolved = resolve_database_dir(requested) if args.db_dir else Path("")
    database_has_files = has_files(resolved) if args.db_dir else False

    if executable and (not args.db_dir or database_has_files):
        status = "PASS"
        message = "geNomad executable is available and database files are present." if args.db_dir else "geNomad executable is available; no database path was checked."
    elif executable:
        status = "WARNING_DB_MISSING"
        message = "geNomad executable is available, but the database directory is missing or empty. Run `genomad download-database <db_dir>` or pass a populated --genomad_db."
    elif database_has_files:
        status = "WARNING_TOOL_MISSING"
        message = "geNomad database files are present, but the genomad executable is not available in PATH. Prebuild/cache the geNomad environment before running --run_genomad true."
    else:
        status = "FAIL"
        message = "geNomad executable is not available and no populated database directory was found. First-run Conda setup and database download may be slow; prefer a cached environment or container validation."

    row = {
        "genomad_available": str(executable is not None).lower(),
        "genomad_executable": executable or "",
        "genomad_version": capture_version(executable),
        "database_requested": str(requested) if args.db_dir else "",
        "database_resolved": str(resolved) if args.db_dir else "",
        "database_has_files": str(database_has_files).lower(),
        "status": status,
        "message": message,
    }
    if args.out:
        write_report(Path(args.out), row)
    print(f"status={status}")
    print(message)
    if args.out:
        print(f"report={args.out}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
