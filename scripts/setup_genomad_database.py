#!/usr/bin/env python3
"""Prepare and audit a geNomad database directory."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "requested_database_dir",
    "resolved_database_dir",
    "auto_download_requested",
    "genomad_available",
    "download_status",
    "status",
    "message",
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


def has_files(path: Path) -> bool:
    if not path.exists() or not path.is_dir():
        return False
    return any(child.is_file() for child in path.rglob("*"))


def resolve_database_dir(path: Path) -> Path:
    nested = path / "genomad_db"
    if nested.exists() and has_files(nested):
        return nested
    return path


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
    parser.add_argument("--auto-download", type=as_bool, default=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    requested = Path(args.db_dir)
    genomad = shutil.which("genomad")
    download_status = "SKIPPED"
    message_parts: list[str] = []

    resolved = resolve_database_dir(requested)
    if not has_files(resolved) and args.auto_download:
        if genomad is None:
            download_status = "FAIL"
            message_parts.append("genomad executable was not found; cannot download database")
        else:
            requested.mkdir(parents=True, exist_ok=True)
            code, output = run_command([genomad, "download-database", str(requested)])
            download_status = "PASS" if code == 0 else "FAIL"
            message_parts.append(first_line(output) or f"genomad download-database exited {code}")
            resolved = resolve_database_dir(requested)
    elif has_files(resolved):
        message_parts.append("geNomad database files are present")
    else:
        message_parts.append("geNomad database is missing and auto-download is disabled")

    if has_files(resolved):
        status = "PASS"
    else:
        status = "FAIL"

    row = {
        "requested_database_dir": str(requested),
        "resolved_database_dir": str(resolved),
        "auto_download_requested": str(bool(args.auto_download)).lower(),
        "genomad_available": str(genomad is not None).lower(),
        "download_status": download_status,
        "status": status,
        "message": "; ".join(message_parts),
    }
    write_row(Path(args.out), row)
    if args.strict and status == "FAIL":
        print(f"geNomad database setup failed; see {args.out}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
