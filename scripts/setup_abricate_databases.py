#!/usr/bin/env python3
"""Set up and optionally refresh ABRicate databases for PanResistome.

PanR2 owns the database bootstrap helper used by PanResistome, but PanResistome
needs an explicit audit trail for remote-user runs. This wrapper records whether
the requested ABRicate databases were present before setup, what setup/update
commands were attempted, and whether each database is available afterward.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "database",
    "requested",
    "present_before",
    "setup_requested",
    "update_requested",
    "setup_status",
    "update_status",
    "present_after",
    "status",
    "message",
]


def run_command(command: list[str], timeout: int = 1800) -> tuple[int, str]:
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


def first_nonempty_line(text: str) -> str:
    for line in str(text or "").splitlines():
        line = line.strip()
        if line:
            return line
    return ""


def available_abricate_databases() -> tuple[set[str], str]:
    if shutil.which("abricate") is None:
        return set(), "ABRicate executable was not found."
    code, output = run_command(["abricate", "--list"], timeout=120)
    if code != 0:
        return set(), output or "abricate --list failed."
    databases: set[str] = set()
    for line in output.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        name = parts[0]
        if name.lower() in {"database", "db", "name"}:
            continue
        databases.add(name)
    return databases, first_nonempty_line(output)


def parse_databases(value: str) -> list[str]:
    seen: set[str] = set()
    databases: list[str] = []
    for item in str(value or "").split(","):
        db = item.strip()
        if db and db not in seen:
            databases.append(db)
            seen.add(db)
    return databases


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def setup_with_panr(databases: list[str], check_only: bool) -> tuple[str, str]:
    if shutil.which("panr") is None:
        if check_only:
            return "SKIPPED", "panr CLI was not found; check-only mode will rely on abricate --list."
        return "FAIL", "panr CLI was not found; cannot run panr setup-db."

    command = ["panr", "setup-db", "--dbs", ",".join(databases)]
    if check_only:
        command.append("--check-only")
    code, output = run_command(command)
    status = "PASS" if code == 0 else "FAIL"
    return status, first_nonempty_line(output) or f"{' '.join(command)} exited {code}."


def update_with_abricate_get_db(databases: list[str]) -> tuple[str, str]:
    updater = shutil.which("abricate-get_db")
    if updater is None:
        return "FAIL", "abricate-get_db was not found; cannot force-refresh ABRicate databases."

    failures: list[str] = []
    messages: list[str] = []
    for db in databases:
        code, output = run_command([updater, "--db", db, "--force"])
        line = first_nonempty_line(output) or f"abricate-get_db --db {db} --force exited {code}."
        messages.append(f"{db}: {line}")
        if code != 0:
            failures.append(db)

    if shutil.which("abricate") is not None:
        code, output = run_command(["abricate", "--setupdb"], timeout=1800)
        line = first_nonempty_line(output) or f"abricate --setupdb exited {code}."
        messages.append(f"setupdb: {line}")
        if code != 0:
            failures.append("abricate --setupdb")

    if failures:
        return "FAIL", "; ".join(messages)
    return "PASS", "; ".join(messages)


def build_rows(args: argparse.Namespace) -> tuple[list[dict[str, str]], int]:
    requested = parse_databases(args.dbs)
    before, before_message = available_abricate_databases()

    setup_status, setup_message = setup_with_panr(requested, args.check_only)
    update_status = "SKIPPED"
    update_message = "ABRicate force-refresh was not requested."
    if args.update:
        update_status, update_message = update_with_abricate_get_db(requested)

    after, after_message = available_abricate_databases()
    rows: list[dict[str, str]] = []
    failures = 0
    for db in requested:
        present_after = db in after
        update_failed = bool(args.update and update_status == "FAIL")
        status = "PASS" if present_after and not update_failed else "FAIL"
        if status == "FAIL":
            failures += 1
        message_parts = [
            f"before: {before_message}",
            f"setup: {setup_message}",
            f"update: {update_message}",
            f"after: {after_message}",
        ]
        rows.append(
            {
                "database": db,
                "requested": "true",
                "present_before": str(db in before).lower(),
                "setup_requested": str(not args.check_only).lower(),
                "update_requested": str(bool(args.update)).lower(),
                "setup_status": setup_status,
                "update_status": update_status,
                "present_after": str(present_after).lower(),
                "status": status,
                "message": " | ".join(message_parts),
            }
        )

    return rows, failures


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dbs", required=True, help="Comma-separated ABRicate database names.")
    parser.add_argument("--out", required=True, help="Output TSV status path.")
    parser.add_argument("--check-only", action="store_true", help="Only check database availability.")
    parser.add_argument("--update", action="store_true", help="Force-refresh requested DBs with abricate-get_db when available.")
    parser.add_argument("--allow-missing", action="store_true", help="Do not exit non-zero if requested databases remain missing.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows, failures = build_rows(args)
    out = Path(args.out)
    write_rows(out, rows)
    if failures and not args.allow_missing:
        missing = ", ".join(row["database"] for row in rows if row["status"] == "FAIL")
        print(f"ABRicate database setup failed for: {missing}; see {out}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
