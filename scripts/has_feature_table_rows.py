#!/usr/bin/env python3
"""Exit 0 when a directory contains at least one meaningful feature table row."""

import csv
import sys
from pathlib import Path


IGNORED_FIELDS = {
    "sample_id",
    "sample",
    "tool",
    "status",
    "source_table",
    "raw_dir",
    "tables_seen",
    "rows_collected",
}


def delimiter_for(path):
    text = path.read_text(errors="ignore")
    first = next((line for line in text.splitlines() if line.strip() and not line.startswith("#")), "")
    return "\t" if first.count("\t") >= first.count(",") else ","


def has_meaningful_row(path):
    try:
        delimiter = delimiter_for(path)
        with path.open(newline="", errors="ignore") as handle:
            lines = [line for line in handle if line.strip() and not line.startswith("#")]
    except OSError:
        return False
    if len(lines) < 2:
        return False
    reader = csv.DictReader(lines, delimiter=delimiter)
    if not reader.fieldnames:
        return False
    for row in reader:
        normalized = {str(key or "").strip().lower(): str(value or "").strip() for key, value in row.items()}
        if normalized.get("rows_collected") == "0":
            continue
        for key, value in normalized.items():
            if key not in IGNORED_FIELDS and value:
                return True
    return False


def main():
    if len(sys.argv) != 2:
        print("usage: has_feature_table_rows.py TABLE_DIR", file=sys.stderr)
        return 2
    table_dir = Path(sys.argv[1])
    if not table_dir.is_dir():
        return 1
    for path in table_dir.rglob("*"):
        if path.is_file() and path.suffix.lower() in {".csv", ".tsv", ".tab"}:
            if has_meaningful_row(path):
                return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
