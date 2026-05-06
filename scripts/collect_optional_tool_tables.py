#!/usr/bin/env python3
"""Collect heterogeneous optional-tool tables into one PanR2-readable TSV."""

import argparse
import csv
from pathlib import Path


def detect_delimiter(path):
    text = path.read_text(errors="ignore")
    first = next((line for line in text.splitlines() if line.strip() and not line.startswith("#")), "")
    return "\t" if first.count("\t") >= first.count(",") else ","


def read_rows(path):
    delimiter = detect_delimiter(path)
    rows = []
    with path.open(newline="", errors="ignore") as handle:
        filtered = [line for line in handle if line.strip() and not line.startswith("#")]
    if not filtered:
        return rows
    reader = csv.DictReader(filtered, delimiter=delimiter)
    if not reader.fieldnames:
        return rows
    for row in reader:
        rows.append({str(k).strip(): str(v).strip() for k, v in row.items() if k is not None})
    return rows


def sample_from_path(path, raw_dir):
    rel = path.relative_to(raw_dir)
    if len(rel.parts) > 1:
        return rel.parts[0]
    name = path.stem
    for suffix in ["_results", "_result", "_summary", "_report", "_typer", "_mobtyper"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def main():
    parser = argparse.ArgumentParser(description="Collect optional bioinformatics tool tables for PanR2.")
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--tool", required=True)
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    status_path = out.parent / f"{args.tool}_collection_status.tsv"

    table_paths = []
    if raw_dir.exists():
        for path in raw_dir.rglob("*"):
            if path.is_file() and path.suffix.lower() in {".csv", ".tsv", ".tab", ".txt"}:
                table_paths.append(path)

    collected = []
    fieldnames = {"sample_id", "tool"}
    for path in sorted(table_paths):
        for row in read_rows(path):
            if not any(str(value).strip() for value in row.values()):
                continue
            row = dict(row)
            row.setdefault("sample_id", sample_from_path(path, raw_dir))
            row.setdefault("tool", args.tool)
            row["source_table"] = str(path)
            fieldnames.update(row.keys())
            collected.append(row)

    ordered = ["sample_id", "tool"] + sorted(name for name in fieldnames if name not in {"sample_id", "tool"})
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(collected)

    with status_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["tool", "raw_dir", "tables_seen", "rows_collected"], delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "tool": args.tool,
            "raw_dir": str(raw_dir),
            "tables_seen": len(table_paths),
            "rows_collected": len(collected),
        })


if __name__ == "__main__":
    main()
