#!/usr/bin/env python3
"""Summarize a PanResistome validation run from generated manifests."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_key_values(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def find_sample_dirs(run_dir: Path) -> list[Path]:
    candidates = []
    for path in run_dir.rglob("panr2_inputs"):
        has_handoff_manifest = (path / "manifest" / "schema_validation_summary.txt").exists()
        has_feature_matrix = (path / "features" / "all_features.tsv").exists()
        if path.is_dir() and (has_handoff_manifest or has_feature_matrix):
            candidates.append(path.parent)
    return sorted(set(candidates))


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t" if path.suffix == ".tsv" else ",")
        rows = list(reader)
    return max(len(rows) - 1, 0)


def summarize_sample(sample_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    manifest = sample_dir / "panr2_inputs" / "manifest"
    features_dir = sample_dir / "panr2_inputs" / "features"

    def add(metric: str, value: object, detail: str = "") -> None:
        rows.append(
            {
                "sample_dir": str(sample_dir),
                "metric": metric,
                "value": str(value),
                "detail": detail,
            }
        )

    metadata_paths = [
        sample_dir / "metadata_output" / "ncbi_enriched.csv",
        sample_dir / "metadata_output" / "ncbi_clean.csv",
        sample_dir / "metadata_output" / "fetchm2_clean.csv",
    ]
    for metadata_path in metadata_paths:
        if metadata_path.exists():
            add("metadata_rows", count_rows(metadata_path), str(metadata_path.relative_to(sample_dir)))
            break

    fasta_count = len(list((sample_dir / "sequence").glob("*.fna"))) + len(
        list((sample_dir / "sequence").glob("*.fa"))
    ) + len(list((sample_dir / "sequence").glob("*.fasta")))
    add("downloaded_fastas", fasta_count)

    qc_rows = read_table(sample_dir / "qc" / "qc_master_report.csv")
    if qc_rows:
        status_counter = Counter(row.get("qc_master_status", row.get("status", "")) or "unknown" for row in qc_rows)
        for status, count in sorted(status_counter.items()):
            add("qc_master_status", count, status)

    schema_values = read_key_values(manifest / "schema_validation_summary.txt")
    for key in [
        "feature_rows",
        "metadata_rows",
        "matched_feature_rows",
        "unmatched_feature_rows",
        "invalid_feature_rows",
        "duplicate_feature_rows",
    ]:
        if key in schema_values:
            add(f"schema_{key}", schema_values[key])

    db_rows = read_table(manifest / "database_setup_status.tsv")
    if db_rows:
        required_failures = [
            row
            for row in db_rows
            if row.get("required_for_profile", "").lower() == "true" and row.get("status") == "FAIL"
        ]
        add("database_setup_required_failures", len(required_failures))
        for row in db_rows:
            add(
                "database_setup_status",
                row.get("status", ""),
                row.get("database_or_tool", ""),
            )

    audit_rows = read_table(manifest / "feature_completeness_audit.tsv")
    if audit_rows:
        for row in audit_rows:
            add(
                "feature_completeness",
                row.get("status", ""),
                f"{row.get('database', '')}: {row.get('feature_rows', '')} rows",
            )

    for feature_table in sorted(features_dir.glob("*.features.tsv")):
        add("feature_table_rows", count_rows(feature_table), feature_table.name)

    all_features = features_dir / "all_features.tsv"
    if all_features.exists():
        add("all_features_rows", count_rows(all_features))

    add("dashboard_exists", (sample_dir / "report" / "index.html").exists())
    return rows


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample_dir", "metric", "value", "detail"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(rows: list[dict[str, str]], path: Path, run_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# PanResistome Validation Summary",
        "",
        f"Run directory: `{run_dir}`",
        "",
        "| Sample directory | Metric | Value | Detail |",
        "| --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| `{row['sample_dir']}` | {row['metric']} | {row['value']} | {row['detail']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="PanResistome run directory or sample directory")
    parser.add_argument("--out-dir", type=Path, help="Directory for validation_summary.csv/md")
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.exists():
        raise SystemExit(f"Run directory does not exist: {run_dir}")

    sample_dirs = find_sample_dirs(run_dir)
    if not sample_dirs and (run_dir / "panr2_inputs").exists():
        sample_dirs = [run_dir]
    if not sample_dirs:
        raise SystemExit(f"No PanR2 handoff directories found under: {run_dir}")

    rows: list[dict[str, str]] = []
    for sample_dir in sample_dirs:
        rows.extend(summarize_sample(sample_dir))

    out_dir = (args.out_dir or run_dir).resolve()
    write_csv(rows, out_dir / "validation_summary.csv")
    write_markdown(rows, out_dir / "validation_summary.md", run_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
