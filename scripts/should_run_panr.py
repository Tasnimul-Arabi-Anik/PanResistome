#!/usr/bin/env python3
"""Decide whether a PanResistome sample directory has enough input for PanR2."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    sep = "\t" if path.suffix.lower() in {".tab", ".tsv"} else ","
    try:
        return pd.read_csv(path, sep=sep)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _first_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def should_run(sample_dir: Path) -> tuple[bool, str]:
    metadata = _read_table(sample_dir / "metadata_output" / "ncbi_clean.csv")
    if metadata.empty:
        return False, "metadata_output/ncbi_clean.csv has no rows after filtering"

    summary_path = _first_existing([
        sample_dir / "abricate" / "ncbi_summary.tab",
        sample_dir / "abricate" / "ncbi_summary.csv",
    ])
    results_path = _first_existing([
        sample_dir / "abricate" / "ncbi_results.tab",
        sample_dir / "abricate" / "ncbi_results.csv",
    ])
    if summary_path is None:
        return False, "ABRicate summary file is missing"
    if results_path is None:
        return False, "ABRicate results file is missing"

    summary = _read_table(summary_path)
    results = _read_table(results_path)
    if summary.empty or "#FILE" not in {str(col).upper() for col in summary.columns}:
        return False, "ABRicate summary has no parseable #FILE header"
    if results.empty:
        return False, "ABRicate results table is empty; no database hits were available for PanR2"
    required_results = {"GENE", "RESISTANCE"}
    upper_results = {str(col).upper() for col in results.columns}
    missing = sorted(required_results - upper_results)
    if missing:
        return False, f"ABRicate results are missing required columns: {','.join(missing)}"
    return True, "PanR2 inputs are present"


def main() -> None:
    parser = argparse.ArgumentParser(description="Check whether PanR2 should run for a PanResistome sample directory.")
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--reason-file", type=Path)
    args = parser.parse_args()
    ok, reason = should_run(args.sample_dir)
    if args.reason_file:
        args.reason_file.parent.mkdir(parents=True, exist_ok=True)
        args.reason_file.write_text(reason + "\n", encoding="utf-8")
    print(reason)
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
