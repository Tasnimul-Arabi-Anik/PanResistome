#!/usr/bin/env python3
"""Check core outputs from a PanResistome comprehensive validation run."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ZERO_SCHEMA_KEYS = [
    "unmatched_feature_rows",
    "invalid_feature_rows",
    "duplicate_feature_rows",
]


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


def data_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def find_sample_dirs(run_dir: Path) -> list[Path]:
    if (run_dir / "panr2_inputs").exists():
        return [run_dir]
    candidates = []
    for path in run_dir.rglob("panr2_inputs"):
        manifest = path / "manifest"
        features = path / "features"
        if (manifest / "schema_validation_summary.txt").exists() or (features / "all_features.tsv").exists():
            candidates.append(path.parent)
    return sorted(set(candidates))


def is_failure_status(status: str, allow_warning_failed: bool = False) -> bool:
    normalized = (status or "").strip().upper()
    if normalized in {"", "PASS", "WARNING", "WARNING_EMPTY", "SKIPPED", "OPTIONAL"}:
        return False
    if allow_warning_failed and normalized == "WARNING_FAILED":
        return False
    return "FAIL" in normalized or normalized == "ERROR"


def check_sample_dir(
    sample_dir: Path,
    *,
    require_checkm2: bool,
    require_genomad: bool,
    allow_mobileelementfinder_warning: bool,
    expect_zero_schema_errors: bool,
) -> list[str]:
    failures: list[str] = []
    manifest = sample_dir / "panr2_inputs" / "manifest"
    features = sample_dir / "panr2_inputs" / "features"

    required_paths = [
        sample_dir / "report" / "index.html",
        features / "all_features.tsv",
        manifest / "schema_validation_summary.txt",
        manifest / "database_setup_status.tsv",
        manifest / "native_runner_merge_audit.tsv",
    ]
    for path in required_paths:
        if not path.exists():
            failures.append(f"missing required output: {path.relative_to(sample_dir)}")

    schema = read_key_values(manifest / "schema_validation_summary.txt")
    if expect_zero_schema_errors:
        for key in ZERO_SCHEMA_KEYS:
            if schema.get(key) not in {"0", "0.0"}:
                failures.append(f"schema {key} expected 0, observed {schema.get(key, 'missing')}")

    db_rows = read_table(manifest / "database_setup_status.tsv")
    for row in db_rows:
        required = row.get("required_for_profile", "").strip().lower() == "true"
        status = row.get("status", "")
        if required and is_failure_status(status):
            failures.append(
                f"required database/tool failed: {row.get('database_or_tool', 'unknown')} status={status}"
            )

    audit_rows = read_table(manifest / "native_runner_merge_audit.tsv")
    modules = {row.get("module", "") for row in audit_rows}
    for module in ["abricate", "integronfinder", "mlst"]:
        if module not in modules:
            failures.append(f"native runner audit missing module: {module}")
    for row in audit_rows:
        module = row.get("module", "")
        allow_warning = allow_mobileelementfinder_warning and module == "mobileelementfinder"
        if is_failure_status(row.get("status", ""), allow_warning_failed=allow_warning):
            failures.append(f"native runner audit failure: {module} status={row.get('status', '')}")

    if require_checkm2:
        if data_rows(sample_dir / "checkm2" / "quality_report.tsv") == 0:
            failures.append("CheckM2 required but quality_report.tsv has no data rows")
        if not (sample_dir / "sequence_qc" / "qc_decisions.tsv").exists():
            failures.append("CheckM2 required but sequence_qc/qc_decisions.tsv is missing")
        enriched = sample_dir / "metadata_output" / "ncbi_enriched.csv"
        if not enriched.exists():
            failures.append("CheckM2 required but metadata_output/ncbi_enriched.csv is missing")
        else:
            lines = enriched.read_text(encoding="utf-8", errors="replace").splitlines()
            header = lines[0] if lines else ""
            if "checkm2_completeness" not in header:
                failures.append("CheckM2 required but ncbi_enriched.csv lacks CheckM2 columns")

    if require_genomad:
        status_rows = read_table(sample_dir / "prophage" / "module_status.tsv")
        if not status_rows:
            failures.append("geNomAD required but prophage/module_status.tsv is missing or empty")
        for row in status_rows:
            if is_failure_status(row.get("status", "")):
                failures.append(f"geNomAD module status failure: {row.get('status', '')}")
        if not (features / "prophage.features.tsv").exists():
            failures.append("geNomAD required but prophage.features.tsv is missing")

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path, help="Run directory or organism sample directory")
    parser.add_argument("--require-checkm2", action="store_true", help="Require populated CheckM2 outputs")
    parser.add_argument("--require-genomad", action="store_true", help="Require geNomAD status and feature outputs")
    parser.add_argument(
        "--allow-mobileelementfinder-warning",
        action="store_true",
        help="Allow MobileElementFinder WARNING_FAILED audit rows as nonfatal",
    )
    parser.add_argument(
        "--expect-zero-schema-errors",
        action="store_true",
        help="Require zero unmatched, invalid, and duplicate feature rows",
    )
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.exists():
        raise SystemExit(f"Run directory does not exist: {run_dir}")

    sample_dirs = find_sample_dirs(run_dir)
    if not sample_dirs:
        raise SystemExit(f"No PanR2 handoff outputs found under: {run_dir}")

    all_failures: list[str] = []
    for sample_dir in sample_dirs:
        failures = check_sample_dir(
            sample_dir,
            require_checkm2=args.require_checkm2,
            require_genomad=args.require_genomad,
            allow_mobileelementfinder_warning=args.allow_mobileelementfinder_warning,
            expect_zero_schema_errors=args.expect_zero_schema_errors,
        )
        if failures:
            all_failures.extend(f"{sample_dir}: {failure}" for failure in failures)
        else:
            print(f"PASS\t{sample_dir}")

    if all_failures:
        for failure in all_failures:
            print(f"FAIL\t{failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
