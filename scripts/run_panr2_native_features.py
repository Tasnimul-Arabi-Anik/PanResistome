#!/usr/bin/env python3
"""Run PanR2-compatible annotation helpers under PanResistome ownership."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path


MODULE_STATUS_FIELDS = [
    "module",
    "enabled",
    "started",
    "completed",
    "status",
    "samples_input",
    "samples_processed",
    "samples_failed",
    "raw_tables_created",
    "feature_rows_created",
    "unique_features_created",
    "output_dir",
    "message",
]


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def count_feature_rows(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    with path.open(newline="", encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    feature_ids = {
        str(row.get("GENE", "") or row.get("feature_id", "")).strip()
        for row in rows
        if str(row.get("GENE", "") or row.get("feature_id", "")).strip()
    }
    return len(rows), len(feature_ids)


def count_fastas(sequence_dir: Path) -> int:
    suffixes = {".fa", ".fna", ".fasta", ".fas"}
    count = 0
    for path in sequence_dir.iterdir() if sequence_dir.exists() else []:
        name = path.name.lower()
        if path.is_file() and (path.suffix.lower() in suffixes or any(name.endswith(f"{suffix}.gz") for suffix in suffixes)):
            count += 1
    return count


def status_row(
    module: str,
    enabled: bool,
    started: str,
    status: str,
    output_dir: Path,
    message: str,
    samples_input: int = 0,
    samples_processed: int = 0,
    samples_failed: int = 0,
    raw_tables_created: int = 0,
    feature_rows_created: int = 0,
    unique_features_created: int = 0,
) -> dict[str, str]:
    return {
        "module": module,
        "enabled": str(enabled).lower(),
        "started": started,
        "completed": utc_now(),
        "status": status,
        "samples_input": str(samples_input),
        "samples_processed": str(samples_processed),
        "samples_failed": str(samples_failed),
        "raw_tables_created": str(raw_tables_created),
        "feature_rows_created": str(feature_rows_created),
        "unique_features_created": str(unique_features_created),
        "output_dir": str(output_dir),
        "message": message,
    }


def write_status(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODULE_STATUS_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sequence-dir", required=True, type=Path)
    parser.add_argument("--abricate-dbs", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--force", type=as_bool, default=False)
    parser.add_argument("--run-integronfinder", type=as_bool, default=False)
    parser.add_argument("--run-mlst", type=as_bool, default=False)
    parser.add_argument("--run-mobileelementfinder", type=as_bool, default=False)
    args = parser.parse_args()

    try:
        from panr2.mlst import run_mlst
        from panr2.runners import run_abricate_databases, run_integronfinder, run_mobileelementfinder
    except Exception as exc:
        raise SystemExit(
            "PanR2 is required in the PanR2 comprehensive environment for native feature runners. "
            f"Import failed: {exc}"
        )

    sample_dir = args.sample_dir.resolve()
    sequence_dir = args.sequence_dir.resolve()
    status_dir = sample_dir / "panr2_native_feature_runners"
    status_path = status_dir / "module_status.tsv"
    rows: list[dict[str, str]] = []
    sample_count = count_fastas(sequence_dir)

    if sample_count == 0:
        rows.append(
            status_row(
                "panr2_native_feature_runners",
                True,
                utc_now(),
                "FAIL",
                status_dir,
                f"No FASTA files found in {sequence_dir}",
            )
        )
        write_status(status_path, rows)
        return 1

    started = utc_now()
    databases = [db.strip() for db in args.abricate_dbs.split(",") if db.strip()]
    try:
        result = run_abricate_databases(
            str(sequence_dir),
            str(sample_dir),
            databases,
            summary_metric="identity",
            force=args.force,
        )
        feature_rows = 0
        unique_features = 0
        for db_dir in result.get("database_dirs", {}).values():
            rows_count, unique_count = count_feature_rows(Path(db_dir) / f"{Path(db_dir).name}_results.tab")
            feature_rows += rows_count
            unique_features += unique_count
        rows.append(
            status_row(
                "abricate",
                True,
                started,
                "PASS",
                sample_dir / "tool_results" / "abricate",
                f"ABRicate completed for databases: {','.join(databases)}",
                sample_count,
                sample_count,
                0,
                len(databases),
                feature_rows,
                unique_features,
            )
        )
    except Exception as exc:
        rows.append(status_row("abricate", True, started, "FAIL", sample_dir / "tool_results" / "abricate", str(exc), sample_count, 0, sample_count))
        write_status(status_path, rows)
        raise

    if args.run_integronfinder:
        started = utc_now()
        try:
            result = run_integronfinder(
                str(sequence_dir),
                str(sample_dir),
                cpu=max(args.threads, 1),
                force=args.force,
            )
            feature_rows, unique_features = count_feature_rows(Path(result["feature_dir"]) / "integronfinder_results.tab")
            rows.append(
                status_row(
                    "integronfinder",
                    True,
                    started,
                    "PASS",
                    Path(result["feature_dir"]),
                    "IntegronFinder completed and was converted to PanR2-compatible tables.",
                    sample_count,
                    sample_count,
                    0,
                    len(result.get("raw_tables", [])),
                    feature_rows,
                    unique_features,
                )
            )
        except Exception as exc:
            rows.append(status_row("integronfinder", True, started, "FAIL", sample_dir / "tool_results" / "integronfinder", str(exc), sample_count, 0, sample_count))
            write_status(status_path, rows)
            raise
    else:
        rows.append(status_row("integronfinder", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "integronfinder", "Not enabled for this profile."))

    if args.run_mlst:
        started = utc_now()
        try:
            result = run_mlst(str(sequence_dir), str(sample_dir), force=args.force)
            raw_path = Path(result["raw_table"])
            raw_rows = max(len(raw_path.read_text(encoding="utf-8", errors="ignore").splitlines()), 0) if raw_path.exists() else 0
            rows.append(
                status_row(
                    "mlst",
                    True,
                    started,
                    "PASS",
                    Path(result["mlst_dir"]),
                    "MLST completed and is available for PanR2 analysis.",
                    sample_count,
                    sample_count,
                    0,
                    1 if raw_path.exists() else 0,
                    raw_rows,
                    raw_rows,
                )
            )
        except Exception as exc:
            rows.append(status_row("mlst", True, started, "FAIL", sample_dir / "tool_results" / "mlst", str(exc), sample_count, 0, sample_count))
            write_status(status_path, rows)
            raise
    else:
        rows.append(status_row("mlst", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "mlst", "Not enabled for this profile."))

    if args.run_mobileelementfinder:
        started = utc_now()
        try:
            result = run_mobileelementfinder(
                str(sequence_dir),
                str(sample_dir),
                threads=max(args.threads, 1),
                force=args.force,
            )
            feature_rows, unique_features = count_feature_rows(Path(result["feature_dir"]) / "mobileelementfinder_results.tab")
            rows.append(
                status_row(
                    "mobileelementfinder",
                    True,
                    started,
                    "PASS",
                    Path(result["feature_dir"]),
                    "MobileElementFinder completed and was converted to PanR2-compatible tables.",
                    sample_count,
                    sample_count,
                    0,
                    len(result.get("raw_csv", [])),
                    feature_rows,
                    unique_features,
                )
            )
        except Exception as exc:
            rows.append(status_row("mobileelementfinder", True, started, "FAIL", sample_dir / "tool_results" / "mobileelementfinder", str(exc), sample_count, 0, sample_count))
            write_status(status_path, rows)
            raise
    else:
        rows.append(status_row("mobileelementfinder", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "mobileelementfinder", "Not enabled for this profile."))

    write_status(status_path, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
