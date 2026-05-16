#!/usr/bin/env python3
"""Run PanR2-compatible annotation helpers under PanResistome ownership."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
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

NATIVE_RUNNER_AUDIT_FIELDS = [
    "module",
    "runner_mode",
    "expected_raw_tables",
    "observed_raw_tables",
    "samples_input",
    "samples_processed",
    "samples_failed",
    "feature_rows",
    "unique_features",
    "status",
    "message",
]


def as_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def worker_count(threads: int, tasks: int) -> int:
    if tasks <= 0:
        return 1
    return max(1, min(max(int(threads or 1), 1), tasks))


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


def write_empty_abricate_style_table(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "#FILE",
        "SEQUENCE",
        "START",
        "END",
        "GENE",
        "COVERAGE",
        "COVERAGE_MAP",
        "GAPS",
        "%COVERAGE",
        "%IDENTITY",
        "DATABASE",
        "ACCESSION",
        "PRODUCT",
        "RESISTANCE",
    ]
    path.write_text("\t".join(header) + "\n", encoding="utf-8")


def is_missing_value(value: str) -> bool:
    text = str(value or "").strip()
    return not text or text in {"-", ".", "?"} or text.lower() in {"na", "n/a", "nan", "none", "null", "unknown"}


def is_placeholder_mlst_feature(value: str) -> bool:
    text = str(value or "").strip()
    if is_missing_value(text):
        return True
    return re.fullmatch(r"[-_:\s]*ST[-_:\s]*", text, flags=re.IGNORECASE) is not None


def count_mlst_feature_rows(path: Path) -> tuple[int, int]:
    """Count biological MLST features in native headerless `mlst` output.

    The `mlst` command writes one row per input FASTA even when the organism has
    no matching PubMLST scheme.  Those unsupported rows are useful run evidence
    but should not be counted as sequence-type features.
    """
    if not path.exists():
        return 0, 0
    feature_ids: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [part.strip() for part in line.split("\t")]
        if len(parts) < 3:
            continue
        st = parts[2]
        if not is_missing_value(st):
            st_feature = f"ST_{st}"
            if not is_placeholder_mlst_feature(st_feature):
                feature_ids.append(st_feature)
        for allele in parts[3:]:
            match = re.fullmatch(r"([^()]+)\(([^()]+)\)", allele)
            if not match:
                continue
            locus, allele_number = match.groups()
            if is_missing_value(locus) or is_missing_value(allele_number):
                continue
            feature_ids.append(f"{locus}_{allele_number}")
    return len(feature_ids), len(set(feature_ids))


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


def audit_row(
    module: str,
    runner_mode: str,
    expected_raw_tables: int | str,
    observed_raw_tables: int | str,
    samples_input: int | str,
    samples_processed: int | str,
    samples_failed: int | str,
    feature_rows: int | str,
    unique_features: int | str,
    status: str,
    message: str,
) -> dict[str, str]:
    return {
        "module": module,
        "runner_mode": runner_mode,
        "expected_raw_tables": str(expected_raw_tables),
        "observed_raw_tables": str(observed_raw_tables),
        "samples_input": str(samples_input),
        "samples_processed": str(samples_processed),
        "samples_failed": str(samples_failed),
        "feature_rows": str(feature_rows),
        "unique_features": str(unique_features),
        "status": status,
        "message": message,
    }


def write_audit(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NATIVE_RUNNER_AUDIT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_mobileelementfinder_failure_outputs(
    sample_dir: Path,
    status_path: Path,
    audit_path: Path,
    rows: list[dict[str, str]],
    audit_rows: list[dict[str, str]],
    runner_mode: str,
    started: str,
    sample_count: int,
    allow_failure: bool,
    error: Exception | str,
) -> str:
    output_dir = sample_dir / "tool_results" / "mobileelementfinder" / "panr2_inputs"
    write_empty_abricate_style_table(output_dir / "mobileelementfinder_results.tab")
    status = "WARNING_FAILED" if allow_failure else "FAIL"
    message = (
        "MobileElementFinder failed but was kept nonfatal; header-only PanR2-compatible output was written. "
        "Inspect native runner status and rerun without --panr2_run_mobileelementfinder if not needed. "
        f"Original error: {error}"
    )
    print(f"Warning: {message}", file=sys.stderr)
    rows.append(status_row("mobileelementfinder", True, started, status, output_dir, message, sample_count, 0, sample_count))
    write_status(status_path, rows)
    audit_rows.append(audit_row("mobileelementfinder", runner_mode, sample_count, 1, sample_count, 0, sample_count, 0, 0, status, message))
    write_audit(audit_path, audit_rows)
    return status


def _run_command(command: list[str], stdout_path: Path | None = None) -> str:
    if stdout_path:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
    else:
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or getattr(completed, "stdout", "").strip()
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}\n{detail}")
    return completed.stdout if not stdout_path else ""


def _capture_command(command: list[str]) -> str:
    try:
        completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        output = (completed.stdout or completed.stderr or "").strip()
        return output if completed.returncode == 0 else f"unavailable: {output}"
    except Exception as exc:
        return f"unavailable: {exc}"


def run_abricate_parallel(sequence_dir: Path, sample_dir: Path, databases: list[str], threads: int, force: bool = False) -> dict:
    from panr2.runners import _parse_abricate_list, find_sequence_files, write_tool_manifest

    executable = shutil.which("abricate")
    if not executable:
        raise FileNotFoundError("ABRicate executable not found: abricate")
    sequence_files = find_sequence_files(str(sequence_dir))
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")
    list_output = _capture_command([executable, "--list"])
    available = _parse_abricate_list(list_output)
    if not available:
        raise ValueError("ABRicate did not report any available databases. Run `panr setup-db` or `abricate --setupdb`.")
    missing = [db for db in databases if db not in available]
    if missing:
        raise ValueError(f"ABRicate database(s) not available: {', '.join(missing)}")

    base_dir = sample_dir / "tool_results" / "abricate"
    base_dir.mkdir(parents=True, exist_ok=True)
    version = _capture_command([executable, "--version"]).splitlines()[0]

    def sequence_label(sequence_file: str) -> str:
        label = Path(sequence_file).name
        for suffix in [".gz", ".fasta", ".fna", ".fa", ".fas"]:
            if label.lower().endswith(suffix):
                label = label[: -len(suffix)]
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", label) or "sample"

    def combine_result_tables(input_paths: list[Path], output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        wrote_header = False
        with output_path.open("w", encoding="utf-8") as output_handle:
            for input_path in input_paths:
                if not input_path.exists():
                    continue
                with input_path.open(encoding="utf-8", errors="ignore") as input_handle:
                    for line in input_handle:
                        if not line.strip():
                            continue
                        is_header = line.startswith("#FILE\t") or line.startswith("FILE\t")
                        if is_header:
                            if not wrote_header:
                                output_handle.write(line)
                                wrote_header = True
                            continue
                        output_handle.write(line)

    def run_sample_database(db: str, sequence_file: str) -> tuple[str, str, Path]:
        db_dir = base_dir / db
        raw_dir = db_dir / "per_sample"
        raw_dir.mkdir(parents=True, exist_ok=True)
        raw_path = raw_dir / f"{sequence_label(sequence_file)}.tab"
        if force or not raw_path.exists():
            _run_command([executable, "--db", db, sequence_file], stdout_path=raw_path)
        return db, sequence_file, raw_path

    workers = worker_count(threads, len(sequence_files))
    per_db_paths: dict[str, list[Path]] = {db: [] for db in databases}
    for db in databases:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(run_sample_database, db, sequence_file): sequence_file for sequence_file in sequence_files}
            for future in as_completed(futures):
                _db, _sequence_file, raw_path = future.result()
                per_db_paths[db].append(raw_path)

    def finalize_database(db: str) -> tuple[str, Path]:
        db_dir = base_dir / db
        db_dir.mkdir(parents=True, exist_ok=True)
        results_path = db_dir / f"{db}_results.tab"
        summary_path = db_dir / f"{db}_summary.tab"
        ordered_paths = [
            db_dir / "per_sample" / f"{sequence_label(sequence_file)}.tab"
            for sequence_file in sequence_files
        ]
        if force or not results_path.exists():
            combine_result_tables(ordered_paths, results_path)
        if force or not summary_path.exists():
            _run_command([executable, "--summary", "--identity", str(results_path)], stdout_path=summary_path)
        return db, db_dir

    output_dirs: dict[str, str] = {}
    for db in databases:
        db, db_dir = finalize_database(db)
        output_dirs[db] = str(db_dir)

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sequence_dir": str(sequence_dir),
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "abricate",
            "executable": executable,
            "version": version,
            "runs": [],
        }],
    }
    for db in databases:
        db_dir = Path(output_dirs[db])
        db_info = available.get(db, {})
        manifest["tools"][0]["runs"].append({
            "database": db,
            "database_sequences": db_info.get("SEQUENCES", ""),
            "database_date": db_info.get("DATE", ""),
            "database_type": db_info.get("DBTYPE", ""),
            "results": str(db_dir / f"{db}_results.tab"),
            "summary": str(db_dir / f"{db}_summary.tab"),
            "status": "completed",
        })
    manifest_paths = write_tool_manifest(str(sample_dir), manifest)
    return {
        "database_dirs": output_dirs,
        "manifest": manifest_paths,
        "raw_tables": [str(path) for paths in per_db_paths.values() for path in paths],
        "parallel_workers": workers,
    }


def run_integronfinder_parallel(sequence_dir: Path, sample_dir: Path, threads: int, force: bool = False) -> dict:
    from panr2.runners import _find_integronfinder_table, convert_integronfinder_outputs, find_sequence_files, write_tool_manifest

    executable = shutil.which("integron_finder")
    if not executable:
        raise FileNotFoundError("IntegronFinder executable not found: integron_finder")
    sequence_files = find_sequence_files(str(sequence_dir))
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")
    raw_dir = sample_dir / "tool_results" / "integronfinder" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    version = _capture_command([executable, "--version"]).splitlines()[0]

    def run_sequence(sequence_file: str) -> str | None:
        prefix = Path(sequence_file).name
        for suffix in [".gz", ".fasta", ".fna", ".fa", ".fas"]:
            if prefix.lower().endswith(suffix):
                prefix = prefix[: -len(suffix)]
        sample_out_dir = raw_dir / prefix
        sample_out_dir.mkdir(parents=True, exist_ok=True)
        expected_table = _find_integronfinder_table(str(sample_out_dir), prefix)
        if force or not expected_table:
            _run_command([executable, sequence_file, "--outdir", str(sample_out_dir), "--cpu", "1"])
            expected_table = _find_integronfinder_table(str(sample_out_dir), prefix)
        return expected_table

    raw_table_paths: list[str] = []
    with ThreadPoolExecutor(max_workers=worker_count(threads, len(sequence_files))) as executor:
        futures = {executor.submit(run_sequence, sequence_file): sequence_file for sequence_file in sequence_files}
        for future in as_completed(futures):
            table = future.result()
            if table:
                raw_table_paths.append(table)
    raw_table_paths = sorted(set(raw_table_paths))
    converted = convert_integronfinder_outputs(sequence_files, raw_table_paths, str(sample_dir))
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sequence_dir": str(sequence_dir),
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "integronfinder",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "integronfinder",
                "database_sequences": "",
                "database_date": "",
                "results": converted["results"],
                "summary": converted["summary"],
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(str(sample_dir), manifest)
    return {"feature_dir": converted["feature_dir"], "manifest": manifest_paths, "raw_tables": raw_table_paths}


def run_mlst_parallel(sequence_dir: Path, sample_dir: Path, threads: int, force: bool = False) -> dict:
    from panr2.runners import find_sequence_files, write_tool_manifest

    executable = shutil.which("mlst")
    if not executable:
        raise FileNotFoundError("MLST executable not found: mlst")
    sequence_files = find_sequence_files(str(sequence_dir))
    if not sequence_files:
        raise FileNotFoundError(f"No FASTA files found in {sequence_dir}")
    raw_dir = sample_dir / "tool_results" / "mlst" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_path = raw_dir / "mlst.tsv"

    def run_sequence(sequence_file: str) -> str:
        return _run_command([executable, sequence_file])

    if force or not raw_path.exists():
        outputs: dict[str, str] = {}
        with ThreadPoolExecutor(max_workers=worker_count(threads, len(sequence_files))) as executor:
            futures = {executor.submit(run_sequence, sequence_file): sequence_file for sequence_file in sequence_files}
            for future in as_completed(futures):
                outputs[futures[future]] = future.result()
        raw_path.write_text("".join(outputs.get(sequence_file, "") for sequence_file in sequence_files), encoding="utf-8")

    version = "available"
    version_output = _capture_command([executable, "--version"])
    if version_output:
        version = version_output.splitlines()[0]
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sequence_dir": str(sequence_dir),
        "sequence_count": len(sequence_files),
        "tools": [{
            "name": "mlst",
            "executable": executable,
            "version": version,
            "runs": [{
                "database": "pubmlst",
                "database_sequences": "",
                "database_date": "",
                "results": str(raw_path),
                "summary": str(raw_path),
                "status": "completed",
            }],
        }],
    }
    manifest_paths = write_tool_manifest(str(sample_dir), manifest)
    return {"mlst_dir": str(raw_dir), "raw_table": str(raw_path), "manifest": manifest_paths}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True, type=Path)
    parser.add_argument("--sequence-dir", required=True, type=Path)
    parser.add_argument("--abricate-dbs", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--mode", choices=["serial", "parallel"], default="serial")
    parser.add_argument("--force", type=as_bool, default=False)
    parser.add_argument("--run-integronfinder", type=as_bool, default=False)
    parser.add_argument("--run-mlst", type=as_bool, default=False)
    parser.add_argument("--run-mobileelementfinder", type=as_bool, default=False)
    parser.add_argument("--mobileelementfinder-allow-failure", type=as_bool, default=True)
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
    audit_path = status_dir / "native_runner_merge_audit.tsv"
    rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
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
        write_audit(audit_path, [
            audit_row(
                "panr2_native_feature_runners",
                args.mode,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                "FAIL",
                f"No FASTA files found in {sequence_dir}",
            )
        ])
        return 1

    started = utc_now()
    databases = [db.strip() for db in args.abricate_dbs.split(",") if db.strip()]
    try:
        if args.mode == "parallel":
            result = run_abricate_parallel(sequence_dir, sample_dir, databases, args.threads, force=args.force)
            abricate_message = (
                f"ABRicate completed for databases: {','.join(databases)} "
                f"(per-database sample-parallel workers={result.get('parallel_workers', worker_count(args.threads, sample_count))})"
            )
        else:
            result = run_abricate_databases(
                str(sequence_dir),
                str(sample_dir),
                databases,
                summary_metric="identity",
                force=args.force,
            )
            abricate_message = f"ABRicate completed for databases: {','.join(databases)}"
        feature_rows = 0
        unique_features = 0
        for db_dir in result.get("database_dirs", {}).values():
            rows_count, unique_count = count_feature_rows(Path(db_dir) / f"{Path(db_dir).name}_results.tab")
            feature_rows += rows_count
            unique_features += unique_count
        raw_tables = len(result.get("raw_tables", [])) or len(databases)
        module_row = status_row(
            "abricate",
            True,
            started,
            "PASS",
            sample_dir / "tool_results" / "abricate",
            abricate_message,
            sample_count,
            sample_count,
            0,
            raw_tables,
            feature_rows,
            unique_features,
        )
        rows.append(module_row)
        expected_raw_tables = sample_count * len(databases) if args.mode == "parallel" else len(databases)
        audit_rows.append(
            audit_row(
                "abricate",
                args.mode,
                expected_raw_tables,
                raw_tables,
                sample_count,
                sample_count,
                0,
                feature_rows,
                unique_features,
                "PASS" if raw_tables >= expected_raw_tables else "WARNING",
                abricate_message,
            )
        )
    except Exception as exc:
        rows.append(status_row("abricate", True, started, "FAIL", sample_dir / "tool_results" / "abricate", str(exc), sample_count, 0, sample_count))
        write_status(status_path, rows)
        audit_rows.append(audit_row("abricate", args.mode, sample_count * len(databases) if args.mode == "parallel" else len(databases), 0, sample_count, 0, sample_count, 0, 0, "FAIL", str(exc)))
        write_audit(audit_path, audit_rows)
        raise

    if args.run_integronfinder:
        started = utc_now()
        try:
            if args.mode == "parallel":
                result = run_integronfinder_parallel(sequence_dir, sample_dir, args.threads, force=args.force)
                integron_message = (
                    "IntegronFinder completed in per-assembly parallel mode and was converted to "
                    f"PanR2-compatible tables (parallel workers={worker_count(args.threads, sample_count)})."
                )
            else:
                result = run_integronfinder(
                    str(sequence_dir),
                    str(sample_dir),
                    cpu=max(args.threads, 1),
                    force=args.force,
                )
                integron_message = "IntegronFinder completed and was converted to PanR2-compatible tables."
            feature_rows, unique_features = count_feature_rows(Path(result["feature_dir"]) / "integronfinder_results.tab")
            raw_tables = len(result.get("raw_tables", []))
            module_row = status_row(
                "integronfinder",
                True,
                started,
                "PASS",
                Path(result["feature_dir"]),
                integron_message,
                sample_count,
                sample_count,
                0,
                raw_tables,
                feature_rows,
                unique_features,
            )
            rows.append(module_row)
            audit_rows.append(
                audit_row(
                    "integronfinder",
                    args.mode,
                    sample_count,
                    raw_tables,
                    sample_count,
                    sample_count,
                    0,
                    feature_rows,
                    unique_features,
                    "PASS" if raw_tables >= sample_count else "WARNING",
                    integron_message,
                )
            )
        except Exception as exc:
            rows.append(status_row("integronfinder", True, started, "FAIL", sample_dir / "tool_results" / "integronfinder", str(exc), sample_count, 0, sample_count))
            write_status(status_path, rows)
            audit_rows.append(audit_row("integronfinder", args.mode, sample_count, 0, sample_count, 0, sample_count, 0, 0, "FAIL", str(exc)))
            write_audit(audit_path, audit_rows)
            raise
    else:
        rows.append(status_row("integronfinder", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "integronfinder", "Not enabled for this profile."))
        audit_rows.append(audit_row("integronfinder", args.mode, 0, 0, sample_count, 0, 0, 0, 0, "SKIPPED", "Not enabled for this profile."))

    if args.run_mlst:
        started = utc_now()
        try:
            if args.mode == "parallel":
                result = run_mlst_parallel(sequence_dir, sample_dir, args.threads, force=args.force)
                mlst_message = f"MLST completed in per-assembly parallel mode (parallel workers={worker_count(args.threads, sample_count)})."
            else:
                result = run_mlst(str(sequence_dir), str(sample_dir), force=args.force)
                mlst_message = "MLST completed and is available for PanR2 analysis."
            raw_path = Path(result["raw_table"])
            feature_rows, unique_features = count_mlst_feature_rows(raw_path)
            raw_tables = 1 if raw_path.exists() else 0
            module_row = status_row(
                "mlst",
                True,
                started,
                "PASS",
                Path(result["mlst_dir"]),
                mlst_message,
                sample_count,
                sample_count,
                0,
                raw_tables,
                feature_rows,
                unique_features,
            )
            rows.append(module_row)
            audit_rows.append(
                audit_row(
                    "mlst",
                    args.mode,
                    1,
                    raw_tables,
                    sample_count,
                    sample_count,
                    0,
                    feature_rows,
                    unique_features,
                    "PASS" if raw_tables >= 1 else "WARNING",
                    mlst_message,
                )
            )
        except Exception as exc:
            rows.append(status_row("mlst", True, started, "FAIL", sample_dir / "tool_results" / "mlst", str(exc), sample_count, 0, sample_count))
            write_status(status_path, rows)
            audit_rows.append(audit_row("mlst", args.mode, 1, 0, sample_count, 0, sample_count, 0, 0, "FAIL", str(exc)))
            write_audit(audit_path, audit_rows)
            raise
    else:
        rows.append(status_row("mlst", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "mlst", "Not enabled for this profile."))
        audit_rows.append(audit_row("mlst", args.mode, 0, 0, sample_count, 0, 0, 0, 0, "SKIPPED", "Not enabled for this profile."))

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
            raw_tables = len(result.get("raw_csv", []))
            module_row = status_row(
                "mobileelementfinder",
                True,
                started,
                "PASS",
                Path(result["feature_dir"]),
                "MobileElementFinder completed and was converted to PanR2-compatible tables.",
                sample_count,
                sample_count,
                0,
                raw_tables,
                feature_rows,
                unique_features,
            )
            rows.append(module_row)
            audit_rows.append(
                audit_row(
                    "mobileelementfinder",
                    args.mode,
                    sample_count,
                    raw_tables,
                    sample_count,
                    sample_count,
                    0,
                    feature_rows,
                    unique_features,
                    "PASS" if raw_tables >= sample_count else "WARNING",
                    "MobileElementFinder completed and was converted to PanR2-compatible tables.",
                )
            )
        except Exception as exc:
            write_mobileelementfinder_failure_outputs(
                sample_dir,
                status_path,
                audit_path,
                rows,
                audit_rows,
                args.mode,
                started,
                sample_count,
                args.mobileelementfinder_allow_failure,
                exc,
            )
            if not args.mobileelementfinder_allow_failure:
                raise
    else:
        rows.append(status_row("mobileelementfinder", False, utc_now(), "SKIPPED", sample_dir / "tool_results" / "mobileelementfinder", "Not enabled for this profile."))
        audit_rows.append(audit_row("mobileelementfinder", args.mode, 0, 0, sample_count, 0, 0, 0, 0, "SKIPPED", "Not enabled for this profile."))

    write_status(status_path, rows)
    write_audit(audit_path, audit_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
