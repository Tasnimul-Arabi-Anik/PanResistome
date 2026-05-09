#!/usr/bin/env python3
"""Write a database/tool setup status report for PanResistome runs."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


FIELDS = [
    "database_or_tool",
    "required_for_profile",
    "checked",
    "status",
    "setup_action",
    "version_or_path",
    "message",
]


TRUE_VALUES = {"1", "true", "yes", "y", "on"}
FALSE_VALUES = {"0", "false", "no", "n", "off", ""}


def as_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean value, got {value!r}")


def run_command(command: list[str], timeout: int = 60) -> tuple[int, str]:
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


def command_version(command: str, args: list[str] | None = None) -> str:
    if shutil.which(command) is None:
        return ""
    code, output = run_command([command] + (args or ["--version"]))
    if code != 0 and not output:
        return ""
    return first_line(output)


def row(name: str, required: bool, checked: bool, status: str, setup_action: str = "", version_or_path: str = "", message: str = "") -> dict[str, str]:
    return {
        "database_or_tool": name,
        "required_for_profile": str(bool(required)).lower(),
        "checked": str(bool(checked)).lower(),
        "status": status,
        "setup_action": setup_action,
        "version_or_path": version_or_path,
        "message": message,
    }


def existing_path(value: str | None) -> str:
    if not value:
        return ""
    path = Path(value)
    return str(path) if path.exists() else ""


def first_dmnd(path: str | None) -> str:
    if not path:
        return ""
    root = Path(path)
    if root.is_file() and root.suffix == ".dmnd":
        return str(root)
    if not root.exists():
        return ""
    try:
        found = next(root.rglob("*.dmnd"))
    except StopIteration:
        return ""
    return str(found)


def file_has_data(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    with path.open(errors="ignore") as handle:
        nonempty = [line for line in handle if line.strip()]
    return len(nonempty) > 1


def read_abricate_databases() -> tuple[set[str], str]:
    if shutil.which("abricate") is None:
        return set(), "ABRicate executable not found in this environment."
    code, output = run_command(["abricate", "--list"])
    if code != 0:
        return set(), output or "abricate --list failed."
    databases: set[str] = set()
    for line in output.splitlines():
        parts = line.strip().split()
        if not parts:
            continue
        name = parts[0].strip()
        if name.lower() in {"database", "db", "name"}:
            continue
        databases.add(name)
    return databases, first_line(output)


def count_status_rows(path: Path, status_value: str = "PASS") -> int:
    if not path.exists():
        return 0
    try:
        with path.open(newline="", errors="ignore") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return sum(1 for record in reader if str(record.get("status", "")).strip().upper() == status_value)
    except csv.Error:
        return 0


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    sample_dir = Path(args.sample_dir)
    rows: list[dict[str, str]] = []

    metadata_path = sample_dir / "metadata_output" / "ncbi_clean.csv"
    rows.append(
        row(
            "fetchm2_metadata",
            True,
            True,
            "PASS" if metadata_path.exists() else "FAIL",
            "metadata_generation",
            str(metadata_path) if metadata_path.exists() else "",
            "FetchM2-compatible metadata table found." if metadata_path.exists() else "Missing metadata_output/ncbi_clean.csv.",
        )
    )

    sequence_dir = sample_dir / ("sequence_filtered" if args.qc_filter and (sample_dir / "sequence_filtered").is_dir() else "sequence")
    fasta_count = len(list(sequence_dir.glob("*.fna"))) if sequence_dir.exists() else 0
    rows.append(
        row(
            "sequence_fasta_inputs",
            True,
            True,
            "PASS" if fasta_count > 0 else "FAIL",
            "fetchm2_sequence_download_or_local_samples",
            str(sequence_dir),
            f"{fasta_count} .fna files available for downstream analysis.",
        )
    )

    checkm2_report = sample_dir / "checkm2" / "quality_report.tsv"
    checkm2_db = existing_path(args.checkm2_db) or first_dmnd(args.checkm2_db_dir)
    if args.run_checkm2:
        if file_has_data(checkm2_report):
            rows.append(
                row(
                    "checkm2_database_and_qc",
                    True,
                    True,
                    "PASS",
                    "supplied_or_auto_downloaded_db",
                    checkm2_db or str(checkm2_report),
                    "CheckM2 quality_report.tsv was produced.",
                )
            )
        else:
            rows.append(
                row(
                    "checkm2_database_and_qc",
                    True,
                    True,
                    "FAIL",
                    "auto_download_enabled" if args.checkm2_auto_download_db else "manual_db_required",
                    checkm2_db,
                    "CheckM2 was enabled but checkm2/quality_report.tsv was missing or empty.",
                )
            )
    else:
        rows.append(row("checkm2_database_and_qc", False, False, "SKIPPED", "not_requested", "", "CheckM2 was disabled."))

    for module, enabled, paths in [
        ("quast", args.run_quast, [sample_dir / "quast" / "report.tsv", sample_dir / "quast" / "analysis" / "panr2_quast_summary.csv"]),
        ("ani", args.run_ani, [sample_dir / "ani" / "analysis" / "panr2_ani_summary.csv"]),
        ("mash", args.run_mash, [sample_dir / "mash" / "analysis" / "closest_mash_neighbor.csv", sample_dir / "mash" / "analysis" / "mash_distance_long.csv"]),
    ]:
        if not enabled:
            rows.append(row(module, False, False, "SKIPPED", "not_requested", "", f"{module} was disabled."))
            continue
        found = next((path for path in paths if path.exists()), None)
        status = "PASS" if found else "FAIL"
        message = f"{module} output found." if found else f"{module} was enabled but expected output was missing."
        if module == "ani" and found:
            ani_status = sample_dir / "ani" / "analysis" / "ani_run_status.tsv"
            if ani_status.exists():
                with ani_status.open(newline="", encoding="utf-8") as handle:
                    status_rows = list(csv.DictReader(handle, delimiter="\t"))
                if status_rows and status_rows[0].get("status", "").startswith("SKIPPED"):
                    status = "WARNING_SKIPPED"
                    message = status_rows[0].get("message") or "ANI summary exists, but pairwise ANI was skipped."
        rows.append(
            row(
                module,
                True,
                True,
                status,
                "module_output_check",
                str(found) if found else "",
                message,
            )
        )

    if args.run_gtdbtk:
        db_path = existing_path(args.gtdbtk_data_path)
        report = sample_dir / "gtdbtk" / "gtdbtk.bac120.summary.tsv"
        rows.append(
            row(
                "gtdbtk_database",
                True,
                True,
                "PASS" if db_path and report.exists() else "FAIL",
                "external_db_required",
                db_path,
                "GTDB-Tk database and summary output found." if db_path and report.exists() else "GTDB-Tk was enabled; provide --gtdbtk_data_path and confirm output.",
            )
        )
    else:
        rows.append(row("gtdbtk_database", False, False, "SKIPPED", "heavy_db_not_default", "", "GTDB-Tk is disabled by default."))

    panr_version = command_version("panr", ["--version"])
    rows.append(
        row(
            "panr2",
            args.run_panr2_comprehensive,
            True,
            "PASS" if panr_version else ("FAIL" if args.run_panr2_comprehensive else "SKIPPED"),
            "conda_environment",
            panr_version,
            "PanR2 CLI is available." if panr_version else "PanR2 CLI was not found.",
        )
    )

    abricate_version = command_version("abricate", ["--version"])
    abricate_databases, abricate_message = read_abricate_databases()
    required_dbs = [item.strip() for item in str(args.panr2_dbs or "").split(",") if item.strip()]
    rows.append(
        row(
            "abricate_tool",
            bool(required_dbs),
            True,
            "PASS" if abricate_version and abricate_databases else ("FAIL" if required_dbs else "SKIPPED"),
            "panr_setup_db_then_list",
            abricate_version,
            abricate_message,
        )
    )
    for db_name in required_dbs:
        rows.append(
            row(
                f"abricate_db:{db_name}",
                True,
                True,
                "PASS" if db_name in abricate_databases else "FAIL",
                "panr_setup_db",
                db_name if db_name in abricate_databases else "",
                f"ABRicate database {db_name!r} is available." if db_name in abricate_databases else f"ABRicate database {db_name!r} is required but unavailable after setup.",
            )
        )

    for command, module, enabled in [
        ("integron_finder", "integronfinder", args.run_integronfinder),
        ("mlst", "mlst", args.run_mlst),
        ("mefinder", "mobileelementfinder", args.run_mobileelementfinder),
        ("defense-finder", "defensefinder", args.run_defensefinder),
    ]:
        version = command_version(command, ["--version"])
        if enabled:
            rows.append(
                row(
                    module,
                    True,
                    True,
                    "PASS" if version or shutil.which(command) else "FAIL",
                    "conda_environment",
                    version or (shutil.which(command) or ""),
                    f"{command} is available." if version or shutil.which(command) else f"{command} was requested but not found.",
                )
            )
        else:
            message = "MobileElementFinder is opt-in because upstream parser failures were observed." if module == "mobileelementfinder" else f"{module} was not requested."
            rows.append(row(module, False, False, "SKIPPED", "not_requested", "", message))

    if args.run_isfinder:
        fasta = existing_path(args.isfinder_db_fasta)
        status_path = sample_dir / "isfinder" / "module_status.tsv"
        pass_rows = count_status_rows(status_path)
        rows.append(
            row(
                "isfinder_authorized_fasta",
                True,
                True,
                "PASS" if fasta and pass_rows else "FAIL",
                "user_supplied_authorized_fasta",
                fasta,
                "Authorized ISfinder FASTA and module PASS status found." if fasta and pass_rows else "ISfinder requires --isfinder_db_fasta and a successful isfinder/module_status.tsv.",
            )
        )
    else:
        rows.append(row("isfinder_authorized_fasta", False, False, "SKIPPED", "legal_restriction", "", "ISfinder is not auto-downloaded or redistributed; supply --isfinder_db_fasta to enable it."))

    if args.run_amrfinderplus:
        status_path = sample_dir / "amrfinderplus" / "tables" / "amrfinderplus_sample_status.tsv"
        pass_rows = count_status_rows(status_path)
        rows.append(
            row(
                "amrfinderplus_database_and_runner",
                True,
                True,
                "PASS" if pass_rows else "FAIL",
                "amrfinder_update_db" if args.amrfinderplus_update_db else "preinstalled_db",
                str(status_path) if status_path.exists() else "",
                f"AMRFinderPlus completed for {pass_rows} samples." if pass_rows else "AMRFinderPlus was enabled but no PASS sample status rows were found.",
            )
        )
    else:
        rows.append(row("amrfinderplus_database_and_runner", False, False, "SKIPPED", "not_requested", "", "AMRFinderPlus was not requested."))

    for module, enabled, required_path, message in [
        ("mobsuite", args.run_mobsuite, sample_dir / "mobsuite" / "tables", "MOB-suite was requested but no tables were found."),
        ("genomad_database", args.run_genomad, Path(args.genomad_db) if args.genomad_db else Path(""), "geNomad requires --genomad_db when enabled."),
        ("kaptive_database", args.run_kaptive, Path(args.kaptive_db) if args.kaptive_db else Path(""), "Kaptive requires --kaptive_db when enabled."),
    ]:
        if not enabled:
            rows.append(row(module, False, False, "SKIPPED", "not_requested", "", f"{module} was disabled."))
            continue
        rows.append(
            row(
                module,
                True,
                True,
                "PASS" if str(required_path) and required_path.exists() else "FAIL",
                "external_db_or_module_output",
                str(required_path) if str(required_path) and required_path.exists() else "",
                f"{module} requirement found." if str(required_path) and required_path.exists() else message,
            )
        )

    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--analysis-profile", default="custom")
    parser.add_argument("--panr2-dbs", default="")
    parser.add_argument("--qc-filter", type=as_bool, default=False)
    parser.add_argument("--run-checkm2", type=as_bool, default=False)
    parser.add_argument("--checkm2-db", default="")
    parser.add_argument("--checkm2-db-dir", default="")
    parser.add_argument("--checkm2-auto-download-db", type=as_bool, default=True)
    parser.add_argument("--run-gtdbtk", type=as_bool, default=False)
    parser.add_argument("--gtdbtk-data-path", default="")
    parser.add_argument("--run-quast", type=as_bool, default=False)
    parser.add_argument("--run-ani", type=as_bool, default=False)
    parser.add_argument("--run-mash", type=as_bool, default=False)
    parser.add_argument("--run-panr2-comprehensive", type=as_bool, default=False)
    parser.add_argument("--run-integronfinder", type=as_bool, default=False)
    parser.add_argument("--run-mlst", type=as_bool, default=False)
    parser.add_argument("--run-mobileelementfinder", type=as_bool, default=False)
    parser.add_argument("--run-defensefinder", type=as_bool, default=False)
    parser.add_argument("--run-isfinder", type=as_bool, default=False)
    parser.add_argument("--isfinder-db-fasta", default="")
    parser.add_argument("--run-amrfinderplus", type=as_bool, default=False)
    parser.add_argument("--amrfinderplus-update-db", type=as_bool, default=True)
    parser.add_argument("--run-mobsuite", type=as_bool, default=False)
    parser.add_argument("--run-genomad", type=as_bool, default=False)
    parser.add_argument("--genomad-db", default="")
    parser.add_argument("--run-kaptive", type=as_bool, default=False)
    parser.add_argument("--kaptive-db", default="")
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = build_rows(args)
    out = Path(args.out)
    write_rows(out, rows)
    failures = [record for record in rows if record["required_for_profile"] == "true" and record["status"] == "FAIL"]
    if args.strict and failures:
        print(f"Database/tool preflight failed; see {out}", file=sys.stderr)
        for failure in failures:
            print(f"- {failure['database_or_tool']}: {failure['message']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
