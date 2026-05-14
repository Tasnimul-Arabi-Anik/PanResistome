#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from panr2_contract import CONTRACT_COLUMNS, FEATURE_COLUMNS, CONTRACT_VERSION, SCHEMA_VERSION, export_contract


def copy_if_exists(src, dst):
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def copytree_if_exists(src, dst):
    if src.exists() and src.is_dir():
        if dst.exists():
            shutil.rmtree(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(src, dst)


def parse_version_line(line, fallback_component):
    if not line or line.startswith("WARNING:"):
        return None
    if line.startswith("[") and line.endswith("]"):
        return {"component": line.strip("[]"), "version": ""}
    if "==" in line:
        component, version = line.split("==", 1)
        return {"component": component.strip(), "version": version.strip()}
    if "=" in line:
        component, version = line.split("=", 1)
        return {"component": component.strip(), "version": version.strip()}
    match = re.match(r"^([A-Za-z][A-Za-z0-9_.+-]*)\s+(?:version\s+)?v?(.+)$", line)
    if match and re.search(r"\d", match.group(2)):
        return {"component": match.group(1).strip(), "version": match.group(2).strip()}
    return {"component": fallback_component, "version": line}


def read_table(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8", errors="ignore") as handle:
        sample = handle.readline()
        if not sample:
            return []
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        handle.seek(0)
        return [
            {str(key): str(value or "") for key, value in row.items() if key is not None}
            for row in csv.DictReader(handle, delimiter=delimiter)
        ]


def git_value(repo_dir, *args):
    if not repo_dir:
        return "unknown"
    try:
        completed = subprocess.run(
            ["git", "-C", str(repo_dir), *args],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    value = completed.stdout.strip()
    return value if completed.returncode == 0 and value else "unknown"


def container_digest(image):
    text = str(image or "")
    if "@sha256:" in text:
        return "sha256:" + text.split("@sha256:", 1)[1]
    return "unknown"


def derive_container_engine(profile_stack):
    profiles = {item.strip().lower() for item in str(profile_stack or "").split(",") if item.strip()}
    for engine in ["docker", "apptainer", "singularity"]:
        if engine in profiles:
            return engine
    return "unknown"


def derive_container_image(profile_stack, explicit_image):
    if explicit_image:
        return explicit_image
    engine = derive_container_engine(profile_stack)
    if engine == "docker":
        return "ghcr.io/tasnimul-arabi-anik/panresistome:experimental"
    if engine in {"apptainer", "singularity"}:
        return "docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental"
    return "unknown"


def tsv_settings(path):
    settings = {}
    for row in read_table(path):
        key = row.get("setting") or row.get("database_or_tool") or row.get("module") or row.get("database")
        if key:
            settings[key] = row.get("value") or row.get("status") or row.get("message") or ""
    return settings


def feature_row_counts(out):
    rows_by_database = {}
    features_dir = out / "features"
    for path in sorted(features_dir.glob("*.features.tsv")):
        rows_by_database[path.name.replace(".features.tsv", "")] = len(read_table(path))
    all_features = features_dir / "all_features.tsv"
    if all_features.exists():
        rows_by_database["all_features"] = len(read_table(all_features))
    return rows_by_database


def write_reproducibility_manifest(sample_dir, out, args):
    manifest_dir = out / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    repo_dir = Path(args.repo_dir).resolve() if args.repo_dir else None
    pipeline_outdir = Path(args.pipeline_outdir) if args.pipeline_outdir else Path("unknown")
    git_commit = args.git_commit or git_value(repo_dir, "rev-parse", "HEAD")
    git_tag = args.git_tag or git_value(repo_dir, "describe", "--tags", "--exact-match")
    git_describe = git_value(repo_dir, "describe", "--tags", "--always", "--dirty")
    engine = args.container_engine or derive_container_engine(args.profile_stack)
    image = derive_container_image(args.profile_stack, args.container_image)
    database_setup_path = manifest_dir / "database_setup_status.tsv"
    report_controls_path = manifest_dir / "report_controls.tsv"
    module_status_path = manifest_dir / "module_status_summary.tsv"
    feature_contract_path = manifest_dir / "feature_contract.json"
    manifest = {
        "manifest_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline": {
            "name": "PanResistome",
            "version": args.pipeline_version or "unknown",
            "git_commit": git_commit,
            "git_tag": git_tag,
            "git_describe": git_describe,
            "contract_version": CONTRACT_VERSION,
            "schema_version": SCHEMA_VERSION,
        },
        "execution": {
            "command_line": args.run_command or "unknown",
            "profile_stack": args.profile_stack or "unknown",
            "launch_dir": args.launch_dir or "unknown",
            "sample_dir": str(sample_dir),
            "pipeline_outdir": args.pipeline_outdir or "unknown",
            "nextflow_session_id": args.nextflow_session_id or os.environ.get("NXF_SESSION_ID", "unknown"),
            "nextflow_run_name": args.nextflow_run_name or "unknown",
        },
        "container": {
            "engine": engine,
            "image": image,
            "digest": args.container_digest or container_digest(image),
        },
        "features": {
            "row_counts_by_table": feature_row_counts(out),
            "all_features": str(out / "features" / "all_features.tsv"),
        },
        "manifests": {
            "software_versions": str(manifest_dir / "software_versions.csv"),
            "database_setup_status": str(database_setup_path),
            "module_status_summary": str(module_status_path),
            "feature_completeness_audit": str(manifest_dir / "feature_completeness_audit.tsv"),
            "schema_validation_summary": str(manifest_dir / "schema_validation_summary.txt"),
            "feature_contract": str(feature_contract_path),
            "report_controls": str(report_controls_path),
            "runtime_summary": str(pipeline_outdir / "pipeline_runtime_summary.tsv") if args.pipeline_outdir else "unknown",
            "runtime_tasks": str(pipeline_outdir / "pipeline_runtime_tasks.tsv") if args.pipeline_outdir else "unknown",
        },
        "report_controls": tsv_settings(report_controls_path),
        "database_setup_status": tsv_settings(database_setup_path),
        "module_status_summary": tsv_settings(module_status_path),
        "notes": [
            "Complete standardized feature TSVs remain authoritative.",
            "Report-facing summaries may be capped when large-dataset safeguards are enabled.",
            "Runtime summary files are generated by the Nextflow completion hook after workflow tasks finish.",
        ],
    }
    path = manifest_dir / "reproducibility_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def main():
    parser = argparse.ArgumentParser(description="Export a PanResistome sample directory into PanR2-ready inputs.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--versions-dir", help="Optional pipeline_versions directory to include in the PanR2 manifest.")
    parser.add_argument("--large-dataset", action="store_true", help="Enable large-dataset report safeguards in the PanR2 handoff bundle.")
    parser.add_argument("--report-mode", choices=["compact", "publication", "exploratory"], default="publication")
    parser.add_argument("--max-features-heatmap", type=int, default=300)
    parser.add_argument("--max-features-network", type=int, default=300)
    parser.add_argument("--max-metadata-columns", type=int, default=80)
    parser.add_argument("--top-n-features-per-database", type=int, default=25)
    parser.add_argument("--skip-heavy-interactive-plots", action="store_true")
    parser.add_argument("--pipeline-version", default="")
    parser.add_argument("--repo-dir", default="")
    parser.add_argument("--run-command", default="")
    parser.add_argument("--profile-stack", default="")
    parser.add_argument("--launch-dir", default="")
    parser.add_argument("--pipeline-outdir", default="")
    parser.add_argument("--container-engine", default="")
    parser.add_argument("--container-image", default="")
    parser.add_argument("--container-digest", default="")
    parser.add_argument("--git-commit", default="")
    parser.add_argument("--git-tag", default="")
    parser.add_argument("--nextflow-session-id", default="")
    parser.add_argument("--nextflow-run-name", default="")
    args = parser.parse_args()
    sample_dir = Path(args.sample_dir)
    out = sample_dir / "panr2_inputs"
    out.mkdir(parents=True, exist_ok=True)

    copy_if_exists(sample_dir / "metadata_output" / "ncbi_clean.csv", out / "metadata" / "ncbi_clean.csv")
    copy_if_exists(sample_dir / "metadata_output" / "ncbi_clean_qc_pass.csv", out / "metadata" / "ncbi_clean_qc_pass.csv")
    copy_if_exists(sample_dir / "metadata_output" / "ncbi_enriched.csv", out / "metadata" / "ncbi_enriched.csv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_clean.csv", out / "metadata" / "fetchm2_clean.csv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_clean.tsv", out / "metadata" / "fetchm2_clean.tsv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_all_assemblies.csv", out / "metadata" / "fetchm2_all_assemblies.csv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_all_assemblies.tsv", out / "metadata" / "fetchm2_all_assemblies.tsv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_clean_compat.csv", out / "metadata" / "fetchm2_clean_compat.csv")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_report.md", out / "metadata" / "fetchm2_report.md")
    copy_if_exists(sample_dir / "metadata_output" / "fetchm2_manifest.json", out / "metadata" / "fetchm2_manifest.json")
    copy_if_exists(sample_dir / "metadata_output" / "sample_map.csv", out / "metadata" / "sample_map.csv")
    copy_if_exists(sample_dir / "metadata_output" / "metadata_completeness.csv", out / "metadata" / "metadata_completeness.csv")
    copy_if_exists(sample_dir / "metadata_output" / "metadata_bias_warning.txt", out / "metadata" / "metadata_bias_warning.txt")
    copy_if_exists(sample_dir / "metadata_output" / "metadata_engine.txt", out / "metadata" / "metadata_engine.txt")
    copy_if_exists(sample_dir / "metadata_analysis" / "metadata_analysis_report.md", out / "metadata_analysis" / "metadata_analysis_report.md")
    copy_if_exists(sample_dir / "metadata_analysis" / "tables" / "field_coverage_summary.csv", out / "metadata_analysis" / "tables" / "field_coverage_summary.csv")
    copy_if_exists(sample_dir / "metadata_analysis" / "tables" / "top_values_by_field.csv", out / "metadata_analysis" / "tables" / "top_values_by_field.csv")
    copy_if_exists(sample_dir / "metadata_analysis" / "tables" / "numeric_summary.csv", out / "metadata_analysis" / "tables" / "numeric_summary.csv")
    copy_if_exists(sample_dir / "audit" / "standardization_summary.csv", out / "metadata_audit" / "standardization_summary.csv")
    copy_if_exists(sample_dir / "audit" / "standardization_audit.md", out / "metadata_audit" / "standardization_audit.md")
    copy_if_exists(sample_dir / "audit" / "production_readiness_gate.md", out / "metadata_audit" / "production_readiness_gate.md")
    copy_if_exists(sample_dir / "audit" / "production_readiness_gate.json", out / "metadata_audit" / "production_readiness_gate.json")
    copy_if_exists(sample_dir / "sequence" / "sequence_download_summary.csv", out / "sequence" / "sequence_download_summary.csv")
    copy_if_exists(sample_dir / "sequence" / "failed_accessions.txt", out / "sequence" / "failed_accessions.txt")
    copy_if_exists(sample_dir / "abricate" / "ncbi_summary.tab", out / "amr" / "ncbi_summary.tab")
    copy_if_exists(sample_dir / "abricate" / "ncbi_results.tab", out / "amr" / "ncbi_results.tab")
    for database_dir in [
        "ncbi",
        "vfdb",
        "plasmidfinder",
        "mobileelementfinder",
        "isfinder",
        "integronfinder",
        "iceberg",
        "mlst",
        "amrfinderplus",
        "defensefinder",
        "prophage",
        "mobsuite",
        "kleborate",
        "kaptive",
        "ectyper",
        "serotypefinder",
        "sccmecfinder",
        "cross_database",
        "temporal",
        "report",
        "panr2_native_feature_runners",
    ]:
        copytree_if_exists(sample_dir / database_dir, out / database_dir)
    copytree_if_exists(sample_dir / "tool_results", out / "tool_results")
    copy_if_exists(
        sample_dir / "panr2_native_feature_runners" / "native_runner_merge_audit.tsv",
        out / "manifest" / "native_runner_merge_audit.tsv",
    )
    copy_if_exists(
        sample_dir / "panr2_native_feature_runners" / "module_status.tsv",
        out / "manifest" / "native_runner_module_status.tsv",
    )
    copy_if_exists(
        sample_dir / "mobsuite" / "mobsuite_database_setup_status.tsv",
        out / "manifest" / "mobsuite_database_setup_status.tsv",
    )
    copy_if_exists(
        sample_dir / "prophage" / "genomad_database_setup_status.tsv",
        out / "manifest" / "genomad_database_setup_status.tsv",
    )
    copy_if_exists(sample_dir / "qc" / "qc_master_report.csv", out / "qc" / "qc_master_report.csv")
    copy_if_exists(sample_dir / "qc" / "excluded_for_panr2.csv", out / "qc" / "excluded_for_panr2.csv")
    copy_if_exists(sample_dir / "ani" / "analysis" / "panr2_ani_summary.csv", out / "ani" / "analysis" / "panr2_ani_summary.csv")
    copy_if_exists(sample_dir / "quast" / "analysis" / "panr2_quast_summary.csv", out / "assembly_qc" / "analysis" / "panr2_quast_summary.csv")

    manifest_rows = []
    candidate_dirs = [
        Path(args.versions_dir) if args.versions_dir else None,
        sample_dir / "pipeline_versions",
        sample_dir.parent / "pipeline_versions",
    ]
    seen_dirs = set()
    for version_dir in [path for path in candidate_dirs if path]:
        version_dir = version_dir.resolve()
        if version_dir in seen_dirs or not version_dir.exists():
            continue
        seen_dirs.add(version_dir)
        for version_file in sorted(version_dir.glob("*_versions.txt")):
            with version_file.open() as handle:
                lines = [line.strip() for line in handle if line.strip()]
            if not lines:
                manifest_rows.append({"component": version_file.stem, "version": "", "source_file": str(version_file)})
                continue
            for line in lines:
                parsed = parse_version_line(line, version_file.stem)
                if parsed:
                    parsed["source_file"] = str(version_file)
                    manifest_rows.append(parsed)
    manifest_path = out / "manifest" / "software_versions.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["component", "version", "source_file"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    schema_path = out / "manifest" / "panr2_feature_contract_columns.txt"
    schema_path.write_text("\n".join(CONTRACT_COLUMNS) + "\n")
    all_schema_path = out / "manifest" / "panr2_feature_contract_all_columns.txt"
    all_schema_path.write_text("\n".join(FEATURE_COLUMNS) + "\n")
    export_contract(
        sample_dir,
        out,
        large_dataset=args.large_dataset,
        report_mode=args.report_mode,
        max_features_heatmap=args.max_features_heatmap,
        max_features_network=args.max_features_network,
        max_metadata_columns=args.max_metadata_columns,
        top_n_features_per_database=args.top_n_features_per_database,
        skip_heavy_interactive_plots=args.skip_heavy_interactive_plots,
    )
    write_reproducibility_manifest(sample_dir, out, args)


if __name__ == "__main__":
    main()
