#!/usr/bin/env python3
import argparse
import csv
import re
import shutil
from pathlib import Path

from panr2_contract import CONTRACT_COLUMNS, FEATURE_COLUMNS, export_contract


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


def main():
    parser = argparse.ArgumentParser(description="Export a PanResistome sample directory into PanR2-ready inputs.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--versions-dir", help="Optional pipeline_versions directory to include in the PanR2 manifest.")
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
    ]:
        copytree_if_exists(sample_dir / database_dir, out / database_dir)
    copytree_if_exists(sample_dir / "tool_results", out / "tool_results")
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
    export_contract(sample_dir, out)


if __name__ == "__main__":
    main()
