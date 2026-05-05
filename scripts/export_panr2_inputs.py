#!/usr/bin/env python3
import argparse
import csv
import shutil
from pathlib import Path


CONTRACT_COLUMNS = [
    "sample_id", "assembly_accession", "database", "feature_id", "feature_category",
    "presence", "identity", "coverage", "contig", "start", "end", "tool",
    "tool_version", "database_version",
]


def copy_if_exists(src, dst):
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main():
    parser = argparse.ArgumentParser(description="Export a PanResistome sample directory into PanR2-ready inputs.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--versions-dir", help="Optional pipeline_versions directory to include in the PanR2 manifest.")
    args = parser.parse_args()
    sample_dir = Path(args.sample_dir)
    out = sample_dir / "panr2_inputs"
    out.mkdir(parents=True, exist_ok=True)

    copy_if_exists(sample_dir / "metadata_output" / "ncbi_clean.csv", out / "metadata" / "ncbi_clean.csv")
    copy_if_exists(sample_dir / "metadata_output" / "ncbi_enriched.csv", out / "metadata" / "ncbi_enriched.csv")
    copy_if_exists(sample_dir / "abricate" / "ncbi_summary.tab", out / "amr" / "ncbi_summary.tab")
    copy_if_exists(sample_dir / "abricate" / "ncbi_results.tab", out / "amr" / "ncbi_results.tab")
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
                if "==" in line:
                    component, version = line.split("==", 1)
                elif "=" in line:
                    component, version = line.split("=", 1)
                else:
                    component, version = version_file.stem, line
                manifest_rows.append({"component": component.strip(), "version": version.strip(), "source_file": str(version_file)})
    manifest_path = out / "manifest" / "software_versions.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["component", "version", "source_file"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    schema_path = out / "manifest" / "panr2_feature_contract_columns.txt"
    schema_path.write_text("\n".join(CONTRACT_COLUMNS) + "\n")


if __name__ == "__main__":
    main()
