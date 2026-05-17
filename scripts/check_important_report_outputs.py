#!/usr/bin/env python3
"""Validate the PanResistome important report shell and linked assets."""

from __future__ import annotations

import argparse
import re
import sys
import zipfile
from pathlib import Path


REQUIRED_SECTION_IDS = [
    "featured",
    "overview",
    "qc",
    "enriched-dataset",
    "prevalence",
    "geography",
    "variations",
    "temporal",
    "cooccurrence",
    "metadata-associations",
    "lineage",
    "diversity",
    "notable-genomes",
    "ordination",
    "concordance",
    "evidence",
    "warnings",
    "downloads",
]

REQUIRED_TABLES = [
    "feature_prevalence.tsv",
    "geographic_distribution_summary.tsv",
    "feature_variation_summary.tsv",
    "temporal_trend_summary.tsv",
    "cooccurrence_pair_summary.tsv",
    "metadata_feature_enrichment.tsv",
    "lineage_summary.tsv",
    "diversity_report_summary.tsv",
    "notable_genomes.tsv",
    "finding_confidence_summary.tsv",
    "warnings_and_limitations.tsv",
    "download_manifest.tsv",
]

REQUIRED_FINAL_TABLES = [
    "notable_genomes.tsv",
    "notable_genome_score_components.tsv",
    "feature_profile_ordination.tsv",
    "database_concordance_summary.tsv",
    "amr_concordance_feature_level.tsv",
    "amr_concordance_by_sample.tsv",
    "evidence_summary.tsv",
    "finding_confidence_summary.tsv",
    "evidence_by_section.tsv",
    "warnings_and_limitations.tsv",
    "warnings_by_section.tsv",
    "module_warning_summary.tsv",
    "report_cap_summary.tsv",
    "important_file_index.tsv",
    "download_manifest.tsv",
]

REQUIRED_DOWNLOADS = [
    "important_tables.zip",
    "important_figures.zip",
    "important_report_assets.zip",
]


def _resolve_paths(input_path: Path) -> tuple[Path, Path]:
    path = input_path.resolve()
    if path.name == "important" and (path / "results.html").exists():
        return path.parent, path
    if (path / "important" / "results.html").exists():
        return path, path / "important"
    if path.name == "results.html" and path.exists():
        return path.parent.parent, path.parent
    raise SystemExit(f"Could not find important/results.html under {input_path}")


def _is_local_link(link: str) -> bool:
    if not link or link.startswith("#"):
        return False
    lowered = link.lower()
    return not (
        lowered.startswith("http://")
        or lowered.startswith("https://")
        or lowered.startswith("mailto:")
        or lowered.startswith("javascript:")
        or lowered.startswith("data:")
    )


def _linked_files(html_text: str) -> set[str]:
    links: set[str] = set()
    for match in re.finditer(r"""(?:href|src)=["']([^"']+)["']""", html_text):
        link = match.group(1).split("#", 1)[0]
        if _is_local_link(link):
            links.add(link)
    return links


def _table_exists(important_dir: Path, table_name: str) -> bool:
    return (
        (important_dir / "tables" / table_name).exists()
        or (important_dir / "key_tables" / table_name).exists()
    )


def _zip_non_empty(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            return bool(archive.namelist())
    except zipfile.BadZipFile:
        return False


def validate(sample_dir: Path) -> list[str]:
    root_dir, important_dir = _resolve_paths(sample_dir)
    errors: list[str] = []
    results_html = important_dir / "results.html"
    html_text = results_html.read_text(encoding="utf-8", errors="ignore")

    required_dirs = [
        root_dir / "basic",
        important_dir / "tables",
        important_dir / "figures",
        important_dir / "downloads",
    ]
    for directory in required_dirs:
        if not directory.exists():
            errors.append(f"Missing required directory: {directory.relative_to(root_dir)}")

    if not (root_dir / "basic" / "enriched_genome_dataset.csv").exists():
        errors.append("Missing basic/enriched_genome_dataset.csv")

    for section_id in REQUIRED_SECTION_IDS:
        if f'id="{section_id}"' not in html_text and f"id='{section_id}'" not in html_text:
            errors.append(f"Missing report section anchor: {section_id}")

    for class_name in [
        "report-header",
        "sidebar",
        "figure-card",
        "table-card",
        "warning-box",
        "download-card",
        "back-to-top",
    ]:
        if class_name not in html_text:
            errors.append(f"Missing report UI class: {class_name}")

    for link in sorted(_linked_files(html_text)):
        target = (important_dir / link).resolve()
        try:
            target.relative_to(root_dir.resolve())
        except ValueError:
            errors.append(f"Linked file escapes output directory: {link}")
            continue
        if not target.exists():
            errors.append(f"Broken relative link: {link}")

    for table_name in REQUIRED_TABLES:
        if not _table_exists(important_dir, table_name):
            errors.append(f"Missing required report table: {table_name}")

    for table_name in REQUIRED_FINAL_TABLES:
        if not (important_dir / "tables" / table_name).exists():
            errors.append(f"Missing required final interpretation table: {table_name}")

    for zip_name in REQUIRED_DOWNLOADS:
        zip_path = important_dir / "downloads" / zip_name
        if not _zip_non_empty(zip_path):
            errors.append(f"Missing or empty download ZIP: downloads/{zip_name}")

    report_caps = important_dir / "tables" / "report_cap_summary.tsv"
    if not report_caps.exists() or report_caps.stat().st_size == 0:
        errors.append("Missing report cap summary; large-mode capping cannot be audited")

    complete_candidates = [
        root_dir / "panr2_inputs" / "features" / "all_features.tsv",
        root_dir / "panr2_inputs" / "manifest" / "feature_contract.json",
        root_dir / "panr2_inputs" / "manifest" / "schema_validation_summary.txt",
    ]
    for path in complete_candidates:
        if not path.exists():
            errors.append(f"Missing preserved complete output: {path.relative_to(root_dir)}")

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_dir", type=Path, help="Sample output directory, important directory, or results.html path")
    args = parser.parse_args(argv)
    errors = validate(args.sample_dir)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("important report QA passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
