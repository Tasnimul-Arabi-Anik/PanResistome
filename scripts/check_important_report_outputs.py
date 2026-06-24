#!/usr/bin/env python3
"""Validate the PanResistome important report shell and linked assets."""

from __future__ import annotations

import argparse
import csv
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
    "report_highlights.tsv",
    "report_highlights_by_section.tsv",
    "warning_priority_summary.tsv",
    "report_visual_index.tsv",
    "report_visual_quality.tsv",
    "finding_confidence_summary.tsv",
    "warnings_and_limitations_summary.tsv",
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
    "report_highlights.tsv",
    "report_highlights_by_section.tsv",
    "warning_priority_summary.tsv",
    "report_visual_index.tsv",
    "report_visual_quality.tsv",
    "warnings_and_limitations_summary.tsv",
    "warnings_and_limitations.tsv",
    "warnings_by_section.tsv",
    "module_warning_summary.tsv",
    "report_cap_summary.tsv",
    "important_file_index.tsv",
    "download_manifest.tsv",
]

REQUIRED_DOWNLOADS = [
    "important_summary_tables.zip",
    "important_tables.zip",
    "important_figures.zip",
    "publication_candidate_figures.zip",
    "important_report_assets.zip",
]

GENERIC_CAPTION_PHRASES = [
    "Report-facing figure with PNG",
    "Report-facing figure with companion downloads",
    "Figure preview with companion PNG",
]

PLACEHOLDER_YEAR_VALUES = {"", "0", "1", "0001", "1900", "1-01-01", "0001-01-01"}
INTERPRETATION_SENSITIVE_SECTIONS = {
    "Geographic Distribution",
    "Co-occurrence / Genomic Context",
    "Lineage / Clonal Structure",
    "Metadata Associations",
    "Temporal Trends",
}


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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8", errors="ignore") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _to_float(value: str) -> float | None:
    try:
        text = str(value or "").strip()
        if not text:
            return None
        return float(text)
    except ValueError:
        return None


def _normalize_feature_id(value: str) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"\s+", "", text)
    text = re.sub(r"[^a-z0-9]+", "", text)
    return text


def _invalid_temporal_value(value: str) -> bool:
    text = str(value or "").strip().lower()
    if text in {"not applicable", "unknown", "missing", "na", "n/a", "none"}:
        return True
    return text in PLACEHOLDER_YEAR_VALUES


def _check_prevalence_consistency(important_dir: Path, errors: list[str]) -> None:
    rows = _read_tsv(important_dir / "tables" / "prevalence_summary_by_database.tsv")
    for row in rows:
        database = row.get("database", "unknown")
        db_positive = _to_float(row.get("positive_genomes", "")) or 0.0
        top_positive = _to_float(row.get("top_feature_positive_genomes", "")) or 0.0
        db_percent = _to_float(row.get("genomes_positive_percent", "")) or 0.0
        top_percent = _to_float(row.get("top_feature_prevalence_percent", "")) or 0.0
        if db_positive + 1e-9 < top_positive:
            errors.append(
                f"Database prevalence inconsistency for {database}: positive_genomes {db_positive:g} < top_feature_positive_genomes {top_positive:g}"
            )
        if db_percent + 1e-9 < top_percent:
            errors.append(
                f"Database prevalence inconsistency for {database}: genomes_positive_percent {db_percent:g} < top_feature_prevalence_percent {top_percent:g}"
            )


def _check_temporal_placeholders(important_dir: Path, errors: list[str]) -> None:
    temporal_tables = [
        important_dir / "key_tables" / "temporal_trend_summary.tsv",
        important_dir / "key_tables" / "temporal_feature_prevalence.tsv",
        important_dir / "key_tables" / "temporal_database_burden.tsv",
    ]
    year_fields = ["collection_year", "first_year", "last_year"]
    for table_path in temporal_tables:
        for row_index, row in enumerate(_read_tsv(table_path), 2):
            for field in year_fields:
                if field in row and _invalid_temporal_value(row.get(field, "")):
                    errors.append(
                        f"Temporal table {table_path.relative_to(important_dir)} contains placeholder {field}={row.get(field, '')!r} on row {row_index}"
                    )
                    break


def _check_highlight_quality(important_dir: Path, errors: list[str]) -> None:
    highlights = _read_tsv(important_dir / "tables" / "report_highlights.tsv")
    for row_index, row in enumerate(highlights, 2):
        if row.get("highlight_type") != "informative_cooccurrence":
            continue
        primary = row.get("primary_feature", "")
        secondary = row.get("secondary_feature", "")
        if primary and secondary and _normalize_feature_id(primary) == _normalize_feature_id(secondary):
            errors.append(
                f"Co-occurrence highlight row {row_index} uses identical normalized features: {primary!r} and {secondary!r}"
            )

    by_section = _read_tsv(important_dir / "tables" / "report_highlights_by_section.tsv")
    sections = {row.get("section", "") for row in by_section if row.get("section", "")}
    if len(by_section) >= 25 and len(sections) < 5:
        errors.append(
            f"Balanced highlights by section contains only {len(sections)} distinct sections; expected at least 5 when enough rows are available"
        )


def _check_visual_quality(important_dir: Path, errors: list[str]) -> None:
    rows = _read_tsv(important_dir / "tables" / "report_visual_quality.tsv")
    required_columns = {
        "asset_quality_label",
        "render_quality_label",
        "axis_label_status",
        "interpretation_quality_label",
        "title_quality_label",
        "caption_quality_label",
        "final_publication_label",
        "default_visibility",
        "publication_candidate",
    }
    if rows:
        missing = required_columns - set(rows[0])
        if missing:
            errors.append(f"Report visual quality table is missing columns: {', '.join(sorted(missing))}")
    for row in rows:
        section = row.get("section", "")
        label = row.get("final_publication_label", "")
        if section in INTERPRETATION_SENSITIVE_SECTIONS and label == "publication_ready":
            errors.append(f"Interpretation-sensitive figure marked publication_ready: {row.get('figure_stem', '')}")
        render_quality = row.get("render_quality_label", "rendered")
        axis_status = row.get("axis_label_status", "not_applicable")
        default_visibility = row.get("default_visibility", "")
        if render_quality != "rendered" and default_visibility in {"featured", "standard"}:
            errors.append(
                f"Figure with render_quality_label={render_quality} promoted to {default_visibility}: {row.get('figure_stem', '')}"
            )
        if render_quality != "rendered" and row.get("publication_candidate", "").lower() == "true":
            errors.append(
                f"Figure with render_quality_label={render_quality} marked publication_candidate: {row.get('figure_stem', '')}"
            )
        if axis_status == "missing_axis_labels" and default_visibility in {"featured", "standard"}:
            errors.append(f"Figure missing axis labels promoted to {default_visibility}: {row.get('figure_stem', '')}")

    for row in rows:
        stem = row.get("figure_stem", "")
        render_quality = row.get("render_quality_label", "rendered")
        if stem.startswith("variation_identity_coverage_") and render_quality == "rendered":
            data_rows = _read_tsv(important_dir / "figures" / f"{stem}.data.tsv")
            if data_rows and not any(
                _to_float(item.get("identity", "")) is not None and _to_float(item.get("coverage", "")) is not None
                for item in data_rows
            ):
                errors.append(f"Variation scatter marked rendered without numeric identity/coverage pairs: {stem}")
        if stem.startswith("cooccurrence_network_") and render_quality == "rendered":
            data_rows = _read_tsv(important_dir / "figures" / f"{stem}.data.tsv")
            if len(data_rows) == 0:
                errors.append(f"Co-occurrence network marked rendered with no edges: {stem}")

    for svg_path in sorted((important_dir / "figures").glob("feature_profile_pcoa_*.svg")):
        svg_text = svg_path.read_text(encoding="utf-8", errors="ignore")
        if "PCoA1" not in svg_text or "PCoA2" not in svg_text:
            errors.append(f"PCoA SVG lacks explicit axis labels: figures/{svg_path.name}")


def validate(sample_dir: Path) -> list[str]:
    root_dir, important_dir = _resolve_paths(sample_dir)
    errors: list[str] = []
    results_html = important_dir / "results.html"
    html_text = results_html.read_text(encoding="utf-8", errors="ignore")

    if "file://" in html_text:
        errors.append("results.html contains file:// links")
    for phrase in GENERIC_CAPTION_PHRASES:
        if phrase in html_text:
            errors.append(f"results.html contains generic caption phrase: {phrase}")
    if "Warning rows represent" in html_text:
        errors.append("results.html contains outdated warning-row percentage wording")

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
        "sidebar-links",
        "figure-card",
        "table-card",
        "analysis-card",
        "warning-box",
        "download-card",
        "back-to-top",
    ]:
        if class_name not in html_text:
            errors.append(f"Missing report UI class: {class_name}")
    if "@media (max-width: 920px)" not in html_text or "grid-template-columns: repeat(auto-fit, minmax(150px, 1fr))" not in html_text:
        errors.append("results.html is missing responsive sidebar grid CSS")

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

    _check_prevalence_consistency(important_dir, errors)
    _check_temporal_placeholders(important_dir, errors)
    _check_highlight_quality(important_dir, errors)
    _check_visual_quality(important_dir, errors)

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
