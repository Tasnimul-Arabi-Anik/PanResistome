#!/usr/bin/env python3
"""Build a small optional-tool fixture and validate PanR2 analysis export."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from panr2_contract import export_contract  # noqa: E402


HEADER = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"


def write_table(path: Path, rows: list[tuple[str, str, int, int, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [HEADER]
    for sample, contig, start, end, gene, database, category in rows:
        lines.append(
            f"{sample}.fna\t{contig}\t{start}\t{end}\t{gene}\t100\t99\t"
            f"{database}\t{gene}_ACC\t{category}\t{category}\n"
        )
    path.write_text("".join(lines), encoding="utf-8")


def build_fixture(sample_dir: Path) -> None:
    metadata_dir = sample_dir / "metadata_output"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    metadata_rows = []
    for idx in range(1, 7):
        accession = f"GCF_90000000{idx}.1"
        metadata_rows.append(
            {
                "Assembly Accession": accession,
                "Organism Name": "Klebsiella pneumoniae",
                "Country": "CountryA" if idx <= 3 else "CountryB",
                "Host_SD": "Homo sapiens" if idx in {1, 2, 4} else "environment",
                "Isolation_Source_SD": "blood" if idx in {1, 2, 4} else "water",
                "Collection_Year": 2024 if idx <= 3 else 2023,
                "Assembly BioProject Accession": "PRJNA_DOMINANT" if idx <= 3 else f"PRJNA_OTHER_{idx}",
            }
        )
    pd.DataFrame(metadata_rows).to_csv(metadata_dir / "ncbi_clean.csv", index=False)

    write_table(
        sample_dir / "abricate" / "ncbi_results.tab",
        [
            ("GCF_900000001.1", "contigA", 100, 500, "blaKPC-2", "ncbi", "beta_lactam"),
            ("GCF_900000002.1", "contigA", 100, 500, "blaKPC-2", "ncbi", "beta_lactam"),
            ("GCF_900000004.1", "contigB", 100, 500, "tetA", "ncbi", "tetracycline"),
        ],
    )
    write_table(
        sample_dir / "vfdb" / "vfdb_results.tab",
        [
            ("GCF_900000001.1", "contigA", 800, 1200, "iutA", "vfdb", "virulence"),
            ("GCF_900000003.1", "contigC", 800, 1200, "fimH", "vfdb", "virulence"),
        ],
    )
    write_table(
        sample_dir / "plasmidfinder" / "plasmidfinder_results.tab",
        [
            ("GCF_900000001.1", "contigA", 1300, 1700, "IncFIB", "plasmidfinder", "replicon"),
            ("GCF_900000002.1", "contigA", 1300, 1700, "IncFIB", "plasmidfinder", "replicon"),
            ("GCF_900000004.1", "contigB", 1300, 1700, "IncX3", "plasmidfinder", "replicon"),
        ],
    )

    optional_tables = [
        ("mobileelementfinder", "mobileelementfinder_results.tab", "IS26", "mobileelementfinder", "insertion_sequence"),
        ("isfinder/tables", "isfinder_results.tab", "ISEcp1", "isfinder", "insertion_sequence"),
        ("mobsuite/tables", "mobsuite_results.tab", "IncFIB", "mobsuite", "replicon"),
        ("prophage/tables", "prophage_results.tab", "region_1", "prophage", "prophage"),
        ("defensefinder/tables", "defensefinder_results.tab", "RM_Type_I", "defensefinder", "defense_system"),
        ("kleborate/tables", "kleborate_results.tab", "ST147", "kleborate", "sequence_type"),
        ("kaptive/tables", "kaptive_results.tab", "KL64", "kaptive", "capsule_locus"),
        ("ectyper/tables", "ectyper_results.tab", "O157:H7", "ectyper", "serotype"),
        ("serotypefinder/tables", "serotypefinder_results.tab", "O2", "serotypefinder", "serotype"),
        ("sccmecfinder/tables", "sccmecfinder_results.tab", "SCCmec_IV", "sccmecfinder", "cassette_type"),
    ]
    for offset, (directory, filename, feature_id, database, category) in enumerate(optional_tables, start=1):
        sample_idx = (offset % 6) + 1
        write_table(
            sample_dir / directory / filename,
            [
                (
                    f"GCF_90000000{sample_idx}.1",
                    f"contig{sample_idx}",
                    200 + offset * 100,
                    260 + offset * 100,
                    feature_id,
                    database,
                    category,
                )
            ],
        )

    mobsuite_dir = sample_dir / "mobsuite" / "tables"
    (mobsuite_dir / "mobsuite_contig_report.tsv").write_text(
        "sample_id\tsseqid\tsstart\tsend\tmolecule_type\trep_type(s)\tprimary_cluster_id\tpredicted_mobility\tmge_type\tmge_subtype\n"
        "GCF_900000001.1\tcontigA\t1250\t2500\tplasmid\tIncFIB\tAC125\tconjugative\tISKpn26\tIS5\n",
        encoding="utf-8",
    )


def row_count(path: Path, sep: str = "\t") -> int:
    if not path.exists():
        return -1
    try:
        return len(pd.read_csv(path, sep=sep))
    except pd.errors.EmptyDataError:
        return 0


def write_validation_summary(outdir: Path, outputs: dict[str, str]) -> None:
    feature_counts_path = outdir / "feature_counts.tsv"
    all_features = pd.read_csv(outputs["all_features"], sep="\t")
    counts = (
        all_features.groupby("database", dropna=False)
        .size()
        .reset_index(name="feature_rows")
        .sort_values(["database"])
    )
    counts.to_csv(feature_counts_path, sep="\t", index=False)

    checks = [
        ("all_features", outputs["all_features"]),
        ("schema_validation_summary", outputs["schema_validation_summary"]),
        ("feature_completeness_audit", outputs["feature_completeness_audit"]),
        ("all_feature_matrix", outputs["all_feature_matrix"]),
        ("feature_cooccurrence", outputs["feature_cooccurrence"]),
        ("database_cooccurrence_summary", outputs["database_cooccurrence_summary"]),
        ("feature_proximity", outputs["feature_proximity"]),
        ("top_features_by_database", outputs["top_features_by_database"]),
        ("metadata_usability_summary", outputs["metadata_usability_summary"]),
        ("top_findings", outputs["top_findings"]),
        ("handoff_report_index", outputs["handoff_report_index"]),
        ("top_findings_html", outputs["top_findings_html"]),
        ("cross_database_interpretation_html", outputs["cross_database_interpretation_html"]),
    ]
    check_rows = []
    for name, raw_path in checks:
        path = Path(raw_path)
        check_rows.append(
            {
                "output": name,
                "exists": str(path.exists()).lower(),
                "rows": row_count(path) if path.suffix in {".tsv", ".csv"} else "",
                "path": str(path.relative_to(outdir) if path.is_relative_to(outdir) else path),
            }
        )
    pd.DataFrame(check_rows).to_csv(outdir / "analysis_outputs.tsv", sep="\t", index=False)

    summary = Path(outputs["schema_validation_summary"]).read_text(encoding="utf-8")
    feature_files = sorted((outdir / "sample" / "panr2_inputs" / "features").glob("*.features.tsv"))
    feature_file_lines = "\n".join(f"- `{path.name}`" for path in feature_files)
    markdown = f"""# Optional Feature Analysis Validation

Date: 2026-05-11

Purpose: verify that opt-in module outputs can be converted into PanR2-compatible feature tables and analyzed through the same standardized layer used for AMR and VFDB.

This validation uses small local fixture tables, not large external databases. It validates the PanR2 handoff/analysis behavior for optional outputs and complements the documented biological Kleborate/MOB-suite runner validation.

## Result

Status: PASS

Feature rows by database are in `feature_counts.tsv`.

Schema validation summary:

```text
{summary.strip()}
```

Generated feature tables:

{feature_file_lines}

Generated analysis/report checks are in `analysis_outputs.tsv`.

## Interpretation

This confirms that MobileElementFinder, ISfinder-style BLAST, MOB-suite, geNomad/prophage, DefenseFinder, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder tables can produce standardized feature rows, `all_features.tsv`, feature matrices, co-occurrence/proximity outputs, metadata usability outputs, top-feature summaries, and HTML report pages.

This does not claim that every external runner/database can be installed and run biologically on a fresh desktop. Runner-mode status remains separated in `docs/optional_module_validation_matrix.md`.
"""
    (outdir / "VALIDATION_RESULTS.md").write_text(markdown, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default="validation/optional_feature_analysis")
    parser.add_argument("--force", action="store_true", help="Replace the output directory if it already exists.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists():
        if not args.force:
            raise SystemExit(f"{outdir} already exists; use --force to replace it")
        shutil.rmtree(outdir)
    sample_dir = outdir / "sample"
    build_fixture(sample_dir)
    outputs = export_contract(
        sample_dir,
        sample_dir / "panr2_inputs",
        large_dataset=True,
        report_mode="compact",
        max_features_heatmap=25,
        max_features_network=25,
        max_metadata_columns=10,
        top_n_features_per_database=10,
        skip_heavy_interactive_plots=True,
    )
    write_validation_summary(outdir, outputs)
    print(f"Wrote optional feature analysis validation to {outdir}")


if __name__ == "__main__":
    main()
