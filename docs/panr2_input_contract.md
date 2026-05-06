# PanR2 Input Contract

PanResistome is the preferred execution layer for heavy genomics tools. PanR2 is the lightweight analysis and reporting layer. Every PanResistome tool module that produces feature-like annotations should export a PanR2-compatible table whenever possible.

## Feature Table Schema

Required columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Local sample name, filename stem, or user-provided sample identifier |
| `assembly_accession` | Stable assembly accession or normalized genome ID used for merging |
| `database` | Feature family or database name, for example `amr`, `vfdb`, `plasmidfinder`, `mobileelementfinder`, `integronfinder`, `mlst`, `defensefinder`, `mobsuite`, `prophage`, `kleborate`, `kaptive`, `ectyper`, `serotypefinder`, `sccmecfinder`, `ani`, or `assembly_qc` |
| `feature_id` | Gene, replicon, system, type, cluster, or feature identifier |
| `feature_category` | Higher-level category where available |
| `presence` | `1` for present, `0` for absent |
| `identity` | Identity, ANI, score, N50, or comparable primary numeric value where relevant |
| `coverage` | Coverage, assembly length, or comparable secondary numeric value where relevant |
| `contig` | Contig or sequence identifier where relevant |
| `start` | Start coordinate where relevant |
| `end` | End coordinate where relevant |
| `tool` | Tool that produced the feature |
| `tool_version` | Tool version if available |
| `database_version` | Database version/date if available |

## Current PanResistome Exports

PanResistome writes `panr2_inputs/` after annotation when `--export_panr2_inputs true`. FetchM2 is the default metadata source, but PanResistome also writes legacy-compatible `ncbi_clean.csv` files so PanR2 and older downstream scripts can read the same run.

For external table inputs whose sample names do not contain `GCF_`/`GCA_` assembly accessions, pass `--panr2_sample_map path/to/sample_map.csv` to PanResistome. The map must contain `sample_id` and `Assembly Accession` columns. When FetchM2 writes `metadata_output/sample_map.csv`, PanResistome passes it to PanR2 automatically.

```text
panr2_inputs/
├── metadata/ncbi_clean.csv
├── metadata/fetchm2_clean.csv
├── metadata/fetchm2_clean.tsv
├── metadata/fetchm2_all_assemblies.csv
├── metadata/fetchm2_all_assemblies.tsv
├── metadata/fetchm2_clean_compat.csv
├── metadata/fetchm2_report.md
├── metadata/fetchm2_manifest.json
├── metadata/sample_map.csv
├── metadata/metadata_completeness.csv
├── metadata/metadata_bias_warning.txt
├── metadata/metadata_engine.txt
├── metadata/ncbi_enriched.csv
├── metadata_analysis/metadata_analysis_report.md
├── metadata_audit/standardization_summary.csv
├── metadata_audit/standardization_audit.md
├── metadata_audit/production_readiness_gate.md
├── sequence/sequence_download_summary.csv
├── sequence/failed_accessions.txt
├── amr/ncbi_summary.tab
├── amr/ncbi_results.tab
├── ani/analysis/panr2_ani_summary.csv
├── assembly_qc/analysis/panr2_quast_summary.csv
├── vfdb/
├── plasmidfinder/
├── mobileelementfinder/
├── integronfinder/
├── mlst/
├── defensefinder/
├── mobsuite/
├── prophage/
├── kleborate/
├── kaptive/
├── ectyper/
├── serotypefinder/
├── sccmecfinder/
├── cross_database/
├── report/
├── qc/qc_master_report.csv
├── qc/excluded_for_panr2.csv
└── manifest/
    ├── software_versions.csv
    └── panr2_feature_contract_columns.txt
```

## Design Rule

Do not make PanR2 install every annotation tool. Add heavy runners and database setup to PanResistome, then export standardized tables for PanR2 to analyze.
