# PanR2 Input Contract

PanResistome is the preferred execution layer for heavy genomics tools. PanR2 is the lightweight analysis and reporting layer. Every PanResistome tool module that produces feature-like annotations should export a PanR2-compatible table whenever possible.

## Feature Table Schema

Required columns:

| Column | Meaning |
| --- | --- |
| `sample_id` | Local sample name, filename stem, or user-provided sample identifier |
| `assembly_accession` | Stable assembly accession or normalized genome ID used for merging |
| `database` | Feature family or database name, for example `amr`, `amrfinderplus`, `vfdb`, `plasmidfinder`, `mobileelementfinder`, `integronfinder`, `mlst`, `defensefinder`, `mobsuite`, `prophage`, `kleborate`, `kaptive`, `ectyper`, `serotypefinder`, `sccmecfinder`, `ani`, or `assembly_qc` |
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

Optional standardized columns may also be present and should be preserved by downstream tools when possible:

```text
feature_name
feature_description
feature_subcategory
mechanism
drug_class
product
sequence_id
strand
source_table
source_file
source_database
raw_feature_id
raw_category
raw_method
evidence_type
confidence
notes
```

## Current PanResistome Exports

PanResistome writes `panr2_inputs/` after annotation when `--export_panr2_inputs true`. FetchM2 is the default metadata source, but PanResistome also writes legacy-compatible `ncbi_clean.csv` files so PanR2 and older downstream scripts can read the same run.

For external table inputs whose sample names do not contain `GCF_`/`GCA_` assembly accessions, pass `--panr2_sample_map path/to/sample_map.csv` to PanResistome. The map must contain `sample_id` and `Assembly Accession` columns. When FetchM2 writes `metadata_output/sample_map.csv`, PanResistome passes it to PanR2 automatically.

```text
panr2_inputs/
├── metadata/ncbi_clean.csv
├── metadata/ncbi_clean_qc_pass.csv
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
├── features/
│   ├── all_features.tsv
│   ├── amr.features.tsv
│   ├── amrfinderplus.features.tsv
│   └── <database>.features.tsv
├── feature_matrices/
│   ├── all_features_presence_absence.tsv
│   └── <database>_presence_absence.tsv
├── metadata_feature_analysis/
│   ├── fetchm2_metadata_audit.tsv
│   ├── metadata_column_eligibility.tsv
│   ├── metadata_normalized_for_analysis.tsv
│   ├── feature_eligibility.tsv
│   ├── database_burden_by_sample.tsv
│   ├── top_findings.tsv
│   ├── top_findings.md
│   └── prevalence_tables/
├── cross_database/
│   ├── feature_cooccurrence.tsv
│   ├── database_cooccurrence_summary.tsv
│   ├── amr_mge_context.tsv
│   ├── amr_plasmid_context.tsv
│   └── amrfinder_abricate_concordance.tsv
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
    ├── panr2_feature_contract_columns.txt
    ├── panr2_feature_contract_all_columns.txt
    ├── schema_validation_report.csv
    ├── schema_validation_summary.txt
    ├── feature_completeness_audit.tsv
    ├── module_status_summary.tsv
    ├── invalid_feature_rows.csv
    ├── duplicate_features.csv
    └── unmatched_features.csv
```

`panr2_inputs/features/*.features.tsv` is the strict contract layer. Raw tool folders are still copied for traceability, but downstream analysis should prefer the standardized feature tables when possible. `schema_validation_report.csv` checks required columns, and `unmatched_features.csv` lists feature rows whose assembly accession cannot be matched to metadata.

## Design Rule

Do not make PanR2 install every annotation tool. Add heavy runners and database setup to PanResistome, then export standardized tables for PanR2 to analyze.
