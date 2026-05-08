# Delftia tsuruhatensis Validation Results

Validated on: 2026-05-08

This validation used the reproducible NCBI Assembly input in this directory:

```text
validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv
```

The input contains 45 current NCBI Assembly records for `Delftia tsuruhatensis`.

## Command

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_current \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --checkm2_db /home/anik/genomics/Tools_Dev/PanResistome/results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
  --threads 4 \
  --fetchm2_download_workers 2 \
  -resume
```

GTDB-Tk was intentionally disabled because of its heavy external reference database.

## Runtime

The final validated rerun completed successfully:

```text
Duration: 1h 57m 17s
CPU hours: 14.0
Succeeded: 10
Cached: 9
```

An earlier full non-cached comprehensive validation completed in approximately 2h55m on the same 45-genome dataset.

## Modules Validated

The following stages completed successfully:

- FetchM2 metadata standardization and sequence download
- sequence QC
- CheckM2 quality assessment
- QUAST assembly statistics
- ANI/skani analysis
- Mash pre-screening
- combined QC decision engine
- AMRFinderPlus database auto-download and indexing
- AMRFinderPlus annotation
- PanR2 comprehensive analysis
- PanR2 handoff export
- final result collection

AMRFinderPlus downloaded and indexed database version `2026-03-24.1` during validation. All 45 AMRFinderPlus sample runs passed.

## QC Summary

All 45 assemblies passed the validated QC gates:

```text
sequence_qc_status: PASS = 45
checkm2_qc_status: PASS = 45
combined_qc_status: PASS = 45
ani_species_consistency_status: PASS = 45
qc_master_status: PASS = 45
```

The main QC report is generated at:

```text
validation_runs/delftia_current/Delftia_tsuruhatensis/qc/qc_master_report.csv
```

## FetchM2 Metadata Coverage

FetchM2 produced rich standardized metadata. Non-missing values in the 45 validated assemblies included:

```text
Country: 38/45
Collection_Year: 42/45
Host_SD: 11/45
Sample_Type_SD: 25/45
Environment_Medium_SD: 10/45
Isolation_Source_SD: 2/45
```

Metadata completeness and bias reports are generated under:

```text
validation_runs/delftia_current/Delftia_tsuruhatensis/qc/
validation_runs/delftia_current/Delftia_tsuruhatensis/metadata_analysis/
```

## PanR2 Feature Contract Validation

The PanR2 feature handoff passed schema validation:

```text
feature_files_checked=6
feature_rows=291
databases_seen=amr,amrfinderplus,integronfinder,mlst,plasmidfinder,vfdb
metadata_accessions=45
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature rows exported:

| Feature table | Rows |
| --- | ---: |
| `all_features.tsv` | 291 |
| `amr.features.tsv` | 19 |
| `amrfinderplus.features.tsv` | 24 |
| `vfdb.features.tsv` | 40 |
| `plasmidfinder.features.tsv` | 2 |
| `integronfinder.features.tsv` | 116 |
| `mlst.features.tsv` | 90 |

The feature audit is generated at:

```text
validation_runs/delftia_current/Delftia_tsuruhatensis/panr2_inputs/manifest/feature_completeness_audit.tsv
```

## Main Outputs

The main dashboard is:

```text
validation_runs/delftia_current/Delftia_tsuruhatensis/report/index.html
```

Important output folders:

```text
validation_runs/delftia_current/Delftia_tsuruhatensis/ncbi/
validation_runs/delftia_current/Delftia_tsuruhatensis/vfdb/
validation_runs/delftia_current/Delftia_tsuruhatensis/plasmidfinder/
validation_runs/delftia_current/Delftia_tsuruhatensis/integronfinder/
validation_runs/delftia_current/Delftia_tsuruhatensis/mlst/
validation_runs/delftia_current/Delftia_tsuruhatensis/cross_database/
validation_runs/delftia_current/Delftia_tsuruhatensis/temporal/
validation_runs/delftia_current/Delftia_tsuruhatensis/panr2_inputs/
```

Representative cross-database outputs:

```text
cross_database/analysis/cross_database_feature_matrix.csv
cross_database/analysis/cross_database_top_associations.csv
cross_database/analysis/amr_integron_associations.csv
cross_database/analysis/amr_plasmid_associations.csv
cross_database/analysis/amr_virulence_associations.csv
cross_database/figures/html_files/cross_database_feature_network.html
cross_database/figures/html_files/global_feature_association_heatmap.html
cross_database/figures/html_files/integrated_feature_presence_heatmap.html
```

Temporal outputs:

```text
temporal/analysis/temporal_feature_trends.csv
temporal/analysis/temporal_burden_trends.csv
temporal/figures/html_files/temporal_top_feature_trends.html
```

## Interpretation Limits

Cross-database co-occurrence in PanR2 is sample/genome-level unless explicit same-contig and coordinate evidence is available. It does not by itself prove physical linkage, plasmid localization, horizontal transfer, or causality.

MobileElementFinder is currently opt-in because real validation exposed upstream parser fragility on some valid assemblies. ISfinder is supported through a user-supplied authorized FASTA because PanResistome does not redistribute or automatically download ISfinder.

