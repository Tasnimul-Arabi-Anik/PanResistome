# Klebsiella pneumoniae Large-Mode Validation Results

Date: 2026-05-09

This validation tested the v0.4.0 large-dataset report safeguards on the committed 100-record `Klebsiella pneumoniae` input. GTDB-Tk was intentionally disabled. The run used the validated v0.3.x execution path with FetchM2, sequence QC, CheckM2, QUAST, ANI, Mash, AMRFinderPlus, PanResistome-native ABRicate/IntegronFinder/MLST runners, PanR2 comprehensive analysis, and PanR2 handoff export.

Large mode is report-facing only: complete feature TSVs remain complete, while report-facing matrices and cross-database summaries are capped for readability.

## Command

Initial command:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_large_mode \
  -profile conda,mamba,desktop_parallel,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true
```

The initial attempt started an automatic CheckM2 database download into the large-mode output directory. Because a validated CheckM2 database already existed from the prior Klebsiella validation, the run was stopped and resumed with the cached database to avoid re-downloading a large file:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_large_mode \
  -profile conda,mamba,desktop_parallel,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --checkm2_db /home/anik/genomics/Tools_Dev/PanResistome/validation_runs/klebsiella_pneumoniae_parallel/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
  -resume
```

## Nextflow Outcome

```text
Completed at: 2026-05-09 20:05:43
Duration: 4h 34m 35s
CPU hours: 26.6
Succeeded processes: 11
Cached processes: 9
Failed processes: 0
```

Completed stages:

```text
FetchM2 metadata/download
Sequence QC
CheckM2 QC
QUAST
ANI
Mash
Combined QC
AMRFinderPlus
PanR2 native feature runners
PanR2 comprehensive analysis
PanR2 handoff export
Result collection
```

## Input, Download, And QC Results

```text
Input records: 100
Downloaded FASTA files: 99
Failed downloads: 1
QC master PASS: 99
QC master FAIL: 0
```

The failed accession was handled as a recorded download failure and was not passed into downstream analysis.

## Feature Contract Results

PanR2 feature contract validation passed cleanly:

```text
feature_files_checked: 6
feature_rows: 12838
databases_seen: amr, amrfinderplus, integronfinder, mlst, plasmidfinder, vfdb
unmatched_feature_rows: 0
invalid_feature_rows: 0
duplicate_feature_rows: 0
```

Complete feature rows by database:

```text
amr: 1234
amrfinderplus: 1350
vfdb: 7821
plasmidfinder: 393
integronfinder: 418
mlst: 1622
all_features.tsv: 12838
```

Complete feature tables were not capped. The full merged table remained:

```text
panr2_inputs/features/all_features.tsv: 12838 data rows
```

## Large-Mode Report Controls

The large profile applied these settings:

```text
large_dataset: true
report_mode: compact
max_features_heatmap: 150
max_features_network: 150
max_metadata_columns: 20
top_n_features_per_database: 50
skip_heavy_interactive_plots: true
samples: 99
unique database-feature pairs: 673
```

Audit file:

```text
panr2_inputs/manifest/report_controls.tsv
```

HTML page:

```text
panr2_inputs/report/report_controls.html
```

## Capped Report-Facing Outputs

Large mode capped the report-facing matrix and network-style summaries:

```text
feature_matrices/all_features_presence_absence.tsv: 100 rows including header, 151 columns
cross_database/feature_cooccurrence.tsv: 11175 data rows, 150 unique features
cross_database/feature_proximity.tsv: 21921 data rows, 150 unique features
```

The 151 matrix columns represent the sample/accession column plus 150 selected feature columns.

Complete proximity evidence was preserved separately:

```text
cross_database/feature_proximity_all.tsv: 23714 data rows, 309 unique features
```

This confirms that large mode makes report-facing summaries manageable without discarding complete feature tables.

## Top-Feature Summary

The top-feature navigation table was generated:

```text
metadata_feature_analysis/top_features_by_database.tsv: 252 data rows
```

Rows by database:

```text
amr: 50
amrfinderplus: 50
integronfinder: 3
mlst: 50
plasmidfinder: 49
vfdb: 50
```

## Native Feature-Runner Audit

```text
ABRicate: PASS, 99 samples processed, 297 raw tables, 9448 feature rows
IntegronFinder: PASS, 99 samples processed, 99 raw tables, 418 feature rows
MLST: PASS, 99 samples processed, 763 native-runner feature rows
MobileElementFinder: SKIPPED, not enabled by default
```

Audit file:

```text
panr2_inputs/manifest/native_runner_merge_audit.tsv
```

## Report Outputs

The following report files were generated:

```text
report/index.html
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/report_controls.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/metadata_quality_and_bias.html
panr2_inputs/report/database_burden_by_metadata.html
panr2_inputs/report/cross_database_interpretation.html
panr2_inputs/report/database_setup_and_contract.html
panr2_inputs/report/bioproject_bias.html
panr2_inputs/report/amrfinder_abricate_concordance.html
```

## Resource Finding

CheckM2 was the dominant memory-sensitive stage. With `desktop_parallel,large`, CheckM2 used the profile cap:

```text
--checkm2_threads 2
```

On a 16 GB desktop-class machine, the DIAMOND phase still applied noticeable memory pressure but completed. For 100+ genome validations on 16 GB systems, use one of:

```bash
--checkm2_threads 1
```

or a cached/precomputed CheckM2 database and `-resume`. The large-mode report controls reduce report size and browser load; they do not reduce runtime for heavy annotation/QC tools.

## Conclusion

PASS. The existing 100-record Klebsiella validation successfully completed under the `large` profile. Large mode preserved complete feature exports, generated report-control audit files, capped report-facing matrices and co-occurrence/proximity summaries, and produced the expected HTML report set.

The next v0.4.0 validation target remains a 300-500 genome run with GTDB-Tk disabled, CheckM2 capped, parallel native feature runners, and large-mode report controls enabled.
