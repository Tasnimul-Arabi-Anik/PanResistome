# Fresh-Clone Validation Results

Validation date: 2026-05-08

Repository state tested:

```text
commit: bf241634d1e7bacf9539a5e891ab0329a863abb7
clone directory: /tmp/panresistome_remote_user_bf24163
```

Machine:

```text
OS: Ubuntu 24.04 / Linux 6.11.0-29-generic
CPU threads visible: 24
Memory visible: 15 GiB
Swap visible: 4.0 GiB
```

## Command

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_fresh \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --threads 4 \
  --fetchm2_download_workers 2
```

No `--checkm2_db` path was supplied.

## Nextflow Result

```text
Completed: 2026-05-08 14:18:53 Asia/Dhaka
Duration: 4h 29m 21s
CPU hours: 30.9
Succeeded processes: 19
Failed processes: 0
Peak running processes: 6
Peak CPUs: 24
Peak memory: 11 GB
```

Completed processes:

```text
FETCHM_ENV_VERSIONS
ABRICATE_ENV_VERSIONS
PANR2_COMPREHENSIVE_ENV_VERSIONS
AMRFINDERPLUS_ENV_VERSIONS
CHECKM2_ENV_VERSIONS
ANI_ENV_VERSIONS
QUAST_ENV_VERSIONS
MASH_ENV_VERSIONS
FETCHM
SEQUENCE_QC
CHECKM2_QC
QUAST_QC
ANI_ANALYSIS
MASH_PRESCREEN
COMBINED_QC
AMRFINDERPLUS_ANALYSIS
PANR2_COMPREHENSIVE
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

## Database And Tool Setup Evidence

The run automatically downloaded the CheckM2 database:

```text
validation_runs/delftia_fresh/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
```

CheckM2 reported successful database download and completed quality prediction.

AMRFinderPlus ran `amrfinder_update` automatically and downloaded:

```text
AMRFinderPlus database version: 2026-03-24.1
AMRFinderPlus version: 4.2.7
```

ABRicate setup/status passed for the required comprehensive databases:

```text
ncbi
vfdb
plasmidfinder
```

Required setup failures:

```text
0
```

GTDB-Tk, MobileElementFinder, DefenseFinder, ISfinder, MOB-suite, geNomad, and Kaptive were skipped because they are not part of the standard public command unless explicitly requested with their required external databases or authorized inputs.

## Data Summary

```text
Input metadata rows: 45
Downloaded FASTA files: 45
QC master PASS: 45
QC master WARNING: 0
QC master FAIL: 0
```

Main outputs:

```text
validation_runs/delftia_fresh/Delftia_tsuruhatensis/report/index.html
validation_runs/delftia_fresh/Delftia_tsuruhatensis/panr2_inputs/features/all_features.tsv
validation_runs/delftia_fresh/Delftia_tsuruhatensis/panr2_inputs/manifest/database_setup_status.tsv
validation_runs/delftia_fresh/Delftia_tsuruhatensis/panr2_inputs/manifest/schema_validation_summary.txt
validation_runs/delftia_fresh/Delftia_tsuruhatensis/panr2_inputs/manifest/feature_completeness_audit.tsv
validation_runs/delftia_fresh/Delftia_tsuruhatensis/panr2_inputs/manifest/module_status_summary.tsv
```

## PanR2 Feature Contract

Schema validation summary:

```text
feature_files_checked=6
feature_rows=291
databases_seen=amr,amrfinderplus,integronfinder,mlst,plasmidfinder,vfdb
metadata_accessions=45
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature rows by database:

```text
amr: 19
amrfinderplus: 24
vfdb: 40
plasmidfinder: 2
integronfinder: 116
mlst: 90
all_features: 291
```

Cross-database and metadata-analysis files were generated, including:

```text
cross_database/feature_cooccurrence.tsv
cross_database/database_cooccurrence_summary.tsv
cross_database/amr_mge_context.tsv
cross_database/amr_plasmid_context.tsv
cross_database/amrfinder_abricate_concordance.tsv
metadata_feature_analysis/database_burden_by_sample.tsv
metadata_feature_analysis/fetchm2_metadata_audit.tsv
metadata_feature_analysis/metadata_column_eligibility.tsv
metadata_feature_analysis/top_findings.md
```

## Runtime Notes

The run is release-passing for the standard public comprehensive command.

PanR2 comprehensive mode currently runs IntegronFinder inside the PanR2 process. This completed for all 45 genomes, but fragmented assemblies created many per-contig files and made this stage the longest PanR2 substage. Future optimization should move per-sample annotation runners into native Nextflow processes and keep PanR2 focused on standardized analysis/reporting.

The generated MLST output for this organism contains unknown sequence-type calls from `mlst`; a follow-up exporter patch suppresses placeholder features such as `ST_-` and `-:ST-` from future PanR2 contract exports.
