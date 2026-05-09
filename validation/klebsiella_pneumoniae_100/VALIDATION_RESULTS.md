# Klebsiella pneumoniae Parallel Validation Results

Date: 2026-05-09

This validation tested the v0.3.0 native feature-runner path on the committed 100-record `Klebsiella pneumoniae` NCBI Assembly input. GTDB-Tk was intentionally disabled. The run used FetchM2 metadata and sequence download, sequence QC, CheckM2, QUAST, ANI, Mash, AMRFinderPlus, PanResistome-native ABRicate/IntegronFinder/MLST runners, PanR2 comprehensive analysis, and PanR2 handoff export.

## Command

```bash
NXF_DISABLE_CHECK_LATEST=true nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_pneumoniae_parallel \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 16 \
  --checkm2_threads 2 \
  --fetchm2_download_workers 2 \
  -resume
```

## Resource Finding

An initial attempt allowed CheckM2 to inherit `--threads 16`. The CheckM2 DIAMOND stage exhausted RAM on a 16 GB desktop-class machine and made the PC unresponsive. The pipeline was updated to add `--checkm2_threads`, defaulting to `min(--threads, 4)`, and the validation was resumed with `--checkm2_threads 2`.

With the CheckM2 cap, the CheckM2 database was reused from the successful automatic download under:

```text
validation_runs/klebsiella_pneumoniae_parallel/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
```

The low-memory CheckM2 run completed successfully. This confirms that high-throughput feature runners can use more workers while CheckM2 remains capped for RAM safety.

## Nextflow Outcome

The resumed workflow completed successfully:

```text
succeededCount: 11
failedCount: 0
cachedCount: 9
peakCpus: 16
peakMemory: 8 GB declared
```

Completed stages included:

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

## Input and Download Results

```text
Input records: 100
Metadata rows after FetchM2/download: 99
Downloaded FASTA files: 99
Failed accession: GCF_055382775.1
QC master status: 99 PASS
```

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

Feature rows by database:

```text
amr: 1234
amrfinderplus: 1350
vfdb: 7821
plasmidfinder: 393
integronfinder: 418
mlst: 1622
all_features.tsv: 12838
```

## Native Feature-Runner Results

The PanResistome-native feature-runner layer completed successfully:

```text
ABRicate: PASS, 99 samples processed, 297 raw tables
IntegronFinder: PASS, 99 samples processed, 99 raw tables
MLST: PASS, 99 samples processed
MobileElementFinder: SKIPPED, not enabled by default
```

## Database/Tool Setup Audit

Required enabled setup checks passed:

```text
fetchm2_metadata: PASS
sequence_fasta_inputs: PASS
checkm2_database_and_qc: PASS
quast: PASS
ani: PASS
mash: PASS
panr2: PASS
abricate_tool: PASS
abricate_db:ncbi: PASS
abricate_db:vfdb: PASS
abricate_db:plasmidfinder: PASS
integronfinder: PASS
mlst: PASS
amrfinderplus_database_and_runner: PASS
```

Intentionally skipped optional modules:

```text
GTDB-Tk
MobileElementFinder
DefenseFinder
ISfinder authorized FASTA
MOB-suite
geNomad
Kaptive
```

## Report Outputs

Main outputs generated:

```text
validation_runs/klebsiella_pneumoniae_parallel/validation_summary.md
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/report/index.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/index.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/panr2_handoff_index.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/top_findings.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/metadata_quality_and_bias.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/database_burden_by_metadata.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/cross_database_interpretation.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/database_setup_and_contract.html
```

## Interpretation

This validation passes the v0.3.0 second-organism goal for the stable default comprehensive path:

- PanResistome ran the heavy tool layer.
- PanR2 consumed standardized feature outputs.
- The parallel native feature-runner backend handled a biologically rich organism.
- FetchM2 metadata, QC, AMR, virulence, plasmid, integron, MLST, AMRFinderPlus, metadata association, temporal, and cross-database report outputs were generated.
- Feature contract validation found zero unmatched, invalid, or duplicate rows.

The important operational lesson is that `--threads` should not be used blindly for every process. On modest desktops, keep `--threads` high enough for safe parallel feature runners, but cap CheckM2 separately with `--checkm2_threads 2` or `4`.
