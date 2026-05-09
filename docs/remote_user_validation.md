# Remote User Validation

This document defines the release-blocking validation path for PanResistome. The goal is to prove that a user starting from a fresh clone can run the standard comprehensive workflow without manually debugging databases or hidden tool dependencies.

The broader release gates are tracked in [`release_reliability_checklist.md`](release_reliability_checklist.md).

## Scope

The recommended public path is:

```text
FetchM2 metadata/download
sequence QC
CheckM2
QUAST
FastANI/skani
Mash
combined QC
AMRFinderPlus
PanR2 comprehensive mode
PanR2 handoff export
HTML report
```

GTDB-Tk is intentionally disabled because it requires a large external database. MobileElementFinder remains opt-in. ISfinder requires a user-supplied authorized FASTA and is not part of the default public run.

## Fresh Clone Setup

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/PanResistome.git
cd PanResistome
```

Required user-installed programs:

```text
Nextflow
Java
Conda
Mamba, strongly recommended
Git
```

The pipeline environments install the bioinformatics tools used by the selected profile. The standard comprehensive command below should not require a user-provided CheckM2 database path, AMRFinderPlus database setup, ABRicate database setup, GTDB-Tk database, MobileElementFinder database, or ISfinder database.

## Quick Offline Check

Run this before a long validation:

```bash
nextflow run main.nf -profile test
```

This test uses tiny local fixtures and does not download genomes or large databases.

## Standard Comprehensive Validation Command

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

Use `--threads 4` on modest desktops. Increase only after a smaller run is stable.

The 2026-05-08 fresh-clone validation of this command on 45 `Delftia tsuruhatensis` assemblies completed in 4h 29m on a 24-thread laptop with 15 GiB RAM visible to the OS. The command used no `--checkm2_db` path and completed all 19 Nextflow processes.

## What The Command Must Do Automatically

The run must:

```text
Download/standardize metadata with FetchM2
Download assemblies with FetchM2 native downloader
Run CheckM2 and auto-download/cache its database when --checkm2_db is omitted
Run QUAST, ANI/skani, and Mash summaries
Run AMRFinderPlus and update/download its database when needed
Run PanR2 setup-db for ABRicate ncbi,vfdb,plasmidfinder
Verify required ABRicate databases before PanR2 analysis
Run PanResistome-native ABRicate, IntegronFinder, and MLST feature runners before PanR2 analysis
Export panr2_inputs/features/*.features.tsv
Export panr2_inputs/features/all_features.tsv
Validate feature contracts and metadata matching
Write the combined HTML dashboard
```

## Expected Main Outputs

```text
validation_runs/delftia_fresh/<organism>/report/index.html
validation_runs/delftia_fresh/<organism>/panr2_inputs/features/all_features.tsv
validation_runs/delftia_fresh/<organism>/panr2_inputs/manifest/database_setup_status.tsv
validation_runs/delftia_fresh/<organism>/panr2_inputs/manifest/schema_validation_summary.txt
validation_runs/delftia_fresh/<organism>/panr2_inputs/manifest/feature_completeness_audit.tsv
validation_runs/delftia_fresh/<organism>/panr2_inputs/manifest/module_status_summary.tsv
validation_runs/delftia_fresh/<organism>/panr2_inputs/manifest/software_versions.csv
validation_runs/delftia_fresh/<organism>/qc/qc_master_report.csv
validation_runs/delftia_fresh/pipeline_versions/
```

After the run completes, generate a compact release-evidence summary:

```bash
scripts/summarize_validation_run.py \
  --run-dir validation_runs/delftia_fresh \
  --out-dir validation_runs/delftia_fresh
```

This writes:

```text
validation_runs/delftia_fresh/validation_summary.csv
validation_runs/delftia_fresh/validation_summary.md
```

## Required Feature Tables

The standard comprehensive validation should normally produce:

```text
panr2_inputs/features/amr.features.tsv
panr2_inputs/features/amrfinderplus.features.tsv
panr2_inputs/features/vfdb.features.tsv
panr2_inputs/features/plasmidfinder.features.tsv
panr2_inputs/features/integronfinder.features.tsv
panr2_inputs/features/mlst.features.tsv
panr2_inputs/features/all_features.tsv
```

`mlst.features.tsv` may be header-only for organisms without a supported PubMLST scheme. That is an acceptable `WARNING_EMPTY` result when the raw MLST command completed and no real ST/allele calls were detected.

No enabled module may silently disappear. If a required table is missing, inspect:

```text
panr2_inputs/manifest/database_setup_status.tsv
panr2_inputs/manifest/feature_completeness_audit.tsv
panr2_inputs/manifest/module_status_summary.tsv
```

## Release-Passing Criteria

A fresh-clone validation is release-passing only if:

```text
The run starts from a clean clone and clean output directory.
No --checkm2_db path is supplied.
GTDB-Tk is disabled.
CheckM2 database setup succeeds or produces a clear actionable failure.
AMRFinderPlus database setup succeeds.
ABRicate ncbi/vfdb/plasmidfinder are present after setup.
database_setup_status.tsv has no FAIL rows where required_for_profile=true.
schema_validation_summary.txt reports zero unmatched, invalid, and duplicate feature rows.
feature_completeness_audit.tsv does not show missing required default feature tables.
report/index.html exists.
software_versions.csv exists.
The run command, runtime, machine specs, feature counts, and known limitations are documented.
```

## Optional Modules

These are intentionally outside the default public command:

```text
GTDB-Tk: large external reference database
MobileElementFinder: opt-in because upstream parser failures were observed on valid assemblies
ISfinder: requires authorized local FASTA; not auto-downloaded or redistributed
MOB-suite: optional plasmid reconstruction/typing
geNomad: requires a geNomad database path
Kleborate/Kaptive/ECTyper: organism-specific; table passthrough is preferred for difficult installs
DefenseFinder: available in PanR2 but not default until its environment is consistently stable
```

## Runtime Caveat

The standard comprehensive command is reliable but not yet fully optimized for wall time. PanResistome now owns the standard feature-runner stage through `--panr2_native_feature_runners true` and passes precomputed ABRicate, IntegronFinder, and MLST directories into PanR2. Large fragmented assemblies can still make IntegronFinder the longest substage; future performance work should split per-assembly execution into finer-grained Nextflow channels while keeping PanR2 focused on standardized analysis/reporting.

An experimental parallel backend is available for validation:

```bash
--panr2_native_feature_runner_mode parallel --threads 16
```

Parallel mode keeps the same output directories but runs one ABRicate database at a time with per-genome workers, then runs per-assembly IntegronFinder/MLST calls concurrently inside the native-runner stage. Keep `serial` as the stable fallback until the parallel Delftia and Klebsiella validations are documented.

On desktops with 16 GB RAM or less, do not let CheckM2 inherit a high global thread count. Use a separate cap:

```bash
--threads 16 --checkm2_threads 2
```

CheckM2 defaults to `min(--threads,4)`, but `--checkm2_threads 2` is safer for 100-genome validations on modest machines.

The 2026-05-08 native feature-runner validation is documented in:

```text
validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md
```

## Documentation After A Passing Run

After a passing run, add or update:

```text
validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md
CHANGELOG.md
README.md
```

Do not commit the full `validation_runs/` directory. It contains downloaded genomes, databases, raw tool outputs, and large intermediate files. Commit only reproducible commands, summary reports, and lightweight example-output documentation.
