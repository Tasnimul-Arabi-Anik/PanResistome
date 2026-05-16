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

GTDB-Tk is intentionally disabled because it requires a large external database. MobileElementFinder remains opt-in and nonblocking when explicitly enabled. ISfinder requires a user-supplied authorized FASTA and is not part of the default public run. DefenseFinder remains table-input/experimental and is excluded from the recommended comprehensive route.

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

Current Acinetobacter Docker/GHCR target:

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 10 \
  --outdir validation_runs/acinetobacter_pittii_10_docker \
  -profile docker \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder false \
  --panr2_run_defensefinder false \
  --threads 8 \
  --checkm2_threads 4 \
  --fetchm2_download_workers 2
```

Equivalent Conda/Mamba target:

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 10 \
  --outdir validation_runs/acinetobacter_pittii_10_conda \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --panr2_run_mobileelementfinder false \
  --panr2_run_defensefinder false \
  --threads 8 \
  --checkm2_threads 4 \
  --fetchm2_download_workers 2
```

Historical Delftia fresh-clone command:

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

On 2026-05-16, a local fixture run reproduced a CheckM2 model-loading regression from the stale environment route:

```text
Saved models could not be loaded ... specific_model_COMP.keras
```

The intended CheckM2 package route is now `checkm2=1.1.0=pyh7e72e81_1` with CPU TensorFlow 2.17 and Python 3.12 in `envs/checkm2.yaml`. On 2026-05-16, a fresh temporary Conda environment and the GHCR Docker image both loaded the CheckM2 model and completed `checkm2 predict` against a local full-genome fixture, producing one real `quality_report.tsv` genome row with no Keras/model-loading error. The Docker-profile Nextflow QC fixture also passed with `checkm2_model_load=PASS`.

The current Docker/GHCR comprehensive route passed on 2026-05-16 with 5 `Acinetobacter pittii` records, no `--checkm2_db` argument, GTDB-Tk/DefenseFinder/MobileElementFinder/ISfinder disabled, and CheckM2/QUAST/ANI/Mash/AMRFinderPlus/geNomAD/native ABRicate/IntegronFinder/MLST/PanR2 comprehensive enabled. The run completed in 2h 21m 30s on a constrained desktop configuration (`--threads 4 --checkm2_threads 1`), produced 5/5 CheckM2 rows, 5/5 combined QC PASS calls, 630 PanR2 feature rows, and zero unmatched, invalid, or duplicate feature rows. See `validation/deployment/ACINETOBACTER_CHECKM2_VALIDATION_STATUS.md`.

## What The Command Must Do Automatically

The run must:

```text
Download/standardize metadata with FetchM2
Download assemblies with FetchM2 native downloader
Run CheckM2 and auto-download/cache its database when --checkm2_db is omitted
Run QUAST, ANI/skani, and Mash summaries
Run AMRFinderPlus and update/download its database when needed
Run PanR2 setup-db for ABRicate ncbi,vfdb,plasmidfinder. By default, `--panr2_update_abricate_db true` force-refreshes requested ABRicate databases with `abricate-get_db --force` when available; set it false only for offline/cached reruns.
Verify required ABRicate databases before PanR2 analysis
Run PanResistome-native ABRicate, IntegronFinder, and MLST feature runners before PanR2 analysis
Export panr2_inputs/features/*.features.tsv
Export panr2_inputs/features/all_features.tsv
Validate feature contracts and metadata matching
Write the combined HTML dashboard
```

## Expected Main Outputs

```text
validation_runs/acinetobacter_pittii_10_docker/<organism>/report/index.html
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/features/all_features.tsv
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/database_setup_status.tsv
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/abricate_database_setup_status.tsv
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/schema_validation_summary.txt
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/feature_completeness_audit.tsv
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/module_status_summary.tsv
validation_runs/acinetobacter_pittii_10_docker/<organism>/panr2_inputs/manifest/software_versions.csv
validation_runs/acinetobacter_pittii_10_docker/<organism>/qc/qc_master_report.csv
validation_runs/acinetobacter_pittii_10_docker/pipeline_versions/
```

After the run completes, generate a compact release-evidence summary:

```bash
scripts/summarize_validation_run.py \
  --run-dir validation_runs/acinetobacter_pittii_10_docker \
  --out-dir validation_runs/acinetobacter_pittii_10_docker
```

This writes:

```text
validation_runs/acinetobacter_pittii_10_docker/validation_summary.csv
validation_runs/acinetobacter_pittii_10_docker/validation_summary.md
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
CheckM2 database setup succeeds, `checkm2/quality_report.tsv` has genome rows, and no `specific_model_COMP.keras`/Keras model-load error appears.
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
MobileElementFinder: opt-in and nonblocking by default; when enabled, failures must produce header-only PanR2 output, a module status row, a native runner merge-audit row, and a visible warning
ISfinder: requires authorized local FASTA; not auto-downloaded or redistributed
MOB-suite: optional plasmid reconstruction/typing
geNomad: requires a geNomad database path
Kleborate/Kaptive/ECTyper: organism-specific; table passthrough is preferred for difficult installs
DefenseFinder: table-input/experimental; leave `--panr2_run_defensefinder false` for broad remote runs
```

## MobileElementFinder Validation Track

After the stable CheckM2-on comprehensive run passes, run a separate Acinetobacter validation with MobileElementFinder enabled:

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 10 \
  --outdir validation_runs/acinetobacter_pittii_10_mobileelementfinder \
  -profile docker \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --panr2_run_mobileelementfinder true \
  --panr2_mobileelementfinder_allow_failure true \
  --panr2_run_defensefinder false \
  --threads 8 \
  --checkm2_threads 4 \
  --fetchm2_download_workers 2
```

Acceptable MobileElementFinder outcomes are either `PASS` with real MGE rows, or `WARNING_FAILED` with clean header-only output and no workflow abort.

## Runtime Caveat

The standard comprehensive command is reliable but not yet fully optimized for wall time. PanResistome now owns the standard feature-runner stage through `--panr2_native_feature_runners true` and passes precomputed ABRicate, IntegronFinder, and MLST directories into PanR2. Large fragmented assemblies can still make IntegronFinder the longest substage; future performance work should split per-assembly execution into finer-grained Nextflow channels while keeping PanR2 focused on standardized analysis/reporting.

The parallel backend is validated for the committed Delftia and Klebsiella validations:

```bash
--panr2_native_feature_runner_mode parallel --threads 16
```

Parallel mode keeps the same output directories but runs one ABRicate database at a time with per-genome workers, then runs per-assembly IntegronFinder/MLST calls concurrently inside the native-runner stage. `serial` remains the conservative fallback.

The validated desktop-scale profile is:

```bash
-profile conda,mamba,desktop_parallel
```

It sets `--threads 16`, `--checkm2_threads 2`, `--fetchm2_download_workers 2`, and `--panr2_native_feature_runner_mode parallel`.

On desktops with 16 GB RAM or less, do not let CheckM2 inherit a high global thread count. Use a separate cap:

```bash
--threads 16 --checkm2_threads 2
```

CheckM2 defaults to `min(--threads,4)`, but `--checkm2_threads 2` is safer for 100-genome validations on modest machines.

The native feature-runner validations are documented in:

```text
validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md
validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md
validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md
```

## Documentation After A Passing Run

After a passing run, add or update:

```text
validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md
CHANGELOG.md
README.md
```

Do not commit the full `validation_runs/` directory. It contains downloaded genomes, databases, raw tool outputs, and large intermediate files. Commit only reproducible commands, summary reports, and lightweight example-output documentation.
