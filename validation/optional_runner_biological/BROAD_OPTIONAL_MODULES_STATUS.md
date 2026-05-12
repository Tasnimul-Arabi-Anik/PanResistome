# Broad Optional Module Validation Status

Date: 2026-05-12

Scope:

```text
1. geNomad/prophage
2. DefenseFinder
3. MobileElementFinder
```

Purpose: validate whether these broad opt-in modules can feed the same PanR2 analysis layer used by AMR/VFDB, and clearly separate table-analysis readiness from real external-runner readiness.

## Summary

Status: PARTIAL PASS

PanR2 analysis path: PASS

Runner biological validation on this machine: PARTIAL PASS

The active shell did not have the external runner commands available:

```text
genomad: not found
defense-finder: not found
mefinder: not found
```

However, cached Nextflow Conda environments were available from prior workflow runs. Those cached environments allowed a small real MobileElementFinder validation on 5 `Klebsiella pneumoniae` genomes. geNomad and DefenseFinder still did not pass real biological runner validation in this round.

This means the current validation can honestly prove that standardized outputs from all three modules are analyzed correctly by PanR2, and can now additionally claim that MobileElementFinder biologically ran on a small real Klebsiella subset and fed PanR2 clean standardized outputs.

## PanR2 Analysis Validation

Command:

```bash
python scripts/validate_optional_feature_analysis.py \
  --outdir validation/optional_feature_analysis \
  --force
```

Result: PASS

Schema validation:

```text
feature_files_checked=13
feature_rows=23
databases_seen=amr,defensefinder,ectyper,isfinder,kaptive,kleborate,mobileelementfinder,mobsuite,plasmidfinder,prophage,sccmecfinder,serotypefinder,vfdb
samples_seen=6
metadata_accessions=6
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Relevant broad optional feature rows:

| Database | Feature rows | Status |
| --- | ---: | --- |
| `mobileelementfinder` | 1 | PASS |
| `prophage` / geNomad-style table | 1 | PASS |
| `defensefinder` | 1 | PASS |

Relevant PanR2 analysis outputs were generated:

```text
panr2_inputs/features/mobileelementfinder.features.tsv
panr2_inputs/features/prophage.features.tsv
panr2_inputs/features/defensefinder.features.tsv
panr2_inputs/features/all_features.tsv
panr2_inputs/feature_matrices/all_features_presence_absence.tsv
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/database_cooccurrence_summary.tsv
panr2_inputs/cross_database/feature_proximity.tsv
panr2_inputs/metadata_feature_analysis/top_features_by_database.tsv
panr2_inputs/metadata_feature_analysis/top_findings.tsv
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/cross_database_interpretation.html
```

Analysis output sizes:

```text
all_features.tsv: 23 rows
all_features_presence_absence.tsv: 6 samples
feature_cooccurrence.tsv: 190 rows
database_cooccurrence_summary.tsv: 78 rows
feature_proximity.tsv: 8 rows
top_features_by_database.tsv: 20 rows
top_findings.tsv: 50 rows
```

Interpretation: if geNomad, DefenseFinder, or MobileElementFinder outputs are present as tables, PanR2 can analyze them like AMR/VFDB: feature contract export, merged feature table, feature matrix, metadata summaries, co-occurrence/proximity, top findings, and HTML handoff pages.

## Runner Status By Module

### geNomad / Prophage

Current pipeline support:

```text
--run_genomad true
--genomad_db /path/to/genomad_db
--genomad_db_dir <outdir>/databases/genomad
--genomad_auto_download_db true
```

What is validated:

```text
geNomad-style table input -> prophage.features.tsv -> PanR2 analysis: PASS
database setup helper unit test with fake genomad executable: PASS
optional-runner smoke path with missing DB/tool produces auditable empty outputs: PASS
```

What is not yet validated:

```text
real geNomad biological run with a real geNomad database
fresh remote-user database download time/disk behavior
prophage feature quality on bacterial genomes
```

Additional local probe:

```text
genomad executable in active shell: false
auto-download requested: false
database status: FAIL, as expected without a DB
status file: validation_runs/genomad_missing_db_probe/genomad_database_setup_status.tsv
```

Interpretation: the geNomad setup helper records an auditable missing-DB state correctly. A true remote-user validation still requires a fresh Conda solve plus `--genomad_auto_download_db true`, or a supplied `--genomad_db`, and should be run first on 2-10 genomes because the geNomad database is large.

Auto-download attempt:

```text
command reached GENOMAD_PROPHAGE: yes
environment creation: did not finish in practical validation window
time before stop: approximately 17 minutes
database download reached: no
biological geNomad outputs generated: no
```

Detailed result: `validation/optional_runner_biological/GENOMAD_AUTO_DOWNLOAD_ATTEMPT_RESULTS.md`

Interpretation: geNomad remains opt-in and should not yet be described as hassle-free automated setup. The next validation should either prebuild/cache the geNomad environment, or validate a container/Apptainer route before attempting the large database download again.

Host-environment bypass smoke:

```text
-profile genomad_host: PASS
GENOMAD_PROPHAGE ran without creating envs/genomad.yaml: yes
runtime with missing host genomad/database: 0.10s
auditable missing-database status: PASS
biological geNomad outputs generated: no, expected for this negative smoke
```

Detailed result: `validation/optional_runner_biological/GENOMAD_HOST_PROFILE_SMOKE_RESULTS.md`

Interpretation: users with a prebuilt host/module/container `genomad` can now bypass only the geNomad Conda solve using `-profile genomad_host` or `--genomad_use_host_env true`. A positive biological validation still needs a working `genomad` executable and populated database.

Recommended next command when a real geNomad DB is available:

```bash
nextflow run main.nf \
  -profile conda,mamba,desktop_parallel,large \
  --local_samples validation_runs/optional_real_100_input \
  --outdir validation_runs/genomad_small_real_validation \
  --sequence_qc_engine python \
  --qc_filter true \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive false \
  --run_abricate false \
  --export_panr2_inputs true \
  --run_genomad true \
  --genomad_db /path/to/genomad_db \
  --threads 8 \
  --capture_versions false
```

Recommended status: keep geNomad opt-in until one real small biological validation passes.

### DefenseFinder

Current pipeline support:

```text
--defensefinder_dir /path/to/defensefinder_tables
--panr2_run_defensefinder true
```

What is validated:

```text
DefenseFinder table input -> defensefinder.features.tsv -> PanR2 analysis: PASS
feature completeness audit coverage: PASS
metadata/cross-database/top-finding outputs from standardized rows: PASS
```

What is not yet validated:

```text
DefenseFinder runner mode as a stable PanResistome-owned process
fresh install/database setup for DefenseFinder
real biological run on bacterial genomes in this environment
```

Additional local probe:

```text
cached defense-finder executable found: true
cached CLI import status: FAIL
failure: ModuleNotFoundError: No module named 'macsypy'
```

Interpretation: this confirms that DefenseFinder should remain a table-input path for now. It should not be advertised as a PanResistome-owned runner until a clean environment and database setup path can be validated.

Recommended near-term path: keep DefenseFinder as stable table input. Add a PanResistome-owned runner only after its database/environment setup can be made reproducible and audited.

Recommended status: stable table input, runner not default.

### MobileElementFinder

Current pipeline support:

```text
--panr2_run_mobileelementfinder true
--panr2_mobileelementfinder_allow_failure true
```

What is validated:

```text
MobileElementFinder table input -> mobileelementfinder.features.tsv -> PanR2 analysis: PASS
auditable nonfatal failure path in native feature-runner layer: implemented
optional module matrix documents upstream parser fragility
```

What is now biologically validated:

```text
direct runner on 5 Klebsiella genomes: PASS
raw MobileElementFinder CSV files: 5/5
PanR2 table export: PASS
Nextflow integration with amr_basic + MobileElementFinder: PASS
Nextflow processes succeeded: 6/6
final feature rows: 753
databases_seen: amr,mobileelementfinder
unmatched_feature_rows: 0
invalid_feature_rows: 0
duplicate_feature_rows: 0
```

Detailed result: `validation/optional_runner_biological/KLEBSIELLA_5_MOBILEELEMENTFINDER_RESULTS.md`

What is still not validated:

```text
upstream parser stability across valid assemblies
large-scale runtime behavior
```

Recommended next command when a working MobileElementFinder installation is available:

```bash
nextflow run main.nf \
  -profile conda,mamba,desktop_parallel \
  --local_samples validation_runs/optional_real_100_input \
  --outdir validation_runs/mobileelementfinder_small_real_validation \
  --sequence_qc_engine python \
  --qc_filter true \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive true \
  --panr2_native_feature_runners true \
  --panr2_run_mobileelementfinder true \
  --panr2_mobileelementfinder_allow_failure true \
  --threads 4 \
  --capture_versions false
```

Recommended status: biologically validated opt-in runner on a small Klebsiella subset; keep opt-in and nonfatal by default until multi-species and larger-scale behavior are validated.

## Practical Conclusion

These three modules should not be added to default comprehensive mode yet.

They are ready for PanR2 analysis when users supply valid outputs, but only geNomad currently has a first-class PanResistome runner/database-helper path. DefenseFinder should remain table-input until a reproducible runner is added. MobileElementFinder should remain opt-in because previous real-data validation exposed upstream parser fragility.

Current honest support level:

| Module | PanR2 analysis | Runner status | Default? |
| --- | --- | --- | --- |
| geNomad/prophage | PASS | Experimental opt-in, needs real DB validation | No |
| DefenseFinder | PASS | Table-input recommended | No |
| MobileElementFinder | PASS | Biologically validated small opt-in runner, nonfatal failure path | No |
