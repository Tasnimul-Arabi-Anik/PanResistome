# Native Feature-Runner Validation Results

Validation date: 2026-05-08

Validation target: `Delftia tsuruhatensis` current NCBI Assembly input

Input file:

```text
validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv
```

Run directory:

```text
validation_runs/delftia_native_runner
```

## Purpose

This validation checks the v0.3.0 native feature-runner architecture:

```text
PanResistome runs standard feature-detection tools first.
PanResistome exports raw outputs and PanR2-compatible inputs.
PanR2 performs standardized analysis and reporting from precomputed directories.
```

The goal was to verify that ABRicate, IntegronFinder, and MLST can be executed under PanResistome ownership before PanR2 reporting, rather than being run as hidden heavy work inside PanR2.

## Command

```bash
NXF_DISABLE_CHECK_LATEST=true nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_native_runner \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --panr2_native_feature_runners true \
  --checkm2_db results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
  --threads 4 \
  --fetchm2_download_workers 2
```

GTDB-Tk was disabled. A cached CheckM2 database path was supplied because automatic CheckM2 database download had already been validated in the fresh-clone validation documented in `FRESH_CLONE_VALIDATION_RESULTS.md`.

## Nextflow Result

The full run completed successfully.

```text
Processes: 20 succeeded / 20 total
Duration: 3h 47m 20s
CPU hours: 20.7
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
PANR2_FEATURE_RUNNERS
PANR2_COMPREHENSIVE
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

## Sample and QC Summary

From `validation_runs/delftia_native_runner/validation_summary.md`:

```text
metadata_rows: 45
downloaded_fastas: 45
qc_master_status: 45 PASS
schema_feature_rows: 201
schema_unmatched_feature_rows: 0
schema_invalid_feature_rows: 0
schema_duplicate_feature_rows: 0
database_setup_required_failures: 0
dashboard_exists: True
```

Database/tool setup status included PASS for FetchM2 metadata, sequence FASTA inputs, CheckM2, QUAST, ANI, Mash, PanR2, ABRicate, ABRicate `ncbi`, ABRicate `vfdb`, ABRicate `plasmidfinder`, IntegronFinder, MLST, and AMRFinderPlus. GTDB-Tk, MobileElementFinder, DefenseFinder, ISfinder, MOB-suite, geNomad, and Kaptive were skipped because they were not requested for this run.

## Native Feature-Runner Evidence

`PANR2_FEATURE_RUNNERS` completed successfully.

Native runner status:

```text
abricate: PASS
integronfinder: PASS
mlst: PASS
mobileelementfinder: SKIPPED
```

Interpretation:

- ABRicate ran under PanResistome for `ncbi`, `vfdb`, and `plasmidfinder`.
- IntegronFinder ran under PanResistome for all 45 genomes.
- MLST ran under PanResistome for all 45 genomes.
- MobileElementFinder remained opt-in and was not enabled.
- PanR2 then consumed precomputed feature-runner directories during reporting.

## Feature Contract Result

Final standardized feature tables:

```text
panr2_inputs/features/amr.features.tsv: 19 rows
panr2_inputs/features/amrfinderplus.features.tsv: 24 rows
panr2_inputs/features/vfdb.features.tsv: 40 rows
panr2_inputs/features/plasmidfinder.features.tsv: 2 rows
panr2_inputs/features/integronfinder.features.tsv: 116 rows
panr2_inputs/features/mlst.features.tsv: 0 rows
panr2_inputs/features/all_features.tsv: 201 rows
```

`mlst.features.tsv` is intentionally header-only in this validation. The raw MLST command completed, but `Delftia tsuruhatensis` did not produce a supported PubMLST scheme/ST call in the native `mlst` output. Placeholder calls such as `-` are retained as run evidence in raw outputs but are not converted into biological features.

Schema validation:

```text
feature_files_checked=6
feature_rows=201
databases_seen=amr,amrfinderplus,integronfinder,plasmidfinder,vfdb
samples_seen=43
metadata_accessions=45
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature completeness audit:

```text
amr: PASS
amrfinderplus: PASS
vfdb: PASS
plasmidfinder: PASS
integronfinder: PASS
mlst: WARNING_EMPTY
```

The `WARNING_EMPTY` MLST status is expected for this organism/input and is preferable to creating misleading `ST_-` placeholder features.

## v0.3.0 Interpretation Outputs

The run produced the expected v0.3.0 analysis and handoff files, including:

```text
panr2_inputs/metadata_feature_analysis/feature_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/database_burden_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/category_burden_by_sample.tsv
panr2_inputs/metadata_feature_analysis/category_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/top_findings.tsv
panr2_inputs/metadata_feature_analysis/top_findings.md
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/database_cooccurrence_summary.tsv
panr2_inputs/cross_database/amr_mge_same_contig.tsv
panr2_inputs/cross_database/amr_plasmid_same_contig.tsv
panr2_inputs/cross_database/amr_integron_same_contig.tsv
panr2_inputs/cross_database/feature_proximity.tsv
panr2_inputs/cross_database/amrfinder_abricate_concordance.tsv
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/metadata_quality_and_bias.html
panr2_inputs/report/database_burden_by_metadata.html
panr2_inputs/report/cross_database_interpretation.html
panr2_inputs/report/database_setup_and_contract.html
```

The top-level PanR2 dashboard was also generated:

```text
report/index.html
```

## Issue Found and Fixed

The validation found an important edge case:

```text
Native MLST raw output can be headerless and can contain one placeholder row per genome when no PubMLST scheme/ST is found.
```

Before the fix, the raw MLST output was detected but did not become an explicit standardized feature table. The contract exporter now:

- Parses headerless native `mlst` output.
- Converts real ST and allele calls into PanR2 `mlst` features.
- Suppresses placeholder calls such as `-`, `ST_-`, and `-:ST-`.
- Writes an explicit header-only `mlst.features.tsv` when MLST ran successfully but produced no biological ST/allele features.

The native feature-runner status code was also updated so placeholder MLST rows are not counted as biological feature rows.

## Performance Notes

The native feature-runner path is functionally valid, but this validation also identified the next scalability target:

- ABRicate and IntegronFinder were the slowest standard feature-runner steps.
- The current native runner still executes these tools serially inside one Nextflow process.
- The next v0.3.0 engineering step should split ABRicate by database/sample and IntegronFinder/MLST by sample into parallel Nextflow processes, while preserving the current native runner as a stable fallback.

## Conclusion

This validation passes the v0.3.0 groundwork target:

```text
PanResistome can run standard feature runners before PanR2,
export standardized feature tables,
generate clean feature-contract validation,
and produce metadata/cross-database interpretation outputs.
```

The next release work should focus on runner parallelization and a second-organism validation using a richer AMR/plasmid/MLST organism such as `Klebsiella pneumoniae`.
