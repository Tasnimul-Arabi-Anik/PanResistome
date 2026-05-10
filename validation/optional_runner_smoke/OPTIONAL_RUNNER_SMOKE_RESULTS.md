# Optional Runner Smoke Validation

Date: 2026-05-10

Purpose: validate the optional-runner orchestration path without running the default heavy comparative-genomics stages.

This smoke test uses the two tiny local fixture assemblies under `tests/fixtures/local_samples/`. It intentionally keeps CheckM2, GTDB-Tk, QUAST, ANI, Mash, AMRFinderPlus, PanR2 comprehensive mode, and legacy ABRicate/PanR disabled so optional modules can be tested quickly.

## Command

```bash
nextflow run main.nf \
  -profile test \
  --local_samples tests/fixtures/local_samples \
  --outdir validation_runs/optional_runner_smoke_local_noconda_v5 \
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
  --run_mobsuite true \
  --run_genomad true \
  --run_kleborate true \
  --run_kaptive true \
  --run_ectyper true \
  --run_isfinder true \
  --isfinder_db_fasta validation/optional_runner_smoke/synthetic_isfinder.fasta \
  --threads 2 \
  --capture_versions false
```

## Result

Status: PASS

Completed processes:

```text
SEQUENCE_QC
COMBINED_QC
ISFINDER_BLAST
MOBSUITE_ANALYSIS
GENOMAD_PROPHAGE
ORGANISM_SPECIFIC_TYPING
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

The first conda-enabled attempt was intentionally not treated as the final validation because this session has no network access and no writable Conda package cache. It failed while trying to create the geNomad environment. The no-conda smoke test above validates the pipeline orchestration, skip/warning behavior, and PanR2 export path without external environment creation.

## Outputs Checked

PanR2 feature tables were created for the enabled optional modules:

| Feature table | Rows |
| --- | ---: |
| `ectyper.features.tsv` | 0 |
| `isfinder.features.tsv` | 0 |
| `kaptive.features.tsv` | 0 |
| `kleborate.features.tsv` | 0 |
| `mobsuite.features.tsv` | 0 |
| `prophage.features.tsv` | 0 |

Schema validation:

```text
feature_files_checked=6
feature_rows=0
metadata_accessions=2
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature completeness audit showed `WARNING_EMPTY` for each enabled optional module with raw output but no biological calls:

```text
ectyper
isfinder
kaptive
kleborate
mobsuite
prophage
```

ISfinder-compatible BLAST status:

```text
status=PASS
samples_input=2
samples_processed=2
samples_failed=0
raw_tables_created=2
feature_rows_created=0
```

The zero ISfinder feature count is expected because the synthetic FASTA is only a legal/local smoke-test database, not a biological ISfinder database.

## Issues Found and Fixed

1. Optional runner processes requested fixed 8 CPUs, which blocked fast 1-CPU test-profile validation. The optional process CPU requests now scale with `--threads`, and test-profile overrides are provided.
2. The legacy non-comprehensive branch always launched ABRicate/PanR after optional modules. `--run_abricate false` now allows optional-runner/table-input validation without default ABRicate.
3. Optional runner outputs with raw files but zero biological rows were incorrectly audited as `FAIL_MISSING_FEATURE_TABLE`. PanR2 export now writes header-only feature tables and records `WARNING_EMPTY`.
4. Kleborate/Kaptive/ECTyper placeholder tables now use valid tab-separated headers.

## Interpretation

This validation does not prove that MOB-suite, geNomad, Kleborate, Kaptive, or ECTyper are biologically useful in this environment because their tools/databases were not available. It does prove that:

```text
optional runner orchestration works,
ISfinder authorized-FASTA execution path works,
missing optional tools/databases produce auditable empty outputs,
PanR2 feature-contract export handles optional empty outputs cleanly,
and optional-only validation can run without default ABRicate/PanR.
```

The next validation step should use real installed optional tools/databases on a 2-10 genome subset.
