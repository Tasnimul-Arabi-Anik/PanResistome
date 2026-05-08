# Native Feature-Runner Parallel Comparison

Date: 2026-05-08

This validation compares the previously documented 45-genome `Delftia tsuruhatensis`
native feature-runner path with the experimental parallel backend.

## Purpose

The v0.3.0 native feature-runner work moves standard feature execution under
PanResistome ownership before PanR2 reporting. The serial path was already
validated. This run tested whether `--panr2_native_feature_runner_mode parallel`
can use a 16-core local workstation while preserving the same standardized
PanR2 feature contract outputs.

## Parallel Strategy

The validated parallel backend uses a conservative database-isolated model:

- ABRicate runs one database at a time.
- Within each ABRicate database, up to `--threads` genomes are processed in
  parallel and then merged into the same database-level result/summary tables
  consumed by PanR2.
- IntegronFinder runs per assembly with one CPU per assembly and up to
  `--threads` concurrent assemblies.
- MLST runs per assembly with up to `--threads` concurrent assemblies.
- MobileElementFinder remains opt-in and is not enabled by force-rerun mode.

This preserves the existing output structure while avoiding unnecessary
competition between multiple ABRicate databases.

## Command

```bash
NXF_DISABLE_CHECK_LATEST=true nextflow run main.nf \
  --local_samples validation_runs/delftia_native_runner \
  --outdir validation_runs/delftia_parallel_feature_only_v2 \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_force_tool_run true \
  --threads 16 \
  --fetchm2_download_workers 2
```

This was a controlled feature-runner validation using the existing Delftia local
sample directory. Heavy QC and database-download stages were intentionally
disabled because they were already validated in the fresh-user run.
AMRFinderPlus feature rows are present because the source local sample directory
already contained validated AMRFinderPlus outputs; AMRFinderPlus execution itself
was not part of this controlled parallel native-runner test.

## Nextflow Result

- Status: PASS
- Processes succeeded: 9/9
- Duration: 11m 20s
- CPU hours: 2.8
- Output directory: `validation_runs/delftia_parallel_feature_only_v2`

## Native Runner Module Status

| Module | Status | Samples processed | Raw tables | Feature rows | Message |
| --- | --- | ---: | ---: | ---: | --- |
| ABRicate | PASS | 45 | 135 | 61 | Per-database sample-parallel workers=16 |
| IntegronFinder | PASS | 45 | 45 | 116 | Per-assembly parallel workers=16 |
| MLST | PASS | 45 | 1 | 0 | Per-assembly parallel workers=16; Delftia has no biological MLST scheme calls |
| MobileElementFinder | SKIPPED | 0 | 0 | 0 | Not enabled |

## PanR2 Feature Contract Result

- `schema_feature_rows`: 201
- `schema_unmatched_feature_rows`: 0
- `schema_invalid_feature_rows`: 0
- `schema_duplicate_feature_rows`: 0
- `database_setup_required_failures`: 0
- Dashboard generated: yes

Feature rows by table:

| Feature table | Rows |
| --- | ---: |
| `amr.features.tsv` | 19 |
| `amrfinderplus.features.tsv` | 24 |
| `vfdb.features.tsv` | 40 |
| `plasmidfinder.features.tsv` | 2 |
| `integronfinder.features.tsv` | 116 |
| `mlst.features.tsv` | 0 |
| `all_features.tsv` | 201 |

## Serial-vs-Parallel Comparison

The parallel output preserved the same core standardized feature calls as the
serial native-runner validation. The following files matched after sorting and
comparing the core contract/provenance-neutral columns `sample_id` through
`database_version`:

- `amr.features.tsv`
- `vfdb.features.tsv`
- `plasmidfinder.features.tsv`
- `integronfinder.features.tsv`
- `all_features.tsv`

Full rows can differ in source/provenance columns because Nextflow work
directories are different between serial and parallel runs.

## Interpretation

This validates the v0.3.0 parallel backend for the 45-genome Delftia feature
runner stage. The next validation step is to run the biologically richer
`Klebsiella pneumoniae` 100-genome validation with the same parallel backend,
GTDB-Tk disabled, and standard QC/annotation enabled.
