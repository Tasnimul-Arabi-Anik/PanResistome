# Troubleshooting

## CheckM2 Makes The Computer Unresponsive

CheckM2 can put heavy memory pressure on DIAMOND. Keep general pipeline parallelism high for feature runners, but cap CheckM2 separately:

```bash
--threads 16 --checkm2_threads 2
```

Use `--checkm2_threads 1` on low-memory machines. The `lowmem`, `desktop_parallel`, and `workstation` profiles set conservative CheckM2 caps.

For 100+ genome desktop validations on 16 GB RAM, prefer `--checkm2_threads 1` when you can tolerate a longer run. `--checkm2_threads 2` completed the Klebsiella 100 large-mode validation, but the DIAMOND phase still created noticeable memory pressure.

## CheckM2 Database Download Fails

PanResistome downloads the CheckM2 database automatically when `--checkm2_db` is not supplied. If a network or storage error occurs:

```bash
--checkm2_db_dir /path/to/reusable/checkm2_db_cache
```

Then resume the run. On shared systems, pre-download the database once and pass:

```bash
--checkm2_db /path/to/CheckM2_database/uniref100.KO.1.dmnd
```

## ABRicate Database Missing

The standard comprehensive profile requires `ncbi`, `vfdb`, and `plasmidfinder`. PanResistome runs `panr setup-db` before reporting and records the explicit setup action in:

```text
<sample>/panr2_inputs/manifest/abricate_database_setup_status.tsv
```

PanResistome force-refreshes the requested ABRicate databases by default. If you want to disable the refresh for offline, cached, or frozen-database reruns, use:

```bash
--panr2_update_abricate_db false
```

The default refresh uses `abricate-get_db --force` when that helper is available in the ABRicate environment, then re-indexes with `abricate --setupdb`. If setup still fails, inspect:

```text
<sample>/panr2_inputs/manifest/database_setup_status.tsv
```

Do not add `isfinder` to `--panr2_abricate_dbs` unless that database is installed and licensed for your use.

## AMRFinderPlus Is Slow

AMRFinderPlus runs one independent job per assembly. By default, PanResistome uses:

```bash
--amrfinderplus_jobs <threads> --amrfinderplus_threads_per_sample 1
```

For a 16-core desktop, this usually gives better throughput than one serial AMRFinderPlus loop. If memory or disk pressure is high, reduce the job count:

```bash
--amrfinderplus_jobs 4
```

For 300+ genome nucleotide runs, AMRFinderPlus can still take multiple hours because each sample may spend several minutes in `tblastn`. The 300-record Klebsiella large-mode validation therefore skipped AMRFinderPlus for the desktop-safe pass and treats AMRFinderPlus as a separate workstation/HPC or overnight benchmark.

Interrupted AMRFinderPlus runs can be resumed more safely: non-empty per-sample TSV outputs are reused instead of recomputed.

Per-sample status and logs are written under:

```text
<sample>/amrfinderplus/tables/amrfinderplus_sample_status.tsv
<sample>/amrfinderplus/raw/*.log
```

## One Genome Fails Download

NCBI records can change or temporarily disappear. The pipeline records failed accessions and continues with downloaded genomes:

```text
<sample>/sequence/failed_accessions.txt
<sample>/sequence/sequence_download_summary.csv
```

If too many fail, rerun later or regenerate the input with current NCBI Assembly records.

## MLST Produces No Calls

Some organisms do not have a matching PubMLST scheme in the `mlst` tool. PanResistome keeps the raw MLST output as run evidence but PanR2 suppresses placeholder features such as `ST_-`.

Expected outputs:

```text
<sample>/tool_results/mlst/raw/mlst.tsv
<sample>/panr2_inputs/features/mlst.features.tsv
```

An empty `mlst.features.tsv` can be valid for unsupported organisms.

## IntegronFinder Is Slow

Use the native parallel runner:

```bash
--panr2_native_feature_runners true \
--panr2_native_feature_runner_mode parallel
```

For a desktop-scale run, combine this with:

```bash
-profile conda,mamba,desktop_parallel
```

On 300+ genomes, IntegronFinder remains a visible runtime cost even with parallel workers. This is expected; inspect `pipeline_runtime_summary.tsv` to decide whether to lower or raise runner parallelism for the machine.

## ANI Is Slow On 300+ Genomes

FastANI all-vs-all can dominate wall time on 300+ genomes. For a first desktop-scale large-mode validation, disable ANI:

```bash
--run_ani false
```

When ANI is enabled in large-dataset mode, PanResistome now protects against accidental expensive all-vs-all runs:

```bash
--run_ani true --large_dataset true --ani_large_run_strategy auto
```

With the default `--ani_max_all_vs_all_genomes 200`, this writes an ANI run-status audit and per-genome placeholder ANI-cluster features instead of launching all-vs-all ANI above the threshold. Use `--ani_large_run_strategy all` only when you intentionally want the full ANI run, or `--ani_large_run_strategy skip` when you want a clear skip regardless of sample count.

Use full ANI later as a separate comparative-genomics pass, preferably with chunking, Mash/sketch pre-screening, representative-only genomes, skani, or an HPC/workstation profile.

## MobileElementFinder Fails

MobileElementFinder is opt-in because upstream parser behavior can vary by installation and assembly. Enable only when needed:

```bash
--panr2_run_mobileelementfinder true
```

If it fails, rerun without it and provide MobileElementFinder tables later as precomputed PanR2 inputs.

## MOB-suite Database Setup Fails

When `--run_mobsuite true` is enabled and no `--mobsuite_db` is supplied, PanResistome initializes a MOB-suite cache under:

```text
<outdir>/databases/mobsuite
```

The setup audit is:

```text
<sample>/mobsuite/mobsuite_database_setup_status.tsv
```

If a shared or preinitialized database is available, prefer:

```bash
--run_mobsuite true --mobsuite_db /path/to/mobsuite_db
```

For offline runs, the database should include the core MOB-suite files and `taxa.sqlite`.

## geNomad Database Download Fails

When `--run_genomad true` is enabled and no `--genomad_db` is supplied, PanResistome downloads/caches the geNomad database under:

```text
<outdir>/databases/genomad
```

The setup audit is:

```text
<sample>/prophage/genomad_database_setup_status.tsv
```

Before launching a geNomad run, use the fast readiness check:

```bash
python scripts/check_genomad_readiness.py \
  --db-dir /path/to/genomad_db \
  --out genomad_readiness.tsv
```

This check does not download the database. It reports whether `genomad` is available in the current environment and whether the database path already contains files.

For restricted networks or shared systems, pre-download the database once and pass:

```bash
--run_genomad true --genomad_db /path/to/genomad_db
```

The first validated auto-download attempt on this machine reached the `GENOMAD_PROPHAGE` process, but remained in first-run Conda/Mamba environment creation for about 17 minutes and did not reach database download. For now, treat geNomad as opt-in and prefer a cached Conda environment, container/Apptainer route, or explicitly supplied `--genomad_db` for real validations.

## Report Has Too Many Features

Use the large-dataset report safeguards first. Complete TSV outputs are still written, but report-facing matrices and co-occurrence summaries are capped:

```bash
--large_dataset true \
--report_mode compact \
--max_features_heatmap 150 \
--max_features_network 150
```

You can also use stricter PanR2 plotting controls:

```bash
--panr2_plot_style compact \
--panr2_label_max_length 30
```

For large validations, inspect the top-level dashboard and top findings before opening large heatmaps:

```text
<sample>/report/index.html
<sample>/panr2_inputs/report/top_findings.html
<sample>/panr2_inputs/manifest/report_controls.tsv
```

For same-contig/proximity evidence, large mode writes a capped report-facing table and a complete table:

```text
<sample>/panr2_inputs/cross_database/feature_proximity.tsv
<sample>/panr2_inputs/cross_database/feature_proximity_all.tsv
```
