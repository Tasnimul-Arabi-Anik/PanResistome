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

The standard comprehensive profile requires `ncbi`, `vfdb`, and `plasmidfinder`. PanResistome runs `panr setup-db` before reporting. If setup still fails, inspect:

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

## MobileElementFinder Fails

MobileElementFinder is opt-in because upstream parser behavior can vary by installation and assembly. Enable only when needed:

```bash
--panr2_run_mobileelementfinder true
```

If it fails, rerun without it and provide MobileElementFinder tables later as precomputed PanR2 inputs.

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
