# HPC execution notes

PanResistome is still validated primarily with local Conda/Mamba profiles. This page documents the intended HPC pattern for v0.4.0 planning; treat scheduler/container profiles as experimental until a cluster validation is recorded.

## Recommended starting point

For a workstation or login-node dry run, first confirm that the test profile works:

```bash
nextflow run main.nf -profile test
```

For real data on shared infrastructure, prefer a capped CheckM2 thread count and the large-report safeguards:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_hpc_dryrun \
  -profile conda,mamba,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --threads 16 \
  --checkm2_threads 4
```

## Resource guidance

- Keep `--checkm2_threads` lower than global `--threads`; CheckM2/DIAMOND can be memory-heavy.
- Keep GTDB-Tk disabled unless the reference database is installed on shared storage and the process has enough memory.
- Use `--large_dataset true` or `-profile large` for large feature matrices. Complete TSVs are still exported, but HTML/report-facing summaries are capped.
- AMRFinderPlus uses per-assembly parallel jobs. Tune `--amrfinderplus_jobs` separately from `--amrfinderplus_threads_per_sample` if a scheduler node has many cores but limited memory.
- Put large databases under stable shared paths rather than inside transient Nextflow work directories.
- Use `-resume` after transient scheduler or network failures.

## Scheduler profile status

No stable SLURM profile is advertised yet. A future profile should define executor settings, queue/account placeholders, process labels, and database mount/cache guidance without changing the default local Conda/Mamba workflows.

## Container profile status

Experimental Docker, Apptainer, and Singularity profiles exist for v0.4.0 deployment testing. On HPC, prefer Apptainer/Singularity with shared database mounts. Validate readiness before submitting a job:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /shared/db/checkm2,/shared/db/genomad \
  --out container_readiness.tsv
```

The profiles are not yet release-grade. They need a small fixture smoke test, then a real-data validation, before replacing the Conda/Mamba route in documentation.

## Outputs to inspect first

After a run, start with:

```text
<outdir>/<organism>/report/index.html
<outdir>/<organism>/panr2_inputs/report/panr2_handoff_index.html
<outdir>/<organism>/panr2_inputs/manifest/report_controls.tsv
<outdir>/<organism>/panr2_inputs/manifest/database_setup_status.tsv
<outdir>/<organism>/panr2_inputs/manifest/schema_validation_summary.txt
<outdir>/pipeline_runtime_summary.tsv
<outdir>/pipeline_runtime_tasks.tsv
```
