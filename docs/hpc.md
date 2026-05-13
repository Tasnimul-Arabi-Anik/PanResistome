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

Experimental Docker, Apptainer, and Singularity profiles exist for v0.4.0
deployment testing. A local Singularity fixture smoke test has passed with
`docker://python:3.11-slim`, and Docker validation has passed for two-genome,
geNomad-enabled, and 100-record Klebsiella workflows using the built
`panresistome:experimental` image. The public GHCR image has also been pulled
successfully and used for two-genome geNomad-enabled biological validation
through both Docker and Singularity CE, and Singularity CE has validated the
GHCR image on a 100-record large-mode workflow. Apptainer biological validation
remains pending. On HPC, prefer Apptainer/Singularity with shared database
mounts. Validate readiness before submitting a job:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /shared/db/checkm2,/shared/db/genomad \
  --out container_readiness.tsv
```

The profiles are still v0.4.0 deployment-stage paths rather than the only
documented route. Docker and Singularity now have biological validation
evidence; Apptainer still needs a small fixture smoke test and real-data
validation before it should be advertised as validated.

Add `--pull-test --pull-test-timeout 7200` on a login/development node when
image pulls are allowed.
Avoid running large database downloads or full biological validations on a login
node; use the pull test only to confirm that the runtime and image are usable.

For Singularity/Apptainer, use a persistent cache path before running Nextflow:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/cache/panresistome/singularity
```

The first GHCR-to-SIF conversion took about 1h15m on the validation host, and a
100-record large-mode Singularity validation completed after the image was
available. A shared cache prevents every run directory from materializing its own
copy of the large image. Apptainer and Singularity profiles default
`--panr2_update_abricate_db false` because the packaged ABRicate databases are
read-only inside SIF images. Use Docker or a future writable database strategy
when a forced ABRicate database refresh is required.

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
