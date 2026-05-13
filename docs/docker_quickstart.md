# Docker Quickstart

This page is the shortest path for users who want PanResistome to manage the
tool stack without solving Conda/Mamba environments locally.

## What Docker Solves

Docker avoids most local tool installation work:

```text
No user-managed Conda/Mamba solving
No manual PanR2/ABRicate/IntegronFinder/MLST installation
No manual ABRicate ncbi/vfdb/plasmidfinder setup in comprehensive mode
Same container image across machines
```

The remaining requirements are:

```text
Docker installed and usable by the current user
Nextflow installed on the host
Internet access for genome/database downloads
Enough disk space for the image, work directory, and results
Mounted paths for large optional databases such as geNomad, GTDB-Tk, CheckM2, or ISfinder
```

## Image

The experimental GHCR image is:

```text
ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

If the GHCR package remains public, users do not need a GitHub login to pull it.
The validation host confirmed that unauthenticated pull completes successfully,
but the image is large and first pull can still be slow on weak networks.

Pull explicitly:

```bash
docker pull ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

If this reports `DENIED`, check GitHub package visibility.

## Docker Permission

Check whether Docker works without `sudo`:

```bash
docker ps
```

If this fails with a Docker socket permission error, ask the machine
administrator to configure Docker access. On Linux this often means adding the
user to the `docker` group and opening a new login session. Docker group
membership grants broad Docker/root-equivalent access, so this should be an
explicit administrative decision.

## Main Workflow

Example 100-record validation-style run:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_100_docker_large \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 8 \
  --fetchm2_download_workers 2
```

This command keeps the validation desktop-safe by disabling CheckM2, QUAST,
ANI, AMRFinderPlus, and GTDB-Tk. It still exercises FetchM2, sequence QC, Mash,
ABRicate ncbi/vfdb/plasmidfinder, IntegronFinder, MLST, PanR2 comprehensive
reporting, large-mode report caps, and PanR2 handoff export.

The same pattern completed locally with the built image:

```text
100 input records
99 genomes downloaded/analyzed
11,488 standardized feature rows
0 unmatched / 0 invalid / 0 duplicate feature rows
12/12 Nextflow processes succeeded
runtime: 22m16s
```

## geNomad With Docker

Download the geNomad database into a mounted host directory:

```bash
mkdir -p /path/to/genomad_db_parent

docker run --rm \
  -v /path/to/genomad_db_parent:/genomad_db \
  ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  bash -lc 'genomad download-database /genomad_db'
```

Run PanResistome with the mounted database:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_genomad_docker \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --container_run_options '-v /path/to/genomad_db_parent:/path/to/genomad_db_parent' \
  --analysis_profile comprehensive \
  --run_gtdbtk false \
  --run_genomad true \
  --genomad_use_host_env true \
  --genomad_db /path/to/genomad_db_parent/genomad_db
```

The Docker geNomad path has been validated for database download, database
mounting, runner execution, and clean PanR2 feature-contract export. The pulled
GHCR image also completed a two-genome geNomad-enabled biological run with 286
standardized feature rows and zero unmatched, invalid, or duplicate feature
rows. The small two-genome validation produced a header-only
`prophage.features.tsv`, so a prophage-rich dataset is still needed to validate
positive prophage calls.

## Output To Check

Start with:

```text
<outdir>/<organism>/report/index.html
<outdir>/<organism>/panr2_inputs/report/panr2_handoff_index.html
<outdir>/<organism>/panr2_inputs/features/all_features.tsv
<outdir>/<organism>/panr2_inputs/manifest/schema_validation_summary.txt
<outdir>/<organism>/panr2_inputs/manifest/database_setup_status.tsv
<outdir>/<organism>/panr2_inputs/manifest/native_runner_merge_audit.tsv
```

For a clean run, feature-contract validation should report:

```text
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

## Singularity/Apptainer Note

HPC users usually prefer Singularity or Apptainer instead of Docker. The same
GHCR image has been validated with Singularity CE for a two-genome
geNomad-enabled workflow. Use a persistent cache before first run:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/cache/panresistome/singularity
```

Then use `-profile singularity,large` and keep optional database paths mounted
with `--container_run_options`. The Singularity/Apptainer profiles default to
frozen ABRicate databases with `--panr2_update_abricate_db false` because SIF
images are read-only.

## Current Caveats

```text
The image is large; first pull can be slow.
GTDB-Tk remains external/opt-in because its database is very large.
ISfinder still requires a user-supplied authorized FASTA.
CheckM2 and AMRFinderPlus work in the pipeline, but the 100-genome Docker scale validation intentionally disabled them for desktop safety.
Normal Docker use requires Docker socket permission or sudo.
```

More detail is available in:

```text
docs/containers.md
validation/deployment/DOCKER_REMOTE_USER_VALIDATION_RESULTS.md
validation/deployment/CONTAINER_PROFILE_SCAFFOLD_RESULTS.md
```
