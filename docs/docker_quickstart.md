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
Mounted paths for large optional databases such as GTDB-Tk or ISfinder; CheckM2 and geNomad can auto-download into the output cache when enabled, but shared predownloaded paths are still useful on restricted networks
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

After group membership is changed, the current shell usually still cannot use
Docker without `sudo`. Log out and back in, or open a fresh login session, then
rerun `docker ps`.

On the PanResistome validation host, user `anik` has been added to the `docker`
group and the Docker socket has been configured with group `docker`. A fresh
login session is still required before normal non-sudo Docker commands inherit
that permission.

## Main Workflow

Recommended stable comprehensive run for a remote Docker/GHCR user:

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

This command keeps GTDB-Tk, DefenseFinder, MobileElementFinder, and ISfinder
out of the first-pass route. It exercises FetchM2, sequence QC, CheckM2, QUAST,
ANI/skani, Mash, AMRFinderPlus, geNomad/prophage, ABRicate
ncbi/vfdb/plasmidfinder, IntegronFinder, MLST, PanR2 comprehensive reporting,
and PanR2 handoff export.

The same route passed a 5-genome `Acinetobacter pittii` Docker/GHCR validation
on 2026-05-16 with no `--checkm2_db` argument, `--threads 4`, and
`--checkm2_threads 1`. That run auto-downloaded CheckM2 and geNomad databases,
produced 5/5 CheckM2 rows, 5/5 combined QC PASS calls, 630 PanR2 feature rows,
and zero unmatched, invalid, or duplicate feature rows.

For a desktop-safe scale run, reduce the heavy modules intentionally:

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

The same pattern completed with both the local built image and the pulled GHCR
image:

```text
100 input records
99 genomes downloaded/analyzed
11,488 standardized feature rows
0 unmatched / 0 invalid / 0 duplicate feature rows
12/12 Nextflow processes succeeded
runtime: 22m16s
```

The GHCR-image validation completed the same workflow in 21m35s with 12/12
Nextflow processes succeeded, 11,488 standardized feature rows, and zero
unmatched, invalid, or duplicate feature rows.

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
geNomad-enabled workflow and a 100-record large-mode workflow. Use a persistent
cache before first run:

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
The 100-genome Docker scale validation intentionally disabled CheckM2 and AMRFinderPlus for desktop safety.
CheckM2 now has GHCR model-load, standalone prediction, Docker-profile QC fixture evidence, and a 5-genome Acinetobacter Docker/GHCR comprehensive PASS with automatic database download.
Normal Docker use requires Docker socket permission or sudo; Docker group changes require a new login session.
GHCR Docker is validated for both two-genome geNomad biological execution and 100-record large-mode scale execution.
```

More detail is available in:

```text
docs/containers.md
validation/deployment/DOCKER_REMOTE_USER_VALIDATION_RESULTS.md
validation/deployment/GHCR_DOCKER_100_VALIDATION_RESULTS.md
validation/deployment/CONTAINER_PROFILE_SCAFFOLD_RESULTS.md
```
