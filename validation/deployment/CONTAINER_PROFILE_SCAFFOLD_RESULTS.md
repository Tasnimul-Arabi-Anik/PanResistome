# Container Profile Scaffold Validation

Date: 2026-05-12

Purpose: add and validate conservative Docker, Apptainer, and Singularity profile scaffolding for v0.4.0 deployment work without claiming full container support before an image/runtime smoke test exists.

## What Was Added

Profiles:

```text
docker
apptainer
singularity
```

Parameters:

```text
--container_image
--container_run_options
```

Readiness helper:

```text
scripts/check_container_readiness.py
```

## Configuration Validation

The profiles resolve with:

```text
conda.enabled=false
process.conda=false
process.container=<experimental image or --container_image>
docker.enabled=true, or apptainer.enabled=true, or singularity.enabled=true
```

This confirms the profiles are wired for container execution rather than mixed Conda/container execution.

## Runtime Availability On This Machine

Initial container runtime availability before local installation:

```text
docker: not found
apptainer: not found
singularity: not found
```

The readiness helper correctly reports this as a pre-run failure:

```text
status=FAIL_RUNTIME_MISSING
message=apptainer executable was not found in PATH.
```

For the real smoke test, Ubuntu `singularity-container` was installed from the
system package manager:

```text
singularity-ce version 4.1.1
```

## Real Singularity Smoke Test

Runtime pull/exec smoke:

```bash
env SINGULARITY_CACHEDIR=/tmp/panresistome_singularity_cache \
  singularity exec docker://alpine:3.19 sh -c 'echo singularity_smoke_ok && uname -m'
```

Result:

```text
singularity_smoke_ok
x86_64
```

PanResistome fixture smoke:

```bash
env NXF_SINGULARITY_CACHEDIR=/tmp/panresistome_singularity_cache \
  nextflow run main.nf \
    -profile test,singularity \
    --container_image docker://python:3.11-slim \
    --outdir validation_runs/container_singularity_test_profile_v3
```

Result:

```text
SEQUENCE_QC: PASS
COMBINED_QC: PASS
COLLECT_RESULTS: PASS
Pipeline completed.
```

This test validates that Nextflow can run PanResistome processes through the
Singularity profile on local fixture data. The profile required an explicit
repository bind mount because processes call scripts through `${baseDir}`.

## Local Docker Production Image Build

Docker was installed and started locally after the initial scaffold validation.

Runtime smoke:

```bash
sudo docker run --rm hello-world
```

Result:

```text
Hello from Docker!
```

Experimental image build:

```bash
sudo docker build -f containers/Dockerfile -t panresistome:experimental .
```

Result:

```text
Successfully built 4fa2a590cc52
Successfully tagged panresistome:experimental
```

The Dockerfile command checks found the expected executable entry points:

```text
abricate
mlst
integron_finder
mefinder
fetchm2
fetchM
seqkit
amrfinder
checkm2
fastANI
skani
quast.py
mash
mob_recon
genomad
kleborate
kaptive
ectyper
panr
```

Runtime sanity check:

```bash
sudo docker run --rm panresistome:experimental bash -lc \
  'set -e; python --version; mefinder --help >/dev/null; genomad --version; mob_recon --version; ectyper --version; panr --version'
```

Result:

```text
Python 3.10.20
geNomad, version 1.12.0
mob_recon 3.1.9
ectyper 2.0.0 running database version 1.0
PanR2 0.1.3-dev
```

Nextflow Docker-profile test using the local image:

```bash
sudo nextflow run main.nf \
  -profile test,docker \
  --container_image panresistome:experimental \
  --outdir /tmp/panresistome_container_test_output \
  -w /tmp/panresistome_container_test_work
```

Result:

```text
SEQUENCE_QC: PASS
COMBINED_QC: PASS
COLLECT_RESULTS: PASS
Pipeline completed.
Results saved to: /tmp/panresistome_container_test_output
```

The first attempt was launched from `/tmp` using an absolute `main.nf` path and
failed because the `test` profile resolved fixture paths relative to `/tmp`.
The passing run was launched from the repository root while keeping work and
output directories in `/tmp`.

Image size:

```text
content size: 7.45 GB
local Docker disk usage after build: 31.4 GB
```

## Real Biological Docker-Profile Validation

A two-genome Klebsiella run was executed through the Docker profile with the
local `panresistome:experimental` image. CheckM2, QUAST, AMRFinderPlus, and
GTDB-Tk were disabled to keep this as a fast container/runtime validation while
still exercising FetchM2 download, sequence QC, ANI/Mash pairwise logic,
PanResistome-native ABRicate/IntegronFinder/MLST runners, PanR2 comprehensive
analysis, handoff export, and final collection.

Command:

```bash
sudo nextflow run main.nf \
  --input /tmp/klebsiella_2_container_ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_2_docker_bio \
  -profile docker,large \
  --container_image panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 4 \
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_klebsiella_2_docker_bio_work
```

Result:

```text
Runtime: 2m58s
Nextflow processes: 14/14 succeeded
Downloaded genomes: 2
QC PASS: 2
Feature rows: 286
Databases seen: amr, mlst, plasmidfinder, vfdb
Unmatched feature rows: 0
Invalid feature rows: 0
Duplicate feature rows: 0
Database setup failures: 0
ABRicate ncbi/vfdb/plasmidfinder setup: PASS
PanR2 handoff report pages: generated
```

Feature rows by table, excluding headers:

```text
amr.features.tsv: 9
mlst.features.tsv: 34
plasmidfinder.features.tsv: 3
vfdb.features.tsv: 240
all_features.tsv: 286
```

Native-runner merge audit:

```text
abricate: PASS, expected_raw_tables=6, observed_raw_tables=6, feature_rows=252
integronfinder: PASS, expected_raw_tables=2, observed_raw_tables=2, feature_rows=0
mlst: PASS, expected_raw_tables=1, observed_raw_tables=1, feature_rows=16
mobileelementfinder: SKIPPED, not requested
```

The run generated the expected PanR2 handoff HTML pages, including:

```text
panr2_inputs/report/index.html
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/metadata_quality_and_bias.html
panr2_inputs/report/database_setup_and_contract.html
panr2_inputs/report/report_controls.html
panr2_inputs/report/bioproject_bias.html
panr2_inputs/report/amrfinder_abricate_concordance.html
panr2_inputs/report/cross_database_interpretation.html
```

The only warning was outside the biological pipeline result: the optional
workflow `onComplete` runtime/resource summary hook tried to execute `python`
from the host environment, but this sudo environment only had `python3`. The
pipeline itself completed and produced valid outputs. The hook has since been
updated to prefer `python3` and fall back to `python`.

## GHCR Public Pull Attempt

An unauthenticated public pull was attempted:

```bash
sudo docker pull ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

The pull started successfully and downloaded layers without requiring a GitHub
login, so authentication was not the observed blocker. The pull was stopped
after roughly ten minutes because the network was too slow for the large image
on this machine. Full remote pull validation remains open.

Important observations:

- The first all-in-one image build exposed a real Biopython conflict: MobileElementFinder requires an older Biopython stack than IntegronFinder. The image now installs MobileElementFinder separately in `mobileelementfinder_env` and keeps PanR2/ABRicate/IntegronFinder/MLST in `panr2_container_env`.
- The MobileElementFinder executable is `mefinder`, not `mobileelementfinder`.
- MobileElementFinder imports `pkg_resources`, so `mobileelementfinder_env` pins `setuptools<81`.
- geNomad installs successfully but is a large layer because it pulls TensorFlow/MMseqs2-related dependencies.
- ECTyper starts, but its optional ~900 MB species-ID sketch download timed out during build. ECTyper should not be documented as fully self-contained until a dedicated ECTyper database/readiness validation passes.
- Docker currently requires `sudo` on this host. A normal Nextflow `-profile docker` run should wait until the user is in the `docker` group and has opened a new login/session; running Nextflow itself with `sudo` is not recommended because it can create root-owned work files.

## Validation Scope

Validated:

```text
Nextflow profile syntax/config resolution
Conda disabled in container profiles
Container image assignment is present
Container readiness helper reports missing runtime/image/database paths
Standard non-container test profile still passes
Singularity runtime installed and executable
Singularity can pull and execute a public Docker image
PanResistome test profile completes through Singularity with docker://python:3.11-slim
Docker runtime installed and can execute hello-world
Experimental PanResistome all-in-one image builds locally
Experimental image command and runtime sanity checks pass for core optional runners
Nextflow `test,docker` profile completes with the local image when launched from the repository root
Two-genome Klebsiella biological run completes through the Docker profile with the local image
```

Not validated:

```text
Apptainer runtime execution
Containerized database mounts
Normal-user Nextflow Docker profile on this host
Full PanResistome production image pull from GHCR
ECTyper species-ID sketch readiness
```

The documented GHCR image started pulling without login on 2026-05-12, but the
full pull was not completed because of local network speed and image size. If
future pulls return `DENIED`, the GitHub package visibility should be checked.
The next deployment blocker is full remote image pull validation, not Nextflow
Docker/Singularity profile wiring.

## Production Image Scaffold

Added:

```text
containers/Dockerfile
.github/workflows/container.yml
.dockerignore
```

The image definition creates the existing PanResistome tool environments inside
one image and exposes their command-line tools on `PATH`, matching the
`conda.enabled=false` behavior of the container profiles. GTDB-Tk is excluded
from the first image target because it remains an opt-in large-database mode.

The readiness helper now supports a real pull/exec check:

```bash
python scripts/check_container_readiness.py \
  --runtime singularity \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --pull-test \
  --out container_readiness.tsv
```

The image build now passes locally as `panresistome:experimental`, and the
local image can run both the Nextflow test profile and a two-genome biological
Klebsiella validation through Docker. The GitHub Actions workflow should still
publish a fully pullable GHCR image before Docker, Apptainer, or Singularity
profiles are advertised as a low-hassle remote-user route.

## Next Required Validation

A full remote-user container validation requires a fully pulled PanResistome
image from GHCR or an Apptainer/Singularity cache on the target machine.

Recommended first command on that host:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /path/to/checkm2,/path/to/genomad \
  --out container_readiness.tsv
```

Then run a small fixture or 5-genome smoke before any large biological validation.
