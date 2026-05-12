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
```

Not validated:

```text
Apptainer runtime execution
Containerized database mounts
Real-data container run
Normal-user Nextflow Docker profile on this host
Full PanResistome production image from GHCR
ECTyper species-ID sketch readiness
```

The documented placeholder image is not currently publicly pullable:

```text
docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental
GET https://ghcr.io/token?... DENIED: requested access to the resource is denied
```

That means the next deployment blocker is image publication/build validation,
not Singularity profile wiring.

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
local image can run the Nextflow test profile through Docker. The GitHub Actions
workflow should still publish a pullable GHCR image before Docker, Apptainer, or
Singularity profiles are advertised as a low-hassle remote-user route.

## Next Required Validation

A full real-data container validation requires a built/pullable PanResistome
image with the required tools installed.

Recommended first command on that host:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /path/to/checkm2,/path/to/genomad \
  --out container_readiness.tsv
```

Then run a small fixture or 5-genome smoke before any large biological validation.
