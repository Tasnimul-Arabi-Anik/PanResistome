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
```

Not validated:

```text
Docker image build or pull
Apptainer runtime execution
Containerized database mounts
Real-data container run
Full PanResistome production image
```

The documented placeholder image is not currently publicly pullable:

```text
docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental
GET https://ghcr.io/token?... DENIED: requested access to the resource is denied
```

That means the next deployment blocker is image publication/build validation,
not Singularity profile wiring.

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
