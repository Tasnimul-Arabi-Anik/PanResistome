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

Container runtimes were not available in `PATH`:

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

## Validation Scope

Validated:

```text
Nextflow profile syntax/config resolution
Conda disabled in container profiles
Container image assignment is present
Container readiness helper reports missing runtime/image/database paths
Standard non-container test profile still passes
```

Not validated:

```text
Docker image build or pull
Apptainer/Singularity image pull
Containerized Nextflow process execution
Containerized database mounts
Real-data container run
```

## Next Required Validation

A real container validation requires a host with Docker, Apptainer, or Singularity installed and a built/pullable PanResistome image.

Recommended first command on that host:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /path/to/checkm2,/path/to/genomad \
  --out container_readiness.tsv
```

Then run a small fixture or 5-genome smoke before any large biological validation.
