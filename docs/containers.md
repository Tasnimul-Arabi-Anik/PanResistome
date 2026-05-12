# Container execution notes

PanResistome currently ships validated Conda/Mamba workflows. Docker and Apptainer/Singularity profiles are present as v0.4.0 deployment scaffolding. A Singularity fixture smoke test has passed, but full production-image and real-data container validation are still pending.

## Design goal

The container strategy should preserve the current architecture:

```text
PanResistome = tool execution, QC, database setup, feature export
PanR2        = standardized feature analysis and reporting
```

Each heavy tool may need its own image or carefully separated environment. A single monolithic image can be convenient, but it can also become difficult to rebuild and debug.

## Database paths

Containerized runs must mount large external databases into the container:

- CheckM2 database
- GTDB-Tk data, if enabled
- geNomad database, if enabled
- Kaptive database, if enabled
- authorized ISfinder FASTA/BLAST database, if enabled

Keep these databases outside `work/` and outside generated result folders.

## Readiness check

Before launching a container run, check the runtime, image string, and host database paths:

```bash
python scripts/check_container_readiness.py \
  --runtime apptainer \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --database-paths /path/to/checkm2,/path/to/genomad \
  --out container_readiness.tsv
```

The checker does not pull images or run Nextflow. It prevents the common failure where a workflow starts before the runtime, image, or database mounts are actually ready.

## Experimental profiles

Profiles:

```text
docker
apptainer
singularity
```

Common parameters:

```text
--container_image
--container_run_options
```

The profiles disable Conda for all processes and assign the supplied image to each process. They also bind the repository path into the container because workflow processes call helper scripts through `${baseDir}`. If no image is supplied, they use the documented experimental placeholder image name. Do not use these profiles for production until the image has been built and validated.

The scaffold and Singularity fixture-smoke validation status is documented at `validation/deployment/CONTAINER_PROFILE_SCAFFOLD_RESULTS.md`.

## Practical first container targets

1. Keep the existing `-profile test` workflow working without containers.
2. Small Singularity fixture smoke test for local fixtures. Completed with `docker://python:3.11-slim`.
3. Build or publish a PanResistome image with the real tool stack.
4. Validate the standard comprehensive command without GTDB-Tk.
5. Document database mounts for CheckM2 and optional heavy databases.
6. Only then advertise Docker/Apptainer profiles as supported.

The current placeholder image is not publicly pullable from GHCR:

```text
docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental
DENIED: requested access to the resource is denied
```

Until that is fixed, pass a known image explicitly with `--container_image` for
smoke testing.

## geNomad-specific note

geNomad is a high-priority container candidate because the 5-genome auto-download validation attempt reached `GENOMAD_PROPHAGE` but remained in first-run Conda/Mamba environment creation for about 17 minutes before database download began. A useful container validation should prove:

- `genomad --version` works inside the image.
- `genomad download-database` can write to a mounted database path.
- The mounted database path can be reused with `--genomad_db`.
- A 2-10 genome run produces raw geNomad output and `prophage.features.tsv`.
- `panr2_inputs/features/all_features.tsv` includes prophage/geNomad rows with zero unmatched, invalid, or duplicate feature rows.

Until that is validated, geNomad should remain opt-in and described as table-analysis-ready but runner-validation-pending.

For a container or host-module geNomad validation, use `-profile genomad_host` so PanResistome does not try to create `envs/genomad.yaml` for the geNomad processes:

```bash
nextflow run main.nf \
  -profile conda,mamba,genomad_host \
  --run_genomad true \
  --genomad_db /databases/genomad
```

This only changes the geNomad processes. Other Conda-backed processes still use the normal selected profile.

## Example future pattern

This is a planning example, not a validated command:

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_container \
  -profile apptainer,large \
  --container_image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_amrfinderplus true \
  --checkm2_db /databases/checkm2/uniref100.KO.1.dmnd
```

Do not assume host database paths are visible inside a container. Mounts must expose the same paths that the Nextflow command passes to PanResistome.
