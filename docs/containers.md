# Container execution notes

PanResistome currently ships validated Conda/Mamba workflows. Docker and
Apptainer/Singularity profiles are present as v0.4.0 deployment scaffolding. A
Singularity fixture smoke test has passed, a local Docker build/runtime smoke
test of the experimental all-in-one image passes, a two-genome Klebsiella
biological run completed through the Docker profile with the local image,
geNomad database download and runner execution work through Docker with a
mounted database path, a 100-record Klebsiella large-mode Docker run completed
with the local image and clean PanR2 feature-contract output, the public GHCR
image has completed an unauthenticated Docker pull plus two-genome
geNomad-enabled biological validation and five-genome geNomAD positive-call
validation with low-memory `--genomad_splits`, and a five-genome Acinetobacter
CheckM2-on comprehensive validation with automatic CheckM2/geNomad database
downloads. Singularity CE has completed GHCR pull/exec validation plus
two-genome geNomad-enabled and 100-record large-mode biological validations.
Apptainer validation remains pending, but the Apptainer profile now follows the
same read-only database defaults as Singularity.

## Current validation status

| Route | Image/source | Validation status | Notes |
| --- | --- | --- | --- |
| Conda/Mamba | Local environments | Validated | Default public route for non-container runs. |
| Docker, local image | `panresistome:experimental` | Validated for two-genome and 100-record Klebsiella biological runs | Proves the image contents and Docker profile can run the biological workflow. |
| Docker, GHCR image | `ghcr.io/tasnimul-arabi-anik/panresistome:experimental` | Validated for unauthenticated pull, two-genome geNomad-enabled biological run, five-genome geNomAD positive-call run, five-genome Acinetobacter CheckM2-on comprehensive run, and 100-record large-mode biological run | The 100-record run completed 12/12 processes in 21m35s with 11,488 standardized feature rows and zero unmatched, invalid, or duplicate feature rows. The focused geNomAD scale run processed 5/5 genomes and produced 23 prophage/geNomAD feature rows. The Acinetobacter comprehensive run processed 5/5 genomes with automatic CheckM2 and geNomad database download, 630 PanR2 feature rows, and zero unmatched, invalid, or duplicate feature rows. |
| Singularity CE, GHCR image | `docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental` | Validated for pull/exec, two-genome geNomad-enabled run, and 100-record large-mode run | Recommended HPC-style route until Apptainer is validated. |
| Apptainer, GHCR image | `docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental` | Pending | Profile exists, but local biological validation has not been recorded because Apptainer is not installed on the validation host. Ubuntu 24.04 packages on this host provide Singularity CE, not the Apptainer runtime. Installing Apptainer requires an explicit administrator decision to add an external package source or install from another trusted route. |
| Non-sudo Docker | Local Docker daemon | Configured on validation host | User `anik` was added to the `docker` group on 2026-05-13 and the socket group was set to `docker`; a fresh login session is required before normal shells can run `docker ps` without `sudo`. |

## Design goal

The container strategy should preserve the current architecture:

```text
PanResistome = tool execution, QC, database setup, feature export
PanR2        = standardized feature analysis and reporting
```

Each heavy tool may need its own image or carefully separated environment. The experimental all-in-one image keeps separate Conda environments for conflicting tools while exposing their command-line entry points on `PATH`. This is convenient for users, but the image is large and slow to build locally.

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

By default, the checker does not pull images or run Nextflow. It prevents the common failure where a workflow starts before the runtime, image, or database mounts are actually ready.

To verify that the runtime can also pull and execute the image, add `--pull-test`:

```bash
python scripts/check_container_readiness.py \
  --runtime singularity \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --pull-test \
  --pull-test-timeout 7200 \
  --out container_readiness.tsv
```

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

The profiles disable Conda for all processes and assign the supplied image to
each process. They also bind the repository path into the container because
workflow processes call helper scripts through `${baseDir}`. If no image is
supplied, they use the documented experimental GHCR image name. Apptainer and
Singularity profiles default `--panr2_update_abricate_db false` because
SIF/container images are read-only and should use the frozen ABRicate databases
already packaged in the image. Docker can refresh ABRicate databases in writable
work directories, but Singularity/Apptainer users should treat the image
databases as frozen unless a writable database strategy is explicitly mounted.

The scaffold, Singularity fixture-smoke, Docker biological, geNomad Docker, and 100-record Docker validation status is documented in `validation/deployment/CONTAINER_PROFILE_SCAFFOLD_RESULTS.md`, `validation/deployment/DOCKER_REMOTE_USER_VALIDATION_RESULTS.md`, `validation/deployment/GHCR_DOCKER_100_VALIDATION_RESULTS.md`, `validation/deployment/SINGULARITY_GHCR_VALIDATION_RESULTS.md`, and `validation/optional_runner_biological/GENOMAD_POSITIVE_CALL_VALIDATION_RESULTS.md`.

## Practical first container targets

1. Keep the existing `-profile test` workflow working without containers.
2. Small Singularity fixture smoke test for local fixtures. Completed with `docker://python:3.11-slim`.
3. Build or publish a PanResistome image with the real tool stack. Local Docker build/runtime sanity has passed, and the GHCR image can be pulled and executed.
4. Validate the standard comprehensive command without GTDB-Tk through a normal user Docker/Apptainer profile. Two-genome and 100-record biological Docker runs with the local image have passed, two-genome geNomad-enabled and 100-record large-mode runs have passed with the pulled GHCR image through Docker, and the same small biological path plus a 100-record large-mode path have passed through Singularity CE.
5. Document database mounts for CheckM2 and optional heavy databases.
6. Only then advertise Docker/Apptainer profiles as supported.

The current GHCR image name is:

```text
ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

For exact reproducibility, prefer pinning an immutable digest after selecting a
validated image:

```bash
docker image inspect ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --format '{{index .RepoDigests 0}}'
```

Then run with the digest string, for example:

```text
--container_image ghcr.io/tasnimul-arabi-anik/panresistome@sha256:<digest>
```

The PanR2 handoff now records the selected image string and any embedded
`sha256` digest in `panr2_inputs/manifest/reproducibility_manifest.json`.

An unauthenticated `docker pull` completed successfully on 2026-05-13, which means a GitHub login was not required for the public experimental image on the validation host. If users see `DENIED`, the GHCR package visibility should be checked in GitHub package settings.

Singularity CE 4.1.1 also pulled and executed the same GHCR image with:

```bash
singularity exec docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental true
```

The first conversion took about 1h15m on the validation host. Subsequent
Nextflow biological validations succeeded once the image was available,
including a 100-record Klebsiella large-mode run. Users should set
`NXF_SINGULARITY_CACHEDIR` to a persistent shared cache on HPC so the large SIF
conversion is not repeated in each work directory.

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/cache/panresistome/singularity
```

On slow networks, Nextflow's default Singularity pull timeout can expire during
the first GHCR-to-SIF conversion. Use a small extra config file when needed:

```nextflow
singularity.pullTimeout = '2h'
```

Then include it with `nextflow -c /path/to/singularity_long_pull.config run ...`.

Normal non-sudo Docker depends on local socket permissions. On Linux, that often
means adding the user to the `docker` group and then opening a new login
session. Docker group membership grants broad Docker/root-equivalent access, so
it should be enabled explicitly by the machine owner or administrator. On the
validation host, user `anik` was added to the `docker` group and a root-checked
socket permission test confirmed `root:docker` ownership with mode `660`; the
active desktop/login session still needs to be refreshed before ordinary
non-sudo Docker commands inherit that group.

## Image build and publication

The experimental all-in-one image definition is:

```text
containers/Dockerfile
```

It pre-creates the existing PanResistome tool environments inside the image and
exposes their `bin` directories on `PATH`, allowing the Nextflow container
profiles to run with `conda.enabled=false`. GTDB-Tk is intentionally excluded
from the first image target because its database and runtime footprint are too
large for the default deployment route.

Build locally on a Docker/Podman-capable machine:

```bash
docker build -f containers/Dockerfile -t panresistome:experimental .
```

Local build/runtime status on 2026-05-12:

```text
Docker install and hello-world: PASS
panresistome:experimental build: PASS
Dockerfile command checks: PASS; current image builds must also load the CheckM2 packaged model
runtime checks: mefinder, genomad, mob_recon, ectyper, panr all start
Nextflow test,docker profile with local image: PASS
image content size: about 7.45 GB
local Docker disk usage during build: about 31.4 GB
```

Two-genome biological Docker-profile validation using the local image also passed on 2026-05-12:

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
14/14 Nextflow processes succeeded
2/2 genomes downloaded and QC-passed
286 standardized feature rows
databases detected: amr, vfdb, plasmidfinder, mlst
unmatched/invalid/duplicate feature rows: 0/0/0
ABRicate ncbi/vfdb/plasmidfinder setup: PASS
PanR2 handoff HTML pages: generated
```

100-record Klebsiella large-mode Docker validation using the local image passed
on 2026-05-13:

```text
12/12 Nextflow processes succeeded
runtime: 22m16s
CPU hours: 5.0
input records: 100
downloaded/analyzed genomes: 99
failed accession: GCF_055382775.1
complete standardized feature rows: 11,488
unmatched/invalid/duplicate feature rows: 0/0/0
ABRicate ncbi/vfdb/plasmidfinder setup: PASS
native runner merge audit: PASS
PanR2 handoff HTML pages: generated
```

geNomad database download and runner validation also passed through Docker with
a mounted database path:

```text
geNomad database v1.9 downloaded and extracted inside Docker
two-genome geNomad-enabled Nextflow run: 16/16 processes succeeded
geNomad process: PASS
prophage.features.tsv: header-only for this tiny dataset
feature-contract validation: zero unmatched/invalid/duplicate rows
```

The header-only prophage table means the runner and PanR2 contract path are
valid, but positive geNomad feature calls still need a prophage-rich validation
dataset.

Two implementation details matter for reproducibility:

- MobileElementFinder is installed in a separate `mobileelementfinder_env`. Its CLI is `mefinder`, not `mobileelementfinder`, and the environment pins `setuptools<81` because the upstream package still imports `pkg_resources`.
- PanR2, ABRicate, IntegronFinder, and MLST are installed in `panr2_container_env` without MobileElementFinder to avoid the Biopython version conflict between MobileElementFinder and IntegronFinder.

ECTyper installed and `ectyper --version` starts, but its optional species-identification sketch download timed out during image build. Do not treat ECTyper species-ID data as fully bundled until a dedicated ECTyper database/readiness validation passes.

GitHub Actions workflow:

```text
.github/workflows/container.yml
```

The workflow publishes:

```text
ghcr.io/tasnimul-arabi-anik/panresistome:experimental
ghcr.io/tasnimul-arabi-anik/panresistome:<git-tag>
ghcr.io/tasnimul-arabi-anik/panresistome:sha-<commit>
```

After the image is pushed, the workflow runs the same command-availability and
runtime sanity checks used locally. The smoke test checks the core CLI entry
points and verifies that `mefinder`, geNomad, MOB-suite, ECTyper, and PanR2 can
start inside the published image.

After each GHCR push, confirm the package visibility is public in GitHub
package settings. If the package is private, Docker/Singularity/Apptainer pulls
will fail with `DENIED: requested access to the resource is denied`.

The validated Singularity readiness pattern is:

```bash
python scripts/check_container_readiness.py \
  --runtime singularity \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --pull-test \
  --pull-test-timeout 7200 \
  --out container_readiness.tsv

env NXF_SINGULARITY_CACHEDIR=/tmp/panresistome_singularity_cache \
  nextflow run main.nf \
    -profile singularity,large \
    --container_image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
    --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
    --outdir validation_runs/klebsiella_100_singularity_large \
    --analysis_profile comprehensive \
    --qc_filter true \
    --run_gtdbtk false \
    --run_checkm2 false \
    --run_quast false \
    --run_ani false \
    --run_mash true \
    --run_amrfinderplus false
```

## geNomad-specific note

geNomad is a high-priority container candidate because the 5-genome auto-download validation attempt reached `GENOMAD_PROPHAGE` but remained in first-run Conda/Mamba environment creation for about 17 minutes before database download began. A useful container validation should prove:

- `genomad --version` works inside the image.
- `genomad download-database` can write to a mounted database path.
- The mounted database path can be reused with `--genomad_db`.
- A 2-10 genome run produces raw geNomad output and `prophage.features.tsv`.
- A prophage-rich dataset produces positive prophage/geNomad rows with zero unmatched, invalid, or duplicate feature rows.

geNomad remains opt-in. The Docker route now validates database download,
database mounting, and runner execution, but a prophage-rich dataset is still
needed before claiming positive biological prophage feature validation.

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
