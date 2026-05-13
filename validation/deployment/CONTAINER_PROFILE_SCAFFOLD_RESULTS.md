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

## 100-Record Docker Large-Mode Validation

A 100-record Klebsiella validation was executed through the Docker profile with
the local `panresistome:experimental` image. CheckM2, QUAST, ANI, AMRFinderPlus,
and GTDB-Tk were disabled for desktop safety; Mash, sequence QC, native feature
runners, PanR2 comprehensive reporting, handoff export, and final collection
were enabled.

Result:

```text
Runtime: 22m16s
CPU hours: 5.0
Nextflow processes: 12/12 succeeded
Input records: 100
Downloaded/analyzed genomes: 99
Failed accession: GCF_055382775.1
Feature rows: 11,488
Unmatched feature rows: 0
Invalid feature rows: 0
Duplicate feature rows: 0
ABRicate ncbi/vfdb/plasmidfinder setup: PASS
PanR2 handoff report pages: generated
Runtime/resource summary: generated
```

Feature rows by table, excluding headers:

```text
amr.features.tsv: 1,234
integronfinder.features.tsv: 418
mlst.features.tsv: 1,622
plasmidfinder.features.tsv: 393
vfdb.features.tsv: 7,821
all_features.tsv: 11,488
```

Native-runner merge audit:

```text
abricate: PASS, expected_raw_tables=297, observed_raw_tables=297, samples_processed=99, feature_rows=9448
integronfinder: PASS, expected_raw_tables=99, observed_raw_tables=99, samples_processed=99, feature_rows=418
mlst: PASS, samples_processed=99, feature_rows=763
mobileelementfinder: SKIPPED, not requested
```

Runtime/resource observations:

```text
FETCHM: 7.78m, peak RSS 0.167 GiB, peak VMEM 3.3 GiB
PANR2_FEATURE_RUNNERS: 11.82m, peak RSS 2.2 GiB, peak VMEM 33.4 GiB
PANR2_COMPREHENSIVE: 1.95m, peak RSS 1.5 GiB, peak VMEM 5.0 GiB
```

This validates the main containerized comparative-genomics path at the
100-record scale using the local image. It does not replace full CheckM2,
AMRFinderPlus, ANI, or GTDB-Tk container validation.

## geNomad Docker Validation

The real geNomad database downloader was run inside Docker with a mounted
writable database directory:

```bash
sudo docker run --rm \
  -v /tmp/panresistome_genomad_db:/genomad_db \
  panresistome:experimental \
  bash -lc 'genomad download-database /genomad_db'
```

Result:

```text
geNomad database v1.9 downloaded from NERSC
Database extracted to /genomad_db/genomad_db
geNomad database ready
```

A two-genome geNomad-enabled Docker run then completed:

```text
Runtime: 5m03s
Nextflow processes: 16/16 succeeded
geNomad process: PASS
Feature rows: 286
prophage.features.tsv: header-only, 0 biological prophage rows
Unmatched/invalid/duplicate feature rows: 0/0/0
geNomad database setup/status: PASS
```

Interpretation: the geNomad Docker runner, database mount, PanR2 handoff, and
feature-contract path are valid. Positive geNomad biological feature calls still
need a prophage-rich validation dataset.

## GHCR Public Pull Status

An unauthenticated public pull was initially attempted:

```bash
sudo docker pull ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

The initial pull started successfully and downloaded layers without requiring a
GitHub login, so authentication was not the observed blocker. The pull was later
rerun to completion on 2026-05-13, and the pulled GHCR image completed a
two-genome geNomad-enabled Docker biological validation. The detailed remote
Docker validation result is recorded in
`validation/deployment/DOCKER_REMOTE_USER_VALIDATION_RESULTS.md`.

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
geNomad database download completes inside Docker with a mounted writable path
Two-genome geNomad-enabled Docker run completes with clean feature-contract validation
100-record Klebsiella Docker large-mode run completes with clean feature-contract validation
```

Not validated:

```text
Apptainer runtime execution
Containerized database mounts
Normal-user Nextflow Docker profile on this host
Full PanResistome production image pull from GHCR
ECTyper species-ID sketch readiness
Positive geNomad biological feature rows on a prophage-rich dataset
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
local image can run the Nextflow test profile, a two-genome biological
Klebsiella validation, a two-genome geNomad-enabled validation, and a
100-record large-mode Klebsiella validation through Docker. The GitHub Actions
workflow should still publish a fully pullable GHCR image before Docker,
Apptainer, or Singularity profiles are advertised as a low-hassle remote-user
route.

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
