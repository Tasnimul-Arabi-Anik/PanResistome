# Docker Remote-User Validation Results

Date: 2026-05-12 to 2026-05-13

Purpose: validate whether the container route can support remote-user style
PanResistome execution without Conda/Mamba solving, and identify the remaining
deployment blockers honestly.

## Summary

The Docker route is now biologically functional with both the local
`panresistome:experimental` image and the pulled GHCR image:

```text
Docker test profile: PASS
Two-genome Klebsiella biological Docker run: PASS
geNomad database download inside Docker: PASS
Two-genome geNomad Docker biological run: PASS
100-record Klebsiella Docker large-mode run: PASS
GHCR unauthenticated pull starts: PASS
GHCR full pull on this machine: PASS
Two-genome GHCR geNomad Docker biological run: PASS
Normal non-sudo Docker on this machine: NOT READY, user is not in docker group
```

## Normal Non-Sudo Docker Access

Current user groups:

```text
anik nogroup
```

Normal Docker commands fail:

```text
permission denied while trying to connect to the docker API at unix:///var/run/docker.sock
```

Interpretation: the pipeline works through Docker, but this host still requires
`sudo` for Docker. A normal user Docker route requires adding the user to the
`docker` group and opening a new login session. That change grants broad
Docker/root-equivalent access and should be made explicitly by the machine
owner or administrator.

## GHCR Pull Status

Command:

```bash
sudo docker pull ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

Result:

```text
Authentication/visibility: PASS, pull started without GitHub login.
Full pull: PASS on this machine.
Image ID: sha256:66da115ec47fbb91d56a0773ebe054c5d3aac0c3aaa200a1b05e228bf655c7ca
Image size: 7.45 GB
```

Interpretation: GHCR package visibility and image pull are validated for the
experimental image. Remote users should not need to log in if the GHCR package
remains public, but the image is large enough that first pull can still be slow
on weak networks.

## Two-Genome geNomad GHCR Docker Biological Run

Command:

```bash
sudo nextflow run main.nf \
  --input /tmp/klebsiella_2_container_ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_2_ghcr_genomad \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --container_run_options '-v /tmp/panresistome_genomad_db:/tmp/panresistome_genomad_db' \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus false \
  --run_genomad true \
  --genomad_use_host_env true \
  --genomad_db /tmp/panresistome_genomad_db/genomad_db \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 4 \
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_klebsiella_2_ghcr_genomad_work
```

Result:

```text
Runtime: 7m21s
CPU hours: 0.6
Nextflow processes: 16/16 succeeded
Downloaded genomes: 2
QC PASS: 2
Feature rows: 286
Feature tables checked: 5
Databases with standardized feature rows: amr, mlst, plasmidfinder, vfdb
Unmatched feature rows: 0
Invalid feature rows: 0
Duplicate feature rows: 0
ABRicate database setup/status: PASS
geNomad database setup/status: PASS
PanR2 handoff report: generated
Runtime/resource summary: generated
```

Runtime/resource observations:

```text
GENOMAD_PROPHAGE: 3.90m, peak RSS 2.7 GiB, peak VMEM 23.5 GiB
PANR2_FEATURE_RUNNERS: 2.02m, peak RSS 0.62 GiB, peak VMEM 10.9 GiB
PANR2_COMPREHENSIVE: 40.80s, peak RSS 0.416 GiB, peak VMEM 3.9 GiB
```

Feature table result:

```text
amr.features.tsv: 9 rows
mlst.features.tsv: 34 rows
plasmidfinder.features.tsv: 3 rows
prophage.features.tsv: header-only, 0 biological prophage rows
vfdb.features.tsv: 240 rows
all_features.tsv: 286 rows
```

Interpretation: the pulled GHCR image can execute the real Nextflow Docker
profile, use a mounted geNomad database, run geNomad, run native feature
runners, export PanR2-compatible standardized feature tables, and generate clean
schema validation. As with the local-image geNomad run, these two Klebsiella
genomes did not yield standardized prophage feature rows, so positive geNomad
feature validation still needs a prophage-rich dataset.

## geNomad Docker Database Download

Command:

```bash
sudo docker run --rm \
  -v /tmp/panresistome_genomad_db:/genomad_db \
  panresistome:experimental \
  bash -lc 'genomad download-database /genomad_db'
```

Result:

```text
Requested: https://portal.nersc.gov/genomad/__data__/genomad_db_v1.9.tar.gz
Downloaded: PASS
Extracted: /genomad_db/genomad_db
Database: geNomad v1.9 ready
```

This validates that the image can run geNomad's real database downloader when
a writable host directory is mounted.

## Two-Genome geNomad Docker Biological Run

Command:

```bash
sudo nextflow run main.nf \
  --input /tmp/klebsiella_2_container_ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_2_docker_genomad \
  -profile docker,large \
  --container_image panresistome:experimental \
  --container_run_options '-v /tmp/panresistome_genomad_db:/tmp/panresistome_genomad_db' \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus false \
  --run_genomad true \
  --genomad_use_host_env true \
  --genomad_db /tmp/panresistome_genomad_db/genomad_db \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 4 \
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_klebsiella_2_docker_genomad_work
```

Result:

```text
Runtime: 5m03s
Nextflow processes: 16/16 succeeded
geNomad process: PASS
Runtime/resource summary: PASS
Downloaded genomes: 2
QC PASS: 2
Feature rows: 286
Feature tables checked: 5
Unmatched feature rows: 0
Invalid feature rows: 0
Duplicate feature rows: 0
geNomad database setup/status: PASS
```

Feature table result:

```text
prophage.features.tsv: header-only, 0 biological prophage rows
```

Interpretation: the geNomad Docker runner, database mount, PanR2 handoff, and
feature-contract path are valid. These two Klebsiella genomes did not yield
standardized prophage feature rows, so a larger or prophage-rich dataset is
still needed to validate positive geNomad biological calls.

## 100-Record Klebsiella Docker Large-Mode Run

Command:

```bash
sudo nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_100_docker_large \
  -profile docker,large \
  --container_image panresistome:experimental \
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
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_klebsiella_100_docker_large_work
```

Result:

```text
Runtime: 22m16s
CPU hours: 5.0
Nextflow processes: 12/12 succeeded
Input records: 100
Downloaded/analyzed genomes: 99
Failed accession: GCF_055382775.1
Feature rows: 11,488
Feature tables checked: 5
Unmatched feature rows: 0
Invalid feature rows: 0
Duplicate feature rows: 0
HTML handoff report pages: generated
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

Database setup status:

```text
ABRicate ncbi: PASS
ABRicate vfdb: PASS
ABRicate plasmidfinder: PASS
IntegronFinder: PASS
MLST: PASS
Mash: PASS
GTDB-Tk: SKIPPED, disabled by default
CheckM2: SKIPPED, disabled for this desktop-safe Docker validation
AMRFinderPlus: SKIPPED, disabled for this desktop-safe Docker validation
```

Runtime/resource observations:

```text
FETCHM: 7.78m, peak RSS 0.167 GiB, peak VMEM 3.3 GiB
PANR2_FEATURE_RUNNERS: 11.82m, peak RSS 2.2 GiB, peak VMEM 33.4 GiB
PANR2_COMPREHENSIVE: 1.95m, peak RSS 1.5 GiB, peak VMEM 5.0 GiB
```

Interpretation: the main containerized comparative-genomics path is validated
at the 100-record Klebsiella scale using the local image. This run intentionally
disabled CheckM2, QUAST, ANI, and AMRFinderPlus to keep the validation focused
on Docker execution, FetchM2, sequence QC, Mash, native feature runners,
large-mode reporting, PanR2 handoff, and feature-contract cleanliness.

## Current Remote-User Judgment

Docker now substantially reduces user setup burden:

```text
No Conda/Mamba environment solving by the user
No manual ABRicate ncbi/vfdb/plasmidfinder setup in normal comprehensive mode
PanR2 reporting works inside the image
MOB-suite, geNomad, Kleborate, Kaptive, ECTyper, and MobileElementFinder entry points start in the image
```

Remaining caveats:

```text
First GHCR pull is large and may be slow.
Normal Docker use requires Docker socket permission or sudo.
GTDB-Tk remains external/opt-in because of database size.
ISfinder remains user-supplied because of licensing.
geNomad positive biological feature calls still need validation on a prophage-rich dataset.
AMRFinderPlus and CheckM2 were not included in the 100-genome Docker scale run.
```

Recommended next deployment target:

```text
1. Complete full GHCR pull on a faster network.
2. Run the same two-genome and 100-genome commands using the GHCR image name instead of the local image.
3. Validate Apptainer/Singularity with the GHCR image on an HPC-like host.
4. Add a prophage-rich 5-10 genome geNomad-positive validation.
```
