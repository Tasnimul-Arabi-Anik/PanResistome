# Singularity GHCR Validation Results

Date: 2026-05-13

Purpose: validate the public GHCR PanResistome image through Singularity CE,
including first pull/conversion behavior, a small geNomad-enabled biological
run, and a 100-record large-mode biological run.

## Runtime

```text
Singularity: singularity-ce version 4.1.1
Image: docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental
Apptainer: not installed on this validation host
```

## Pull And Conversion

Command:

```bash
singularity exec docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental true
```

Result:

```text
Status: PASS
Approximate first GHCR-to-SIF conversion time: 1h15m
Observed Singularity cache growth before extraction/SIF creation: about 7.3 GB
Cached SIF size used by Nextflow: about 4.3 GB
```

The first readiness-helper attempt timed out at the original 300 second pull
test limit. The helper now supports `--pull-test-timeout` so users can allow
longer first-time Singularity/Apptainer image conversion:

```bash
python scripts/check_container_readiness.py \
  --runtime singularity \
  --image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --pull-test \
  --pull-test-timeout 7200 \
  --out container_readiness.tsv
```

After the image was cached, the readiness helper passed:

```text
status=PASS
singularity is available: singularity-ce version 4.1.1; container pull/exec test passed.
```

## Singularity Profile Defaults

The first biological Singularity attempt showed an important container-specific
issue: ABRicate databases were present in the image, but the default forced
refresh tried to update them inside the read-only SIF filesystem.

The `apptainer` and `singularity` profiles now default to:

```text
params.panr2_update_abricate_db = false
env.MPLCONFIGDIR = /tmp
```

This keeps in-image ABRicate databases frozen/read-only by default and avoids
matplotlib config writes under a read-only home path.

## Two-Genome geNomad Biological Validation

Command:

```bash
nextflow run main.nf \
  --input /tmp/klebsiella_2_container_ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_2_singularity_genomad_profile_default \
  -profile singularity,large \
  --container_image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --container_run_options '-B /tmp/panresistome_genomad_db:/tmp/panresistome_genomad_db' \
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
  -w /tmp/panresistome_klebsiella_2_singularity_genomad_profile_default_work
```

Result:

```text
Status: PASS
Runtime after image availability: 2m08s
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

As with the Docker geNomad validation, `prophage.features.tsv` was header-only
for these two Klebsiella genomes. Positive geNomad feature calls still need a
prophage-rich validation dataset.

## 100-Record Large-Mode Biological Validation

Command:

```bash
env NXF_SINGULARITY_CACHEDIR=/tmp/panresistome_klebsiella_2_singularity_genomad_profile_default_work/singularity \
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir /tmp/panresistome_klebsiella_100_singularity_large \
  -profile singularity,large \
  --container_image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
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
  -w /tmp/panresistome_klebsiella_100_singularity_large_work
```

Result:

```text
Status: PASS
Runtime: 21m12s
CPU hours: 4.8
Nextflow processes: 12/12 succeeded
Input records: 100
Downloaded/analyzed genomes: 99
Failed accession: GCF_055382775.1
QC PASS: 99
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

Runtime/resource observations:

```text
FETCHM: 7.62m, peak RSS 0.184 GiB, peak VMEM 3.3 GiB
PANR2_FEATURE_RUNNERS: 10.93m, peak RSS 2.2 GiB, peak VMEM 33.4 GiB
PANR2_COMPREHENSIVE: 1.93m, peak RSS 1.5 GiB, peak VMEM 4.9 GiB
```

## Interpretation

Singularity CE is now validated with the public GHCR image for both:

```text
small geNomad-enabled biological workflow
100-record large-mode biological workflow
```

The main practical caveats are:

```text
First GHCR-to-SIF conversion is slow for the large all-in-one image.
Use NXF_SINGULARITY_CACHEDIR on shared/persistent storage.
ABRicate database force-refresh should remain disabled under read-only SIF execution.
Apptainer still needs validation on a host where it is installed.
```
