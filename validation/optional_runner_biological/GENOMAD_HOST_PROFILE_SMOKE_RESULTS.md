# geNomad Host-Environment Profile Smoke Test

Date: 2026-05-12

Purpose: verify that PanResistome can bypass Nextflow geNomad Conda environment creation when the user supplies a prebuilt host/module/container geNomad environment.

This addresses the bottleneck observed in `GENOMAD_AUTO_DOWNLOAD_ATTEMPT_RESULTS.md`, where the run reached `GENOMAD_PROPHAGE` but spent about 17 minutes creating `envs/genomad.yaml` before geNomad database download began.

## Implementation Tested

New option:

```text
--genomad_use_host_env true
```

New profile:

```text
-profile genomad_host
```

Effect:

```text
GENOMAD_ENV_VERSIONS and GENOMAD_PROPHAGE do not create envs/genomad.yaml.
They use the host/container PATH instead.
Other Conda-backed processes still use the selected Conda/Mamba profile.
```

## Smoke Command

This was a negative smoke test: no `genomad` executable and no geNomad database were available in the host PATH. The expected behavior was a fast auditable failure, not a Conda solve.

```bash
nextflow run main.nf \
  --local_samples validation_runs/broad_optional_small_input \
  --outdir validation_runs/genomad_host_missing_tool_smoke \
  -profile conda,mamba,genomad_host \
  --analysis_profile custom \
  --sequence_qc_engine python \
  --qc_filter false \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive false \
  --panr2_native_feature_runners false \
  --run_abricate false \
  --run_genomad true \
  --genomad_auto_download_db false \
  --threads 1 \
  --capture_versions false \
  --export_panr2_inputs false
```

## Result

Status: PASS for host-profile behavior

```text
Nextflow processes succeeded: 5/5
GENOMAD_PROPHAGE runtime: 0.10s
geNomad Conda environment creation attempted: no
geNomad database setup status: FAIL, expected
genomad_available: false
auto_download_requested: false
```

Audit output:

```text
prophage/genomad_database_setup_status.tsv
prophage/module_status.tsv
pipeline_runtime_summary.tsv
```

`module_status.tsv` recorded:

```text
module=genomad
enabled=true
status=FAIL
message=geNomad database setup failed; inspect prophage/genomad_database_setup_status.tsv.
```

## Interpretation

The `genomad_host` profile solves the immediate first-run environment bottleneck for users who already have geNomad available through a host module, prebuilt Conda environment, or container image.

It does not validate geNomad biological output by itself. A future positive validation still needs:

```text
genomad executable available in PATH
populated --genomad_db or successful genomad download-database
2-10 genome biological run
prophage.features.tsv generated
all_features.tsv includes prophage/geNomad rows
schema validation has zero unmatched/invalid/duplicate rows
```

Recommended next geNomad command on a prepared host/container environment:

```bash
python scripts/check_genomad_readiness.py \
  --db-dir /path/to/genomad_db \
  --out genomad_readiness.tsv

nextflow run main.nf \
  -profile conda,mamba,genomad_host \
  --local_samples validation_runs/broad_optional_small_input \
  --outdir validation_runs/genomad_host_positive_validation \
  --analysis_profile custom \
  --run_genomad true \
  --genomad_db /path/to/genomad_db \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive false \
  --run_abricate false \
  --export_panr2_inputs true \
  --threads 2
```
