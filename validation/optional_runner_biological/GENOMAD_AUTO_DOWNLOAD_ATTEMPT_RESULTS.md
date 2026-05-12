# geNomad Auto-Download Validation Attempt

Date: 2026-05-12

Purpose: test whether a remote-style user can enable geNomad with automatic database download on a small biological subset without manually preparing the geNomad database.

Input:

```text
Source validation input: validation_runs/broad_optional_small_input
Subset: 5 Klebsiella pneumoniae assemblies
```

Command attempted:

```bash
nextflow run main.nf \
  --local_samples validation_runs/broad_optional_small_input \
  --outdir validation_runs/genomad_auto_small_nextflow \
  -profile conda,mamba,lowmem \
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
  --run_genomad true \
  --genomad_auto_download_db true \
  --threads 2 \
  --capture_versions false \
  --export_panr2_inputs true \
  -resume
```

## Result

Status: DID NOT COMPLETE

The workflow reached the correct `GENOMAD_PROPHAGE` process, but it did not reach geNomad database download or genome analysis. The run was stopped after approximately 17 minutes because it was still creating the first geNomad Conda environment:

```text
mamba env create --prefix work/conda/env-305d9c7cec00ca07-... --file envs/genomad.yaml
```

Completed processes before stop:

```text
SEQUENCE_QC: cached/completed
COMBINED_QC: completed
GENOMAD_PROPHAGE: waiting on geNomad environment creation
```

Runtime summary contained only QC tasks:

```text
SEQUENCE_QC: cached
COMBINED_QC: completed
```

No geNomad database, raw geNomad output, prophage feature table, or PanR2 prophage export was generated.

## Interpretation

This is a real remote-user bottleneck. The geNomad auto-download path cannot be called validated yet because the first-run Conda/Mamba environment creation did not finish in a practical short validation window. The pipeline wiring is correct enough to reach `GENOMAD_PROPHAGE`, but the current user experience is not yet low-hassle for fresh machines.

Current honest support level:

```text
PanR2 prophage/geNomad table analysis: PASS
geNomad missing-DB audit helper: PASS
geNomad auto-download biological runner validation: NOT PASSED
Primary blocker observed: first-run geNomad Conda environment creation time
```

## Recommendation

Keep geNomad opt-in and do not advertise it as a hassle-free automated runner yet.

Recommended next engineering options:

1. Use `python scripts/check_genomad_readiness.py` before launching geNomad runs.
2. Use `-profile genomad_host` or `--genomad_use_host_env true` when geNomad is available from a host module, prebuilt Conda environment, or container image.
3. Prefer container/Apptainer validation for geNomad before claiming remote-user ease.
4. If Conda remains the supported route, test geNomad environment creation separately and document expected solve/install time.
5. Re-run biological validation only after the geNomad environment is cached or containerized, then test DB auto-download and 2-10 genomes.
6. Keep `--genomad_db` support for users/HPC systems with prebuilt database paths.

A future passing validation should document all of these:

```text
geNomad environment creation time
geNomad database download time and size
raw geNomad outputs per genome
prophage.features.tsv rows
all_features.tsv includes prophage/genomad rows
schema validation has zero unmatched/invalid/duplicate rows
report pages generated
```
