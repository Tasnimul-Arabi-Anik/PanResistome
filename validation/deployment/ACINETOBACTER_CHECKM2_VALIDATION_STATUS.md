# Acinetobacter CheckM2 Comprehensive Validation Status

Date: 2026-05-16

## Status

Conda and Docker/GHCR CheckM2 runtime smokes passed on 2026-05-16. The Docker
Nextflow QC-only fixture also passed after the version-capture smoke was patched
to use the Python executable next to the resolved `checkm2` binary. A 5-genome
Acinetobacter Docker/GHCR comprehensive run then passed with CheckM2, QUAST,
ANI, Mash, AMRFinderPlus, geNomAD, ABRicate ncbi/vfdb/plasmidfinder,
IntegronFinder, MLST, and PanR2 comprehensive analysis enabled. GTDB-Tk,
DefenseFinder, MobileElementFinder, and ISfinder remained disabled by design.

The remaining validation gap is the separate opt-in MobileElementFinder track.

## Reproduced Failure

Local fixture command:

```bash
conda run -n checkm2 checkm2 predict \
  --threads 1 \
  --input /tmp/panresistome_checkm2_smoke/input \
  --output-directory /tmp/panresistome_checkm2_smoke/out \
  -x fna \
  --force \
  --lowmem \
  --database_path results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
```

Observed failure:

```text
Saved models could not be loaded: SavedModel file does not exist at:
.../checkm2/models/specific_model_COMP.keras/{saved_model.pbtxt|saved_model.pb}
```

The failing local environment had CheckM2 1.1.0 with Python 3.8, TensorFlow
2.4.1, Keras 2.4.3, NumPy 1.19.2, scikit-learn 0.23.2, and h5py 2.10.0.

## Package Route Selected

`envs/checkm2.yaml` now pins:

```text
checkm2=1.1.0=pyh7e72e81_1
python=3.12
tensorflow=2.17.*=cpu*
tensorflow-base=2.17.*=cpu*
keras=3.*
scikit-learn=1.6.1
diamond=2.1.11
```

A mamba dry-run solve for this CPU-only route completed successfully. A fresh
temporary install at `/tmp/panresistome_checkm2_env_completion` also completed
successfully and reported `checkm2 --version` as `1.1.0`.

Runtime model-load smoke:

```bash
conda run -p /tmp/panresistome_checkm2_env_completion \
  python -c "from checkm2 import modelProcessing; modelProcessing.modelProcessor(1); print('checkm2_model_load=PASS')"
```

Observed result:

```text
checkm2_model_load=PASS
```

Prediction smoke:

```bash
conda run -p /tmp/panresistome_checkm2_env_completion checkm2 predict \
  --threads 1 \
  --input /tmp/panresistome_checkm2_completion/input \
  --output-directory /tmp/panresistome_checkm2_completion/out \
  -x fna \
  --force \
  --lowmem \
  --database_path results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
```

Input genome:

```text
results_real_validation_checkm2_auto_v2/Klebsiella_oxytoca/sequence/GCF_001022115.1_ASM102211v1_genomic.fna
```

Observed result:

```text
CheckM2 finished successfully.
quality_report.tsv rows: 2
GCF_001022115.1_ASM102211v1_genomic completeness: 100.0
GCF_001022115.1_ASM102211v1_genomic contamination: 1.11
Completeness model: Neural Network (Specific Model)
```

No `specific_model_COMP.keras` or Keras model-loading error occurred in the
fresh environment.

An initial sandboxed `predict` attempt failed before annotation because Python
multiprocessing could not create its manager socket. The same command passed
when rerun outside the sandbox, so the failure was a harness restriction rather
than a CheckM2 package/runtime failure.

## Docker/GHCR Fixture Pass

Docker socket access was granted locally with an ACL, then the existing GHCR
image was validated directly:

```bash
docker run --rm ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  checkm2 --version
```

Observed result:

```text
1.1.0
```

Model-load smoke:

```bash
docker run --rm ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  /opt/conda/envs/checkm2_env/bin/python -c \
  "from checkm2 import modelProcessing; modelProcessing.modelProcessor(1); print('checkm2_model_load=PASS')"
```

Observed result:

```text
checkm2_model_load=PASS
```

Standalone Docker prediction smoke:

```bash
docker run --rm --user 1000:1000 \
  -v /tmp/panresistome_checkm2_completion/input:/input:ro \
  -v /tmp/panresistome_checkm2_docker_completion/out:/output \
  -v /home/anik/genomics/Tools_Dev/PanResistome/results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd:/db/uniref100.KO.1.dmnd:ro \
  ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  checkm2 predict \
    --threads 1 \
    --input /input \
    --output-directory /output \
    -x fna \
    --force \
    --lowmem \
    --database_path /db/uniref100.KO.1.dmnd
```

Observed result:

```text
CheckM2 finished successfully.
quality_report.tsv rows: 2
GCF_001022115.1_ASM102211v1_genomic completeness: 100.0
GCF_001022115.1_ASM102211v1_genomic contamination: 1.11
Completeness model: Neural Network (Specific Model)
```

Nextflow Docker QC-only fixture:

```bash
nextflow run main.nf \
  -profile docker \
  --local_samples /tmp/panresistome_docker_pipeline_local_samples \
  --outdir validation_runs/docker_checkm2_fixture_fixed \
  --run_checkm2 true \
  --checkm2_db results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
  --stop_after_qc true \
  --qc_filter true \
  --sequence_qc_engine python \
  --threads 2 \
  --checkm2_threads 1 \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

Observed result:

```text
Succeeded: 7
Duration: 7m 29s
pipeline_versions/checkm2_env_versions.txt includes checkm2_model_load=PASS
checkm2/quality_report.tsv rows: 2
sequence_qc/qc_decisions.tsv combined_qc_status: PASS
```

This Docker pipeline run initially exposed an integration issue: in the
all-in-one image, plain `python` resolves to the first environment on `PATH`,
not necessarily the CheckM2 environment. `CHECKM2_ENV_VERSIONS` now resolves
`command -v checkm2` and uses the sibling `python` executable for the model-load
smoke, which is compatible with both Conda and Docker profiles.

## Remaining Pipeline Fixture Pass

Before the Conda/Mamba CheckM2 route is called equally verified at the Nextflow
pipeline level, run one or two fixture genomes with CheckM2 enabled through the
Conda profile:

```bash
nextflow run main.nf \
  -profile conda,mamba \
  --local_samples tests/fixtures/local_samples \
  --outdir validation_runs/checkm2_fixture \
  --run_checkm2 true \
  --checkm2_db results_real_validation_checkm2_auto_v2/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
  --stop_after_qc true \
  --qc_filter true \
  --threads 2 \
  --checkm2_threads 1
```

Required outputs:

```text
checkm2/quality_report.tsv with real genome rows
metadata_output/ncbi_enriched.csv with CheckM2 columns
sequence_qc/qc_decisions.tsv
no specific_model_COMP.keras or Keras model-loading error
```

## Acinetobacter Stable Comprehensive Target

Observed PASS command on 2026-05-16:

```bash
nextflow run main.nf -profile docker \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 5 \
  --outdir validation_runs/acinetobacter_pittii_5_docker_comprehensive \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --genomad_jobs 1 \
  --genomad_threads_per_sample 1 \
  --genomad_splits 8 \
  --genomad_sensitivity 3.0 \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder false \
  --panr2_run_defensefinder false \
  --threads 4 \
  --checkm2_threads 1 \
  --fetchm2_download_workers 2
```

Observed run summary:

```text
Completed: 2026-05-16
Initial full-run duration: 2h 21m 30s
Final post-fix resume duration: 17m 44s
Final Nextflow status: Succeeded 12, Cached 11
Input records: 5
Downloaded/analyzed genomes: 5
CheckM2 database path: validation_runs/acinetobacter_pittii_5_docker_comprehensive/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd
CheckM2 database size: 2.9G
geNomAD database size: 1.4G
```

Required output evidence:

```text
checkm2/quality_report.tsv: 5 real genome rows
metadata_output/ncbi_enriched.csv: CheckM2 columns present
sequence_qc/qc_decisions.tsv: 5/5 combined_qc_status PASS
pipeline_versions/checkm2_env_versions.txt: checkm2_model_load=PASS
panr2_native_feature_runners/module_status.tsv:
  abricate PASS, integronfinder PASS, mlst PASS, mobileelementfinder SKIPPED
panr2_native_feature_runners/native_runner_merge_audit.tsv:
  abricate PASS, integronfinder PASS, mlst PASS, mobileelementfinder SKIPPED
prophage/module_status.tsv: geNomAD PASS, 5/5 samples processed, 27 feature rows
panr2_inputs/manifest/database_setup_status.tsv:
  required CheckM2, QUAST, ANI, Mash, AMRFinderPlus, geNomAD, ABRicate,
  IntegronFinder, MLST, and PanR2 checks PASS
panr2_inputs/manifest/schema_validation_summary.txt:
  feature_files_checked=7
  feature_rows=630
  unmatched_feature_rows=0
  invalid_feature_rows=0
  duplicate_feature_rows=0
panr2_inputs/features/all_features.tsv: 630 data rows
report/index.html: present
```

Feature rows by table:

```text
amr.features.tsv: 32
amrfinderplus.features.tsv: 28
integronfinder.features.tsv: 0
mlst.features.tsv: 81
plasmidfinder.features.tsv: 0
prophage.features.tsv: 27
vfdb.features.tsv: 462
all_features.tsv: 630
```

This run used no `--checkm2_db` argument; CheckM2 auto-downloaded and reused
its database under the output directory. No `specific_model_COMP.keras` or
Keras model-loading error occurred.

The first full pass exposed a downstream audit gap for zero-hit PlasmidFinder:
raw ABRicate PlasmidFinder files existed, but no header-only
`plasmidfinder.features.tsv` was exported. The final resume records the fix:
header-only native ABRicate tables are exported into PanR2 feature tables, and
PlasmidFinder now appears as `WARNING_EMPTY` instead of
`FAIL_MISSING_FEATURE_TABLE`.

Larger follow-up target:

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 10 \
  --outdir validation_runs/acinetobacter_pittii_10_docker \
  -profile docker \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder false \
  --panr2_run_defensefinder false \
  --threads 8 \
  --checkm2_threads 4 \
  --fetchm2_download_workers 2
```

## MobileElementFinder Follow-Up Target

After the stable comprehensive run passes, rerun or resume with:

```bash
--panr2_run_mobileelementfinder true \
--panr2_mobileelementfinder_allow_failure true
```

Accept either real MobileElementFinder rows with `PASS`, or
`WARNING_FAILED` with header-only PanR2-compatible output, a module status row,
a native runner merge-audit row, and no workflow abort.
