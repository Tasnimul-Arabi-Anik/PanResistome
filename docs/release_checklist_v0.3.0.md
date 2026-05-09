# PanResistome v0.3.0 Release Checklist

Use this checklist before tagging v0.3.0. Do not move an existing tag; create a new tag only after these gates are satisfied.

## Required Gates

- [x] CI passes for the release commit.
- [x] `pytest` passes locally or in CI.
- [x] `nextflow run main.nf -profile test` passes.
- [x] `python -m py_compile` passes for pipeline helper scripts.
- [x] `bash -n scripts/bootstrap.sh` passes.
- [x] `git diff --check` passes.
- [x] `validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md` exists.
- [x] `validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md` exists.
- [x] `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md` exists.
- [x] `validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md` exists.
- [x] `docs/validation_matrix.md` is updated.
- [x] `validation_runs/`, `work/`, downloaded FASTA files, and downloaded databases are not committed.
- [x] CheckM2 thread safety is documented, including `--checkm2_threads 2` for 16 GB desktop runs.
- [x] The standard comprehensive command still works with GTDB-Tk off.
- [x] Required ABRicate databases (`ncbi`, `vfdb`, `plasmidfinder`) are setup-checked before PanR2 reporting.
- [x] AMRFinderPlus auto-update behavior is documented.
- [x] PanR2 feature-contract validation reports zero unmatched, invalid, and duplicate rows in the committed validation reports.
- [x] PanR2 dashboards/reports are generated in the committed validation reports.
- [x] `panr2_inputs/manifest/database_setup_status.tsv` is documented as the database/tool audit.
- [x] `panr2_inputs/manifest/native_runner_merge_audit.tsv` is documented as the native-runner merge audit.
- [x] README and roadmap state that PanResistome runs heavy tools and PanR2 analyzes standardized outputs.
- [x] Optional or restricted modules remain clearly marked opt-in.

## Release Decision

Tag v0.3.0 only if the required gates pass and no release-blocking parser, install, or report-generation issue is open. If only documentation or small parser fixes remain, keep them as v0.3.x patch candidates after v0.3.0.
