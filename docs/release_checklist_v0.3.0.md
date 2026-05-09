# PanResistome v0.3.0 Release Checklist

Use this checklist before tagging v0.3.0. Do not move an existing tag; create a new tag only after these gates are satisfied.

## Required Gates

- [ ] CI passes for the release commit.
- [ ] `pytest` passes locally or in CI.
- [ ] `nextflow run main.nf -profile test` passes.
- [ ] `python -m py_compile` passes for pipeline helper scripts.
- [ ] `bash -n scripts/bootstrap.sh` passes.
- [ ] `git diff --check` passes.
- [ ] `validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md` exists.
- [ ] `validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md` exists.
- [ ] `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md` exists.
- [ ] `validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md` exists.
- [ ] `docs/validation_matrix.md` is updated.
- [ ] `validation_runs/`, `work/`, downloaded FASTA files, and downloaded databases are not committed.
- [ ] CheckM2 thread safety is documented, including `--checkm2_threads 2` for 16 GB desktop runs.
- [ ] The standard comprehensive command still works with GTDB-Tk off.
- [ ] Required ABRicate databases (`ncbi`, `vfdb`, `plasmidfinder`) are setup-checked before PanR2 reporting.
- [ ] AMRFinderPlus auto-update behavior is documented.
- [ ] PanR2 feature-contract validation reports zero unmatched, invalid, and duplicate rows in the committed validation reports.
- [ ] PanR2 dashboards/reports are generated in the committed validation reports.
- [ ] `panr2_inputs/manifest/database_setup_status.tsv` is documented as the database/tool audit.
- [ ] `panr2_inputs/manifest/native_runner_merge_audit.tsv` is documented as the native-runner merge audit.
- [ ] README and roadmap state that PanResistome runs heavy tools and PanR2 analyzes standardized outputs.
- [ ] Optional or restricted modules remain clearly marked opt-in.

## Release Decision

Tag v0.3.0 only if the required gates pass and no release-blocking parser, install, or report-generation issue is open. If only documentation or small parser fixes remain, keep them as v0.3.x patch candidates after v0.3.0.
