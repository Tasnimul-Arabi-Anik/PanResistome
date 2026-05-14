# PanResistome v0.4.0 Release Checklist

Use this checklist before tagging v0.4.0. Do not move existing tags.

## Required Gates

- [x] CI passes for the release commit after push.
- [x] `pytest -q` passes locally.
- [x] `python -m py_compile scripts/*.py` passes.
- [x] `bash -n scripts/bootstrap.sh` passes.
- [x] `nextflow config -profile docker,large,genomad_host -flat` resolves.
- [x] `nextflow run main.nf -profile test` passes.
- [x] `git diff --check` passes.
- [x] Version fields are bumped to `0.4.0` in `VERSION`, `main.nf`, `nextflow.config`, and README.
- [x] `CHANGELOG.md` promotes current changes to `0.4.0 - 2026-05-13` and leaves `Unreleased` empty.
- [x] `docs/validation_matrix.md` includes the 300-genome large-mode validation and 5-genome geNomAD validation.
- [x] `validation/klebsiella_pneumoniae_300/LARGE_MODE_CHECKM2_OFF_VALIDATION_RESULTS.md` exists.
- [x] `validation/deployment/GHCR_DOCKER_100_VALIDATION_RESULTS.md` exists.
- [x] `validation/deployment/SINGULARITY_GHCR_VALIDATION_RESULTS.md` exists.
- [x] `validation/optional_runner_biological/GENOMAD_POSITIVE_CALL_VALIDATION_RESULTS.md` includes the 5-genome scale result.
- [x] Optional modules remain opt-in unless their default-path validation is complete.
- [x] ISfinder and GTDB-Tk still require explicit user-supplied database paths.
- [x] Complete TSV outputs remain preserved when `--large_dataset true` is enabled.
- [x] No generated `validation_runs/`, work directories, downloaded FASTA files, downloaded databases, Docker images, or SIF files are tracked.

## Release Decision

Tag v0.4.0 only after the release commit is pushed and CI passes. If CI fails, fix only the failing release gate and retest before tagging.

Final release status: `v0.4.0` was tagged at commit `20cbb68c0efa99a0f4fe07bf8367a998cc4ede21`. The release commit CI passed, and the tag-triggered container workflow completed the GHCR image build/push and published-image smoke test successfully: <https://github.com/Tasnimul-Arabi-Anik/PanResistome/actions/runs/25817372496>.
