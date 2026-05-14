# PanResistome v0.5.0 Release Checklist

Release theme: reproducible, contract-governed, lineage-aware comparative
genomics.

## Release Gates

- [x] `pytest -q` passes.
- [x] `python -m py_compile scripts/*.py` passes.
- [x] `bash -n scripts/bootstrap.sh` passes.
- [x] `nextflow config -profile docker,large -flat` resolves.
- [x] `nextflow config -profile conda,mamba,large -flat` resolves.
- [x] `nextflow run main.nf -profile test` passes.
- [x] `git diff --check` passes.
- [x] Compact 10-genome Docker/GHCR biological validation is documented.
- [x] `reproducibility_manifest.json` is generated and documented.
- [x] `feature_contract.json` is generated and documented.
- [x] Lineage-aware outputs are generated and documented.
- [x] Pan-feature diversity outputs are generated and documented.
- [x] Statistical summary outputs are generated and documented.
- [x] `validation_runs/`, `work/`, downloaded FASTA files, downloaded
  databases, and container image artifacts are not tracked.
- [ ] CI passes for the v0.5.0 release commit after push.
- [ ] Annotated `v0.5.0` tag is created only after release-commit CI passes.

## Validation Evidence

- `validation/klebsiella_pneumoniae_10/V0_5_0_DOCKER_VALIDATION_RESULTS.md`
- `docs/validation_matrix.md`
- `docs/feature_contract_spec.md`
- `docs/panr2_input_contract.md`

## Notes

- The compact v0.5.0 Docker/GHCR validation intentionally disables CheckM2,
  QUAST, ANI, and GTDB-Tk. It validates the new v0.5.0 interpretation layer, not
  every heavyweight QC/taxonomy path.
- ANI-based lineage warnings were not exercised in the compact validation
  because ANI was disabled by design.
- Existing release tags must not be moved or retagged.
