# PanResistome v0.6.0 Release Checklist

Release theme: deployment and broader biological validation.

## Release Gates

- [x] `pytest -q` passes.
- [x] `python -m py_compile scripts/*.py` passes.
- [x] `bash -n scripts/bootstrap.sh` passes.
- [x] `nextflow config -profile singularity,large -flat` resolves.
- [x] `nextflow config -profile conda,mamba,large -flat` resolves.
- [x] `nextflow run main.nf -profile test` passes.
- [x] `git diff --check` passes.
- [x] ANI-enabled 10-genome Singularity/GHCR validation is documented.
- [x] Positive IntegronFinder biological evidence is documented.
- [x] Native IntegronFinder/MobileElementFinder handoff export regression test is present.
- [x] `validation_runs/`, `work/`, downloaded FASTA files, downloaded
  databases, and container image artifacts are not tracked.
- [ ] CI passes for the v0.6.0 release commit after push.
- [ ] Annotated `v0.6.0` tag is created only after release-commit CI passes.

## Validation Evidence

- `validation/klebsiella_pneumoniae_10/ANI_LINEAGE_SINGULARITY_VALIDATION_RESULTS.md`
- `validation/integronfinder_positive/VALIDATION_RESULTS.md`
- `docs/validation_matrix.md`
- `docs/roadmap_v0.6.0.md`

## Notes

- The ANI-lineage validation intentionally disables CheckM2, QUAST, GTDB-Tk,
  and AMRFinderPlus. It validates Singularity/GHCR deployment, skani ANI,
  ANI-cluster lineage context, native feature runners, and strict PanR2 handoff
  export.
- The first uncached GHCR-to-SIF conversion can exceed Nextflow's default
  Singularity pull timeout on slow networks. `docs/containers.md` documents the
  validated `singularity.pullTimeout = '2h'` workaround.
- Existing release tags must not be moved or retagged.
