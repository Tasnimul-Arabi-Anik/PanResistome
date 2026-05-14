# PanResistome v0.5.0 Roadmap

Theme: reproducibility, schema governance, and deployment hardening.

## Primary targets

1. Keep v0.4.x stable with patch-only fixes unless a real release issue appears.
2. Treat `reproducibility_manifest.json` and `feature_contract.json` as the
   next stable audit layer for every PanR2 handoff export.
3. Validate Apptainer when a host with Apptainer is available.
4. Validate the experimental `slurm` profile on a real cluster before marking it
   stable.
5. Run a bounded 100-genome Docker/GHCR comprehensive validation with CheckM2
   and AMRFinderPlus enabled if memory and runtime allow.

## Not immediate targets

- New species-specific default databases.
- New default AMR databases.
- Default GTDB-Tk execution.
- 300-genome full comprehensive validation on a desktop.

## Release gates for v0.5.0

- `pytest -q`, `py_compile`, `bash -n`, `nextflow config`, and `-profile test`
  pass.
- `reproducibility_manifest.json` is produced in fixture and biological runs.
- `feature_contract.json` validates against the documented contract.
- Container image digest guidance is documented.
- Any new deployment profile is clearly labeled validated or experimental.
