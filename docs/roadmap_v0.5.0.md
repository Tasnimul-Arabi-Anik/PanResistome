# PanResistome v0.5.0 Roadmap

Theme: reproducible, contract-governed, lineage-aware comparative genomics.

## Primary targets

1. Keep v0.4.x stable with patch-only fixes unless a real release issue appears.
2. Treat `reproducibility_manifest.json` and `feature_contract.json` as the
   next stable audit layer for every PanR2 handoff export.
3. Add lineage-aware interpretation so metadata-feature findings can be flagged
   when they are dominated by one MLST ST, one ANI cluster, or one BioProject.
4. Add pan-feature diversity summaries, including feature richness, database
   diversity, Jaccard distance, core/accessory/rare classification, and
   pan-feature accumulation.
5. Add a statistical-summary page that audits how many metadata-feature tests,
   top findings, q-value hits, and warning categories were produced.
6. Validate Apptainer when a host with Apptainer is available.
7. Validate the experimental `slurm` profile on a real cluster before marking it
   stable.
8. Run a compact Docker/GHCR biological validation that exercises the
   reproducibility manifest, feature contract manifest, lineage warnings,
   diversity summaries, statistical summary, and report pages.

## Not immediate targets

- New species-specific default databases.
- New default AMR databases.
- Default GTDB-Tk execution.
- 300-genome full comprehensive validation on a desktop.
- Complex phylogeny-aware regression or mixed models before simpler lineage
  warning layers are validated.

## Release gates for v0.5.0

- `pytest -q`, `py_compile`, `bash -n`, `nextflow config`, and `-profile test`
  pass.
- `reproducibility_manifest.json` is produced in fixture and biological runs.
- `feature_contract.json` validates against the documented contract.
- Lineage-aware, diversity, and statistical-summary outputs are produced in the
  test profile and one compact biological run.
- Container image digest guidance is documented.
- Any new deployment profile is clearly labeled validated or experimental.
