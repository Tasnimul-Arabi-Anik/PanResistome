# PanResistome v0.6.0 Roadmap

Theme: deployment and broader biological validation.

v0.5.0 established reproducibility manifests, feature-contract governance,
lineage-aware interpretation, diversity summaries, statistical summaries, and a
compact Docker/GHCR biological validation. v0.6.0 should focus on validating the
same architecture under additional deployment and biological evidence paths,
without adding new default databases.

## Release goals

- Add at least one new validation tier beyond v0.5.0:
  - Apptainer/Singularity validation, or
  - SLURM validation, or
  - ANI-enabled lineage validation, or
  - positive IntegronFinder biological validation.
- Keep existing v0.4.0 and v0.5.0 commands backward-compatible.
- Keep GTDB-Tk, ISfinder, geNomAD, MOB-suite, MobileElementFinder,
  DefenseFinder, and organism-specific modules opt-in.
- Do not add new default databases.

## Validation targets

### ANI-enabled lineage validation

Run a compact 10-genome `Klebsiella pneumoniae` workflow with ANI enabled and
CheckM2/GTDB-Tk disabled. The goal is to confirm that ANI-derived lineage
context can be written into PanR2 lineage outputs when pairwise ANI is available.

Expected outputs:

- `ani/analysis/ani_run_status.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_summary.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_adjusted_warnings.tsv`
- `panr2_inputs/report/lineage_context.html`

Status: completed on 2026-05-14 with Singularity/GHCR. See
`validation/klebsiella_pneumoniae_10/ANI_LINEAGE_SINGULARITY_VALIDATION_RESULTS.md`.

### Positive IntegronFinder validation

Use `validation/integronfinder_positive/ncbi_dataset.tsv`, selected from
previous real positive IntegronFinder calls, to confirm that a compact rerun can
produce biological IntegronFinder features quickly when unrelated heavy stages
are disabled.

Expected outputs:

- `panr2_inputs/features/integronfinder.features.tsv` with biological rows.
- `panr2_inputs/manifest/native_runner_merge_audit.tsv` with
  `integronfinder` status `PASS`.
- `panr2_inputs/manifest/schema_validation_summary.txt` with zero unmatched,
  invalid, and duplicate feature rows.

Status: completed as real positive-call evidence from the 100-record
Klebsiella GHCR/Docker validation and re-exercised in the 10-genome
ANI-lineage Singularity/GHCR validation. See
`validation/integronfinder_positive/VALIDATION_RESULTS.md`.

### Deployment validation

If available locally, prefer non-root Singularity or Apptainer validation over
sudo Docker. Docker remains validated through GHCR, but local non-root Docker use
depends on host group permissions.

## Non-goals

- No new default databases.
- No species-specific default modules.
- No full 300-genome comprehensive run as a blocker.
- No GTDB-Tk default execution.
- No unvalidated claim that SLURM or Apptainer is stable on all clusters.
