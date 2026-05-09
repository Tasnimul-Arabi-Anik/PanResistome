# PanResistome v0.4.0 Roadmap

Theme: large-dataset scalability and deployment groundwork.

## Current v0.4.0 groundwork

- `--large_dataset true` enables compact report defaults for large feature matrices.
- `--report_mode compact|publication|exploratory` controls handoff HTML density.
- `--max_features_heatmap`, `--max_features_network`, `--max_metadata_columns`, and `--top_n_features_per_database` limit report-facing summaries while complete TSV exports remain unchanged.
- `panr2_inputs/manifest/report_controls.tsv` records the applied report limits.
- `-profile large` enables large-dataset safeguards. Combine it with `desktop_parallel` or `workstation` when parallel native feature runners are desired.
- `docs/hpc.md` and `docs/containers.md` document deployment planning without claiming unvalidated container/HPC support.

## Near-term targets

1. Validate large-dataset mode with a synthetic enlarged feature table and the existing Nextflow test profile.
2. Run a 300-500 genome validation with GTDB-Tk disabled, CheckM2 capped, parallel native runners, and `-profile large`.
3. Improve report navigation for compact/publication/exploratory modes.
4. Add runtime/resource summary tables from Nextflow trace output when available.
5. Draft experimental Docker/Apptainer/SLURM profiles only after a container smoke test.

## Not v0.4.0 blockers

- New AMR databases such as CARD/RGI or ResFinder/PointFinder.
- Default MOB-suite/geNomad execution.
- Default organism-specific typing runners.
- GTDB-Tk default execution.

These should remain optional until dependency setup and validation are stable.
