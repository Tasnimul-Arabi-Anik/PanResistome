# PanResistome v0.4.0 Roadmap

Theme: large-dataset scalability and deployment groundwork.

## Current v0.4.0 groundwork

- `--large_dataset true` enables compact report defaults for large feature matrices.
- `--report_mode compact|publication|exploratory` controls handoff HTML density.
- `--max_features_heatmap`, `--max_features_network`, `--max_metadata_columns`, and `--top_n_features_per_database` limit report-facing summaries while complete TSV exports remain unchanged.
- `panr2_inputs/manifest/report_controls.tsv` records the applied report limits.
- `-profile large` enables large-dataset safeguards. Combine it with `desktop_parallel` or `workstation` when parallel native feature runners are desired.
- `docs/hpc.md` and `docs/containers.md` document deployment planning without claiming unvalidated container/HPC support.
- The existing 100-record `Klebsiella pneumoniae` validation input has passed with `-profile conda,mamba,desktop_parallel,large`; see `validation/klebsiella_pneumoniae_100/LARGE_MODE_VALIDATION_RESULTS.md`.
- AMRFinderPlus now supports bounded per-assembly parallelism through `--amrfinderplus_jobs` and `--amrfinderplus_threads_per_sample`.
- Completed runs write `pipeline_runtime_summary.tsv` and `pipeline_runtime_tasks.tsv` from the Nextflow trace.
- A 300-record BioProject-diverse `Klebsiella pneumoniae` large-mode validation passed with CheckM2, GTDB-Tk, ANI, and AMRFinderPlus disabled for desktop safety; see `validation/klebsiella_pneumoniae_300/LARGE_MODE_CHECKM2_OFF_VALIDATION_RESULTS.md`.
- The 300-record validation produced 299 downloaded genomes, 299 QC PASS calls, 36,638 standardized feature rows, zero unmatched/invalid/duplicate feature rows, and a compact large-mode report.
- `--ani_large_run_strategy auto|all|skip` and `--ani_max_all_vs_all_genomes` now protect large-dataset runs from accidental all-vs-all ANI launches while preserving ANI audit/status outputs.
- AMRFinderPlus resumed runs now reuse existing per-sample TSV outputs and print bounded progress updates.
- `docs/optional_module_validation_matrix.md` now separates stable default modules, stable table-input paths, experimental runners, and restricted database workflows.
- Synthetic PanR2 contract tests cover optional table-input exports for MobileElementFinder, ISfinder, MOB-suite, prophage/geNomad, DefenseFinder, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder.

## Near-term targets

1. Run targeted 2-10 genome smoke validations for optional runners before using them in large organism runs: MobileElementFinder, MOB-suite, geNomad, Kleborate/Kaptive, and ECTyper.
2. Add a faster full ANI strategy for 300+ genomes: chunked FastANI, skani, Mash-prescreened pairs, or representative-only ANI.
3. Treat AMRFinderPlus at 300+ genomes as a separate workstation/HPC benchmark unless nucleotide `tblastn` runtime is acceptable.
4. Further split or chunk native feature runners if ABRicate/IntegronFinder wall time becomes limiting above 300 genomes.
5. Validate large-dataset mode with a synthetic enlarged feature table if a faster regression fixture is needed.
6. Improve report navigation for compact/publication/exploratory modes.
7. Draft experimental Docker/Apptainer/SLURM profiles only after a container smoke test.

## Not v0.4.0 blockers

- New AMR databases such as CARD/RGI or ResFinder/PointFinder.
- Default MOB-suite/geNomad execution.
- Default organism-specific typing runners.
- GTDB-Tk default execution.

These should remain optional until dependency setup and validation are stable.
