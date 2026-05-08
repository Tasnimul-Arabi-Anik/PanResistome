# Changelog

## Unreleased

### Added
- Optional NCBI AMRFinderPlus runner (`--run_amrfinderplus`) and precomputed table input (`--amrfinderplus_dir`) with version capture and standardized PanR2 contract export.
- Analysis profile presets through `--analysis_profile` for `qc_only`, `amr_basic`, `amr_vp`, `amr_vp_mge`, and `comprehensive` runs.
- Strict `panr2_inputs/features/*.features.tsv` contract exports plus schema validation reports, unmatched-feature reports, and all-feature merged tables.
- Optional standardized PanR2 feature columns for feature names/descriptions, mechanism, drug class, source table/file, raw IDs, evidence type, confidence, and notes.
- Feature completeness audit, module status summary, invalid/duplicate feature reports, metadata audit, metadata eligibility, feature eligibility, prevalence tables, database burden summaries, feature matrices, and cross-database context outputs with exact Fisher screening p-values and BH-FDR q-values in the PanR2 handoff bundle.
- Database/tool setup audit at `panr2_inputs/manifest/database_setup_status.tsv`, with strict required-database failure behavior for comprehensive PanR2 runs.
- MLST results are now standardized into `panr2_inputs/features/mlst.features.tsv`, including sequence-type and allele-level feature rows.
- Optional PanResistome ISfinder-compatible BLAST module (`--run_isfinder`, `--isfinder_db_fasta`) that builds a local BLAST database from an authorized ISfinder FASTA and exports PanR2-readable `isfinder/tables/*_results.tab` files.
- Reproducible current NCBI Assembly validation input for `Delftia tsuruhatensis` under `validation/delftia_tsuruhatensis_current/`.
- `scripts/generate_ncbi_assembly_input.py` for creating FetchM2/PanResistome-compatible validation inputs from NCBI Assembly E-utilities.
- README module stability table clarifying stable, active-development, and experimental modules for public users.
- `docs/remote_user_validation.md` documenting the fresh-clone standard comprehensive validation path and release-passing criteria.
- `pytest.ini` limiting test discovery to repository tests so local runs do not collect Nextflow work-directory Conda package tests.

### Changed
- Optional downstream annotation modules now run after the combined QC decision step, so `--qc_filter true` can use `sequence_filtered/` before AMRFinderPlus, MOB-suite, geNomad, and organism-specific typing.
- AMRFinderPlus database setup now runs by default before AMRFinderPlus execution, and the process records per-sample status instead of silently swallowing missing-database failures.
- CheckM2 QC now declares the same CPU count as `--threads` and allows longer wall time for full 45-genome laptop validation runs.
- MobileElementFinder is now opt-in through `--panr2_run_mobileelementfinder true`; real Delftia validation showed that the upstream `mefinder` BLAST JSON parser can fail on otherwise valid assemblies, so the public comprehensive default avoids a brittle hard failure.
- The FetchM2/QC Conda environment no longer installs PanR2 from GitHub; PanR2 is kept in the dedicated PanR2 comprehensive environment and the legacy `PANR` process now uses that environment.

### Fixed
- Relative `--checkm2_db` values are now resolved from the launch directory before CheckM2 runs in a Nextflow work directory.
- Comprehensive PanR2 profiles no longer request the ABRicate `isfinder` database by default because standard ABRicate setup does not always provide it; users can still opt in with `--panr2_abricate_dbs ...isfinder` when that database is installed.
- `--qc_filter true` now fails at the combined QC decision step with a clear report when zero FASTA files remain for downstream analysis, instead of failing later inside PanR2.

### Validated
- A 45-genome `Delftia tsuruhatensis` run completed with FetchM2, sequence QC, CheckM2, QUAST, ANI, Mash, combined QC, AMRFinderPlus database auto-download, AMRFinderPlus, comprehensive PanR2, and PanR2 handoff export, with GTDB-Tk disabled and 4 threads.
- The validated PanR2 handoff contained AMR, AMRFinderPlus, VFDB, PlasmidFinder, IntegronFinder, and MLST feature tables in `panr2_inputs/features/all_features.tsv` with zero unmatched, invalid, or duplicate feature rows.

## 0.2.1 - 2026-05-06

### Added
- FetchM2 dependency pin updated to `fetchm2==0.1.7`, including downstream contract files such as `sample_map.csv`, `metadata_completeness.csv`, `metadata_bias_warning.txt`, and `fetchm2_manifest.json`.
- Comprehensive PanR2 mode now pins PanR2 release `v0.1.3`, adding PanR2 analysis support for MOB-suite, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder table inputs.
- Automatic CheckM2 database download support through `--checkm2_auto_download_db` and `--checkm2_db_dir` when `--checkm2_db` is not supplied.
- Optional MOB-suite plasmid reconstruction/typing process with PanR2 handoff under `mobsuite/`.
- Optional geNomad prophage/viral-region process with PanR2 handoff under `prophage/` when a geNomad database is provided.
- Optional organism-specific typing process for Kleborate, Kaptive, and ECTyper, plus SerotypeFinder/SCCmecFinder table passthrough to PanR2.
- Explicit precomputed table passthrough parameters for DefenseFinder, MOB-suite, prophage/geNomad-style tables, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder so PanR2 can analyze these feature families without forcing every external runner into the default environment.
- `--panr2_sample_map` passthrough plus automatic use of FetchM2 `metadata_output/sample_map.csv` when available for matching external tool sample IDs to assembly accessions.
- A `mamba` Nextflow profile and bootstrap guidance for faster optional heavy-tool environment resolution.
- FetchM2 representative GCF-preferred `fetchm2_clean.csv` behavior is used by default while preserving full uncollapsed rows in `fetchm2_all_assemblies.csv`.
- `--fetchm2_keep_assembly_duplicates` to pass FetchM2's `--keep-assembly-duplicates` option when users intentionally want paired GCA/GCF rows in `fetchm2_clean.csv`.
- `--fetchm2_download_engine` with native `fetchm2 seq` as the default downloader and the PanResistome downloader retained as an explicit/fallback path.
- FetchM2 is now the default metadata engine, adding richer standardized host, source, environment, geography, disease/health, collection-year, and metadata-audit fields while retaining a reversible legacy FetchM mode.
- FetchM2 compatibility adapter that preserves `fetchm2_clean.csv/.tsv`, writes legacy-compatible `ncbi_clean.csv`, records `metadata_engine.txt`, and organizes FetchM2 metadata analysis/audit outputs under each sample directory.
- FetchM2-specific filters for standardized sample type, isolation source, environment medium, and collection-year ranges.
- PanR2 skip precheck so runs with no valid downstream ABRicate/PanR2 inputs end cleanly and write `panr_output/panr2_input_status.txt`.
- Optional QUAST assembly-structure QC module with PanR2-ready assembly QC export.
- Optional FastANI/skani pairwise ANI module with ANI matrix, closest-genome table, near-duplicate clusters, ANI outlier report, and PanR2-ready ANI summary.
- Optional Mash sketch/distance pre-screen module for fast large-dataset screening.
- Combined QC decision engine writing `qc/qc_master_report.csv`, pass/fail/warning sample lists, and `qc/excluded_for_panr2.csv`.
- Optional representative-genome selection using ANI duplicate clusters with `--representative_only`.
- `panr2_inputs/` handoff export directory with metadata, AMR tables, QC reports, ANI summaries, QUAST summaries, manifests, and the formal feature contract columns.
- Formal PanR2 input contract documentation under `docs/panr2_input_contract.md`.
- New Conda environments for ANI (`envs/ani.yaml`), QUAST (`envs/quast.yaml`), and Mash (`envs/mash.yaml`).
- Offline `test` profile using tiny local fixtures.
- Python sequence-stat fallback for fixture-based CI tests.
- CI checks for offline sequence QC and metadata enrichment outputs.

### Fixed
- Relative `--checkm2_db_dir` values are now resolved from the launch directory instead of the Nextflow task work directory.
- CheckM2 automatic database downloads now retry transient download failures before failing the run.
- GitHub Actions CI now installs pinned Nextflow 24.10.4 through `nf-core/setup-nextflow@v2` instead of a brittle raw installer command.

### Validated
- Fresh remote-style run on `test.tsv` with FetchM2 0.1.7, native sequence download, CheckM2 automatic database download, QUAST, FastANI, Mash, QC filtering, comprehensive PanR2 mode, PanR2 handoff export, and GTDB-Tk disabled completed successfully.
- The validated run downloaded 10/10 selected assemblies, downloaded and verified the CheckM2 database under `<outdir>/databases/checkm2`, produced CheckM2 quality reports for 10 assemblies, and completed all 17 Nextflow processes.
- Comprehensive PanR2 mode generated NCBI AMR, VFDB, PlasmidFinder, MobileElementFinder, IntegronFinder, MLST, cross-database, temporal, QC, citation, software-version, and dashboard outputs.
- GitHub Actions CI passed for PanR2 `v0.1.3` and PanResistome after the hardened Nextflow setup update.

## 0.2.0 - 2026-05-01

### Added
- Sequence QC after FetchM using `seqkit stats`.
- CheckM2 completeness/contamination QC with low-memory mode enabled by default.
- Optional GTDB-Tk taxonomy consistency QC, disabled by default.
- QC-enriched `metadata_output/ncbi_enriched.csv`.
- Pass-only metadata and FASTA outputs for downstream filtering with `--qc_filter true`.
- Conda environment version reports under `pipeline_versions/`.
- Portable CheckM2 1.1.0 CPU environment using Python 3.12 and TensorFlow CPU 2.17.
- `scripts/bootstrap.sh` preflight helper for new users.
- `test_small.tsv` for lightweight full downstream validation.

### Changed
- Default metadata wording, environment version capture, README, and PanR2 handoff documentation now refer to FetchM2 as the preferred metadata source.
- ABRicate empty-output handling now writes parseable header-only summary/result files instead of malformed empty files.
- `--local_samples` now ignores support directories such as `pipeline_versions`, preventing local validation runs from treating version folders as sample inputs.
- Empty sequence-QC header writing is now safe for column names containing percent signs.
- `envs/fetchm.yaml` now installs PanR2 from the GitHub source so PanResistome uses the current PanR2 reporting layer during integrated runs.
- Default thread count is now 8.
- CheckM2 process resources default to 8 CPUs and 8 GB RAM.
- FetchM process now writes to an internal work output folder so any `--outdir` value works.

### Validated
- FetchM2 0.1.5 smoke validation on `test.tsv` with native `fetchm2 seq`; `fetchm2_clean.csv` contained 40 representative GCF rows, `fetchm2_all_assemblies.csv` contained 80 full rows, and 2/2 selected assemblies downloaded successfully.
- FetchM2 0.1.5 broader validation on `test.tsv` with native sequence download, CheckM2, QUAST, FastANI, Mash, QC filtering, ABRicate, PanR2, and PanR2 handoff export enabled; all 17 Nextflow processes completed successfully, 16/16 selected assemblies downloaded, `qc_master_status` reported 16 PASS and 24 FAIL rows, and `panr2_inputs/manifest/software_versions.csv` captured `fetchm2==0.1.5`.
- Earlier FetchM2 adapter validation on `test.tsv` with CheckM2, QUAST, FastANI, Mash, QC filtering, ABRicate, PanR2, and PanR2 handoff export enabled completed successfully before the FetchM2 0.1.5 representative-row update.
- `test.tsv` through FetchM, sequence QC, and CheckM2 with GTDB-Tk disabled.
- `test_small.tsv` through FetchM, sequence QC, CheckM2, ABRicate, PanR2, and result collection with `--qc_filter true`.
