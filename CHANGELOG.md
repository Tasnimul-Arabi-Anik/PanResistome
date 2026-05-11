# Changelog

## Unreleased

### Added
- v0.4.0 large-dataset report safeguards through `--large_dataset`, `--report_mode`, feature/metadata report caps, top-feature-per-database summaries, and `panr2_inputs/manifest/report_controls.tsv`.
- A `large` profile that enables compact report mode, report-facing feature caps, and large-dataset safeguards while preserving complete TSV exports.
- Draft HPC and container execution notes under `docs/hpc.md` and `docs/containers.md`.
- Large-mode validation evidence for the existing 100-record `Klebsiella pneumoniae` input under `validation/klebsiella_pneumoniae_100/LARGE_MODE_VALIDATION_RESULTS.md`.
- Parallel AMRFinderPlus sample execution through `--amrfinderplus_jobs` and `--amrfinderplus_threads_per_sample`, with per-sample status preserved in `amrfinderplus_sample_status.tsv`.
- Runtime/resource summary tables generated from the Nextflow trace as `pipeline_runtime_summary.tsv` and `pipeline_runtime_tasks.tsv`.
- A 300-record BioProject-diverse `Klebsiella pneumoniae` validation input and desktop-safe large-mode validation report under `validation/klebsiella_pneumoniae_300/`.
- Large-run ANI controls through `--ani_large_run_strategy auto|all|skip` and `--ani_max_all_vs_all_genomes`, with `ani/analysis/ani_run_status.tsv` documenting whether all-vs-all ANI ran or was skipped.
- An optional module validation matrix at `docs/optional_module_validation_matrix.md`, separating stable default modules, stable table-input paths, experimental runners, restricted database workflows, and remaining validation gaps.
- Fast optional-runner smoke validation evidence under `validation/optional_runner_smoke/`, using local fixtures and unrelated heavy/default modules disabled.
- First real biological optional-runner validation evidence under `validation/optional_runner_biological/`, using two complete Klebsiella genomes with Kleborate and MOB-suite enabled.

### Changed
- PanR2 handoff export now applies configured feature caps to presence/absence matrices and co-occurrence/proximity summaries, preserves complete proximity evidence as `feature_proximity_all.tsv`, and surfaces report-control settings in `panr2_inputs/report/report_controls.html`.
- AMRFinderPlus now runs independent assemblies concurrently by default using `--threads` as the sample-job count and one AMRFinderPlus thread per sample, avoiding the previous serial per-genome bottleneck.
- AMRFinderPlus per-sample execution now reuses existing non-empty TSV outputs, making interrupted large runs safer to resume without recomputing completed assemblies.
- AMRFinderPlus now prints periodic per-sample progress and records per-sample runtime in `amrfinderplus_sample_status.tsv`.
- In large-dataset mode, ANI `auto` strategy skips all-vs-all ANI above 200 genomes by default instead of silently launching a long FastANI/skani all-vs-all run; users can force full ANI with `--ani_large_run_strategy all`.
- PanR2 contract export now preserves non-ABRicate tool names for optional ABRicate-style table inputs such as ISfinder, MOB-suite, prophage/geNomad, DefenseFinder, and organism-specific typing tables.
- Optional runner CPU requests for MOB-suite, geNomad, and organism-specific typing now scale with `--threads`, with 1-CPU test-profile overrides for fast smoke validation.
- `--run_abricate false` can now bypass the legacy ABRicate/PanR branch when PanR2 comprehensive mode is disabled, allowing optional-runner/table-input smoke validation to export PanR2 inputs directly.
- PanR2 contract export now creates header-only feature tables with `WARNING_EMPTY` audit status for enabled optional modules that produced raw output but no biological feature rows.
- Kleborate, Kaptive, and ECTyper placeholder outputs now use valid tab-separated headers.
- Kleborate runner mode now uses the required `--preset kpsc` argument for Kleborate v3 and collects output tables from the Kleborate output directory.
- PanR2 contract export now converts real Kleborate output into standardized features for ST, virulence/resistance scores, K/O loci, siderophore markers, wzi, and AMR markers.
- MOB-suite runner mode now preserves per-sample stdout/stderr and writes `mobsuite/module_status.tsv`, making broken local installations visible in the PanR2 module status summary.
- MOB-suite environment creation now installs `mob-suite==3.1.9` with pip inside a Conda-managed environment to avoid Bioconda post-link `mob_init` fragility during environment creation.
- Added `--mobsuite_db` to pass a preinitialized MOB-suite database directory to `mob_recon --database_directory`.
- MOB-suite runner mode now uses a task-local writable `HOME` for ETE taxonomy cache creation and the database preflight distinguishes incomplete core MOB-suite databases from missing `taxa.sqlite`.
- MOB-suite runner mode now passes `--force` to `mob_recon`, and PanR2 contract export now converts real MOB-suite biomarker, plasmid-cluster, mobility, host-range, and MGE rows into standardized `mobsuite.features.tsv` rows.

### Tested
- Synthetic PanR2 contract tests now cover optional table-input exports for MobileElementFinder, ISfinder, MOB-suite, prophage/geNomad, DefenseFinder, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder.
- Optional-runner smoke validation completed on two local fixture assemblies with ISfinder-compatible BLAST, MOB-suite, geNomad/prophage, Kleborate, Kaptive, ECTyper, PanR2 export, and result collection enabled while CheckM2, GTDB-Tk, QUAST, ANI, Mash, AMRFinderPlus, PanR2 comprehensive mode, and legacy ABRicate/PanR were disabled.
- Kleborate biological validation completed on two complete `Klebsiella pneumoniae` assemblies, producing 25 standardized `kleborate.features.tsv` rows with zero unmatched, invalid, or duplicate feature rows.
- MOB-suite biological validation completed on the same two assemblies after preinitializing `taxa.sqlite`, processing 2/2 samples and producing 253 standardized `mobsuite.features.tsv` rows plus 25 Kleborate rows with zero unmatched, invalid, or duplicate feature rows.

### Validated
- The 100-record `Klebsiella pneumoniae` input completed with `-profile conda,mamba,desktop_parallel,large`, GTDB-Tk disabled, CheckM2 capped, AMRFinderPlus enabled, and parallel native feature runners. The run produced 12,838 complete standardized feature rows, 99 QC PASS calls, zero unmatched/invalid/duplicate feature rows, capped report-facing matrices/summaries, and the expected HTML report set.
- The 300-record `Klebsiella pneumoniae` large-mode validation completed with CheckM2, GTDB-Tk, ANI, and AMRFinderPlus disabled for desktop safety. It downloaded 299 genomes, produced 299 QC PASS calls, generated 36,638 standardized feature rows across AMR, VFDB, PlasmidFinder, IntegronFinder, and MLST, kept zero unmatched/invalid/duplicate feature rows, and generated compact large-mode reports.
- The same 300-record validation identified two scale-sensitive optional steps: FastANI all-vs-all and AMRFinderPlus nucleotide `tblastn`. Both should be treated as separate workstation/HPC or optimized-strategy targets for 300+ genome runs.

## 0.3.1 - 2026-05-09

### Added
- BioProject/study-bias reporting in the PanR2 handoff bundle, including `metadata_feature_analysis/bioproject_bias_report.tsv` and `report/bioproject_bias.html`.
- Confidence and warning labels for `metadata_feature_analysis/top_findings.tsv`, including support sample counts, BioProject dominance flags, and interpretation labels.
- Richer AMRFinderPlus-vs-ABRicate concordance summaries with normalized gene-symbol matching, class-level possible matches, and a dedicated `report/amrfinder_abricate_concordance.html` page.
- Explicit `evidence_level` and interpretation-warning columns in same-contig/proximity cross-database outputs.
- Tests covering BioProject dominance warnings, top-finding labels, AMRFinderPlus/ABRicate concordance, and proximity evidence columns.
- Klebsiella interpretation documentation covering BioProject bias, metadata usability, top-finding labels, concordance categories, and proximity evidence levels.

## 0.3.0 - 2026-05-09

### Added
- PanResistome-native PanR2 feature-runner stage controlled by `--panr2_native_feature_runners` for ABRicate, IntegronFinder, MLST, and opt-in MobileElementFinder execution before PanR2 reporting.
- Experimental native feature-runner backend selection with `--panr2_native_feature_runner_mode serial|parallel`; parallel mode runs ABRicate one database at a time with per-genome workers, then runs IntegronFinder/MLST per assembly concurrently within the PanResistome native-runner stage.
- Controlled 45-genome Delftia parallel-backend validation documented in `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md`; the run completed 9/9 processes with 16 threads, produced 201 standardized feature rows, preserved the same core feature calls as the serial native-runner validation, and kept zero unmatched, invalid, or duplicate feature rows.
- Separate `--checkm2_threads` parameter for CheckM2-only CPU control; it defaults to `min(--threads, 4)` so high native feature-runner parallelism does not automatically force high-memory CheckM2 DIAMOND runs.
- `scripts/run_panr2_native_features.py` to run standard PanR2-compatible annotation helpers under PanResistome ownership and write module status rows.
- Same-contig and coordinate-proximity cross-database outputs: `amr_mge_same_contig.tsv`, `amr_plasmid_same_contig.tsv`, `amr_integron_same_contig.tsv`, and `feature_proximity.tsv`.
- Metadata interpretation outputs for `feature_metadata_associations.tsv`, `database_burden_metadata_associations.tsv`, `category_burden_by_sample.tsv`, and `category_metadata_associations.tsv`.
- Lightweight PanR2 handoff HTML pages for top findings, metadata quality and bias, database burden by metadata, cross-database interpretation, and database setup/feature-contract status.
- `docs/roadmap_v0.3.0.md` and `validation/klebsiella_pneumoniae_100/README.md` to define the scale and interpretation validation target.
- Native feature-runner validation evidence under `validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md`.
- BioProject-diverse 100-record `Klebsiella pneumoniae` validation input under `validation/klebsiella_pneumoniae_100/`.
- `lowmem`, `desktop_parallel`, and `workstation` profiles for safer resource defaults, including independent CheckM2 thread caps.
- Native feature-runner merge audit at `panr2_inputs/manifest/native_runner_merge_audit.tsv`, recording expected vs observed raw table counts for ABRicate, IntegronFinder, MLST, and MobileElementFinder.
- v0.3.0 release-hardening docs: `docs/validation_matrix.md`, `docs/release_checklist_v0.3.0.md`, `docs/troubleshooting.md`, and `docs/example_klebsiella_interpretation.md`.

### Changed
- Comprehensive mode now passes precomputed ABRicate, IntegronFinder, and MLST directories to PanR2 by default instead of asking PanR2 to execute those standard runners internally.
- Native MLST runner status now counts only biological ST/allele features, not unsupported-organism placeholder rows emitted by the `mlst` command.
- CheckM2 now prints the effective CheckM2-only thread cap before execution and warns when the cap is high enough to risk desktop memory pressure.

### Fixed
- PanR2 feature-contract export now parses headerless native `mlst` output, suppresses placeholder calls such as `ST_-`, and writes a header-only `mlst.features.tsv` when MLST ran successfully but no valid ST/allele features were detected.
- `scripts/generate_ncbi_assembly_input.py` now accepts the documented `--limit`, `--out`, and `--prefer-refseq` options, writes organism-specific validation documentation instead of Delftia-specific boilerplate, and preserves an existing curated `README.md` by writing `INPUT_GENERATION.md`.
- `scripts/generate_ncbi_assembly_input.py` can now fetch a larger candidate pool and select validation inputs round-robin across BioProjects with `--diverse-bioproject`, reducing single-study dominance in second-organism validation inputs.

### Validated
- A 45-genome `Delftia tsuruhatensis` native feature-runner run completed all 20 Nextflow processes with GTDB-Tk disabled, CheckM2/QUAST/ANI/Mash enabled, AMRFinderPlus enabled, `--panr2_native_feature_runners true`, and 4 threads.
- The native feature-runner validation produced PanResistome-owned ABRicate, IntegronFinder, and MLST raw outputs before PanR2 reporting, then generated clean PanR2 handoff outputs with 201 standardized feature rows and zero unmatched, invalid, or duplicate feature rows.
- The same validation generated the v0.3.0 metadata interpretation, cross-database same-contig/proximity, feature matrix, handoff HTML, and top-level PanR2 report outputs.
- A 100-record `Klebsiella pneumoniae` parallel validation completed with GTDB-Tk disabled, CheckM2 capped at 2 threads, QUAST/ANI/Mash enabled, AMRFinderPlus enabled, and the parallel native feature-runner backend. The run downloaded 99 genomes, produced 99 QC PASS calls, generated 12,838 standardized feature rows across AMR, AMRFinderPlus, VFDB, PlasmidFinder, IntegronFinder, and MLST, and kept zero unmatched, invalid, or duplicate feature rows.

## 0.2.2 - 2026-05-08

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
- `scripts/summarize_validation_run.py` for writing compact CSV/Markdown validation summaries from generated PanR2 manifests.
- README module stability table clarifying stable, active-development, and experimental modules for public users.
- `docs/remote_user_validation.md` documenting the fresh-clone standard comprehensive validation path and release-passing criteria.
- `docs/release_reliability_checklist.md` defining the 20 release gates for a reliable fresh-user comprehensive workflow.
- Fresh-clone validation evidence under `validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md`.
- `pytest.ini` limiting test discovery to repository tests so local runs do not collect Nextflow work-directory Conda package tests.

### Changed
- Optional downstream annotation modules now run after the combined QC decision step, so `--qc_filter true` can use `sequence_filtered/` before AMRFinderPlus, MOB-suite, geNomad, and organism-specific typing.
- AMRFinderPlus database setup now runs by default before AMRFinderPlus execution, and the process records per-sample status instead of silently swallowing missing-database failures.
- CheckM2 QC now declares the same CPU count as `--threads` and allows longer wall time for full 45-genome laptop validation runs.
- MobileElementFinder is now opt-in through `--panr2_run_mobileelementfinder true`; real Delftia validation showed that the upstream `mefinder` BLAST JSON parser can fail on otherwise valid assemblies, so the public comprehensive default avoids a brittle hard failure.
- The FetchM2/QC Conda environment no longer installs PanR2 from GitHub; PanR2 is kept in the dedicated PanR2 comprehensive environment and the legacy `PANR` process now uses that environment.
- PanR2 contract export now suppresses placeholder MLST features such as `ST_-` and `-:ST-` so unknown sequence-type calls do not dominate top-finding summaries.

### Fixed
- Relative `--checkm2_db` values are now resolved from the launch directory before CheckM2 runs in a Nextflow work directory.
- Comprehensive PanR2 profiles no longer request the ABRicate `isfinder` database by default because standard ABRicate setup does not always provide it; users can still opt in with `--panr2_abricate_dbs ...isfinder` when that database is installed.
- `--qc_filter true` now fails at the combined QC decision step with a clear report when zero FASTA files remain for downstream analysis, instead of failing later inside PanR2.

### Validated
- A 45-genome `Delftia tsuruhatensis` run completed with FetchM2, sequence QC, CheckM2, QUAST, ANI, Mash, combined QC, AMRFinderPlus database auto-download, AMRFinderPlus, comprehensive PanR2, and PanR2 handoff export, with GTDB-Tk disabled and 4 threads.
- The validated PanR2 handoff contained AMR, AMRFinderPlus, VFDB, PlasmidFinder, IntegronFinder, and MLST feature tables in `panr2_inputs/features/all_features.tsv` with zero unmatched, invalid, or duplicate feature rows.
- A true fresh-clone validation at commit `bf241634d1e7bacf9539a5e891ab0329a863abb7` completed all 19 processes with no user-supplied CheckM2 database path, CheckM2 database auto-download, AMRFinderPlus database auto-update, ABRicate `ncbi/vfdb/plasmidfinder` setup verification, 45/45 QC PASS, 291 PanR2 feature rows, and a generated `report/index.html` dashboard.

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
