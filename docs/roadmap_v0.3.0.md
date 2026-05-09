# PanResistome v0.3.0 Roadmap

Theme: scalable metadata-aware comparative genomics.

PanResistome `v0.2.2` established that a fresh user can run the standard comprehensive workflow end to end. The `v0.3.0` target is to make the workflow more scalable, easier to interpret, and closer to manuscript-ready comparative genomics.

## 1. Architecture

- Keep PanResistome as the heavy execution layer.
- Keep PanR2 as the standardized feature-table analysis and reporting layer.
- Prefer `panr2_inputs/features/*.features.tsv` and `panr2_inputs/features/all_features.tsv` as the strict downstream contract.
- Preserve raw tool outputs for auditability.
- Keep fragile or restricted modules opt-in.

Implemented groundwork after `v0.2.2`:

- `--panr2_native_feature_runners true` runs ABRicate, IntegronFinder, MLST, and optional MobileElementFinder under PanResistome before PanR2.
- `--panr2_native_feature_runner_mode serial|parallel` preserves the validated serial path while exposing an experimental parallel backend that runs each ABRicate database with per-genome workers, then runs per-assembly IntegronFinder/MLST execution.
- PanR2 receives precomputed `--abricate-dir`, `--vfdb-dir`, `--plasmidfinder-dir`, `--integronfinder-dir`, and `--mlst-dir` inputs rather than being asked to run those tools internally.
- A 45-genome `Delftia tsuruhatensis` validation passed with native ABRicate, IntegronFinder, and MLST runners enabled before PanR2 reporting. Results are documented in `validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md`.
- The parallel native-runner backend passed the 45-genome Delftia comparison and the 100-record Klebsiella validation. See `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md` and `validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md`.

## 2. Scalability

Current target:

- Split large per-assembly annotation work into finer-grained Nextflow channels.
- Keep the current sample-directory process as a stable fallback.
- Add runtime summaries that identify the slowest modules and samples.
- Add large-dataset plotting safeguards.

Near-term candidates:

- Per-database and eventually per-assembly ABRicate channels.
- Per-assembly IntegronFinder channel.
- Per-assembly MLST channel.
- Optional representative-only downstream analysis after ANI duplicate clustering.

The native-runner Delftia validation showed that the architecture is functionally valid, and the parallel backend now provides a validated speed path for ABRicate, IntegronFinder, and MLST while preserving the serial fallback.

Implemented resource-safety work:

- `--checkm2_threads` caps CheckM2 independently from general `--threads`.
- `lowmem`, `desktop_parallel`, and `workstation` profiles provide safer defaults for common local machines.
- `PANR2_FEATURE_RUNNERS` writes `native_runner_merge_audit.tsv` so expected and observed raw table counts can be audited after serial or parallel runs.

## 3. Metadata Interpretation

PanResistome/PanR2 should answer:

- Which metadata groups carry higher AMR, VFDB, plasmid, MGE, integron, MLST, or AMRFinderPlus burden?
- Which features are enriched in a country, host, sample type, isolation source, environment, year, or BioProject?
- Which findings are likely too sparse or biased to trust?

Implemented groundwork:

- `feature_metadata_associations.tsv`
- `database_burden_metadata_associations.tsv`
- `category_burden_by_sample.tsv`
- `category_metadata_associations.tsv`
- `top_findings.html`
- `metadata_quality_and_bias.html`
- `database_burden_by_metadata.html`

## 4. Cross-Database Biology

Co-occurrence must be interpreted in tiers:

- Level 1: same genome/sample.
- Level 2: same contig.
- Level 3: same contig within 10 kb.
- Level 4: overlapping or adjacent coordinates.

Implemented groundwork:

- `amr_mge_same_contig.tsv`
- `amr_plasmid_same_contig.tsv`
- `amr_integron_same_contig.tsv`
- `feature_proximity.tsv`
- `cross_database_interpretation.html`

These outputs do not prove transfer, expression, phenotype, or plasmid localization. They provide stronger context than sample-level co-occurrence and should be interpreted with assembly fragmentation and annotation limits in mind.

## 5. Second-Organism Validation

The next validation organism should be biologically richer than `Delftia tsuruhatensis`.

Preferred:

- `Klebsiella pneumoniae`

Alternatives:

- `Escherichia coli`
- `Salmonella enterica`
- `Acinetobacter baumannii`
- `Pseudomonas aeruginosa`
- `Staphylococcus aureus`

Validation target:

- Around 100 genomes.
- GTDB-Tk disabled initially.
- CheckM2 enabled.
- QUAST, ANI/skani, and Mash enabled.
- AMRFinderPlus enabled.
- Comprehensive PanR2 enabled.
- Optional Kleborate/Kaptive table-input or runner validation when databases are available.

Completed v0.3.0 validation:

- The 45-genome Delftia parallel comparison passed.
- The 100-record Klebsiella validation passed with 99 downloaded genomes, 99 QC PASS genomes, 12,838 feature rows, and zero unmatched, invalid, or duplicate PanR2 feature rows.

Next validation targets:

1. Repeat Klebsiella without `-resume` for final release evidence if resources allow.
2. Add a larger 300-500 genome stress test after v0.3.0.
3. Add container/HPC validation after deployment profiles are implemented.

## 6. Deployment

Conda/Mamba remains the validated public path. Future deployment work should add:

- Apptainer/Singularity profile.
- Docker profile.
- SLURM profile.
- Low-memory profile.
- Large-dataset profile.
- Containerized validation evidence before advertising those profiles as stable.

## 7. Release Rules

- Patch releases: bug fixes, documentation, parser fixes, and validation evidence.
- Minor releases: new validated workflow capabilities.
- Tags must remain immutable.
- `docs/release_reliability_checklist.md` must be updated before each release.
