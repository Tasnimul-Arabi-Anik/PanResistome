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

The native-runner Delftia validation showed that the current architecture is functionally valid but still serial inside `PANR2_FEATURE_RUNNERS`. ABRicate and IntegronFinder were the most important bottlenecks, so parallelizing those runners is the next engineering priority before larger organism validations.

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

Recommended order:

1. Re-run the 45-genome `Delftia tsuruhatensis` native-runner validation with `--panr2_native_feature_runner_mode parallel --threads 16`.
2. Compare serial and parallel outputs in `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md`.
3. Run the planned `Klebsiella pneumoniae` validation to test richer AMR, plasmid, VFDB, MLST, and metadata-feature associations.

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
