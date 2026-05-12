# Optional Module Validation Matrix

This matrix separates validated default behavior from optional and experimental modules. The rule for optional modules is:

```text
No module should be advertised as stable runner mode until it has:
raw output preservation,
module status or audit reporting,
PanR2-compatible feature export,
feature-completeness audit coverage,
and at least one documented validation route.
```

## Status Key

| Status | Meaning |
| --- | --- |
| Stable | Validated in real or release-gate runs and suitable for normal use. |
| Stable table input | Precomputed tables are accepted and exported into the PanR2 contract. Runner mode may still be experimental. |
| Experimental runner | Runner exists but dependency/database setup or upstream behavior is not yet validated enough for default use. |
| Restricted | Requires user-supplied authorized data or database paths; PanResistome must not download or redistribute the database. |
| Planned | Documented roadmap item, not a supported run path yet. |

## Validation Matrix

| Module | Runner status | Table-input status | Database required | Validated route | PanR2 feature export | Default? | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ABRicate NCBI AMR | Stable | Stable | Bundled/setup ABRicate DB | Delftia fresh-clone, native runner, Klebsiella 100/300 | `amr.features.tsv` | Comprehensive mode | Primary assembled-genome AMR screen. |
| ABRicate VFDB | Stable | Stable | Bundled/setup ABRicate DB | Delftia/Klebsiella validation | `vfdb.features.tsv` | Comprehensive mode | Main virulence feature screen. |
| ABRicate PlasmidFinder | Stable | Stable | Bundled/setup ABRicate DB | Delftia/Klebsiella validation | `plasmidfinder.features.tsv` | Comprehensive mode | Replicon screen, not plasmid reconstruction. |
| AMRFinderPlus | Stable optional | Stable | Auto-updated AMRFinderPlus DB | Delftia/Klebsiella validation; 300-genome run documents runtime limit | `amrfinderplus.features.tsv` | No | Strong AMR module, but nucleotide mode can be slow on 300+ genomes. |
| IntegronFinder | Stable in native runner | Stable | Tool environment | Delftia/Klebsiella native runner validation | `integronfinder.features.tsv` | Comprehensive mode | Runtime can grow with fragmented assemblies. |
| MLST | Stable in native runner | Stable | Tool schemes bundled by `mlst` | Delftia/Klebsiella native runner validation | `mlst.features.tsv` | Comprehensive mode | Unsupported organisms produce header-only/no-call feature tables. |
| MobileElementFinder | Experimental runner | Stable table input | CGE tool/database environment | Synthetic PanR2 contract test; broad optional-module table validation; runner kept opt-in | `mobileelementfinder.features.tsv` | No | Upstream parser failures were observed on valid assemblies. Failures are nonfatal by default through `--panr2_mobileelementfinder_allow_failure true`, with auditable header-only outputs. |
| ISfinder-compatible BLAST | Restricted runner | Stable table input | User-supplied authorized FASTA | Unit test for BLAST converter; optional-runner smoke with synthetic local FASTA | `isfinder.features.tsv` | No | PanResistome does not download or redistribute ISfinder. Smoke validates orchestration, not biological ISfinder calls. |
| MOB-suite | Biologically validated 100-genome runner | Stable table input | MOB-suite environment/database plus ETE `taxa.sqlite` | Synthetic PanR2 contract test; optional-runner smoke; Klebsiella 2-genome and 100-genome biological validation | `mobsuite.features.tsv` | No | The repo env now installs `mob-suite` through pip to avoid Bioconda post-link fragility. The 100-genome parallel validation processed 100/100 genomes, produced 16,332 PanR2 feature rows, and reduced the MOB-suite stage from ~1h25m to 55.63m. |
| geNomad/prophage | Experimental runner with auto-download helper | Stable table input | geNomad DB, downloaded when enabled and no DB path is supplied | Synthetic PanR2 contract test; broad optional-module table validation; optional-runner smoke without DB | `prophage.features.tsv` | No | PanResistome supports `--genomad_auto_download_db true` into `--genomad_db_dir`, but full biological runner validation still needs a local/fresh-downloaded DB. |
| DefenseFinder | Table-input recommended | Stable table input | DefenseFinder environment/database | Synthetic PanR2 contract test; broad optional-module table validation | `defensefinder.features.tsv` | No | Not part of default comprehensive mode until dependency/database setup is stable in a PanResistome-owned runner. |
| Kleborate | Biologically validated 100-genome runner | Stable table input | Kleborate environment | Synthetic PanR2 contract test; Klebsiella 2-genome and 100-genome biological validation | `kleborate.features.tsv` | No | Kleborate produced 1,158 PanR2 feature rows across the 100-genome Klebsiella validation with clean schema validation; per-genome parallel validation reduced the typing stage to 3.78m. |
| Kaptive | Experimental runner | Stable table input | User-supplied Kaptive DB | Synthetic PanR2 contract test; optional-runner smoke without DB | `kaptive.features.tsv` | No | Requires explicit `--kaptive_db` for runner mode. |
| ECTyper | Experimental runner | Stable table input | ECTyper environment/database | Synthetic PanR2 contract test; optional-runner smoke without installed tool | `ectyper.features.tsv` | No | Relevant mainly for E. coli validation. |
| SerotypeFinder | Planned runner | Stable table input | CGE database if runner is added later | Synthetic PanR2 contract test | `serotypefinder.features.tsv` | No | Table-input only for now. |
| SCCmecFinder | Planned runner | Stable table input | CGE database if runner is added later | Synthetic PanR2 contract test | `sccmecfinder.features.tsv` | No | Table-input only for now. |
| GTDB-Tk | Stable heavy optional | Metric/taxonomy output | User-supplied GTDB-Tk DB | Partial/local validation only | Metrics/handoff, not feature-like by default | No | Remains disabled by default because of database size. |
| QUAST | Stable optional | Metrics | None or optional reference | Delftia/Klebsiella validation | `assembly_metrics.tsv` | No | Assembly structure QC. |
| ANI/skani/FastANI | Stable optional with large-run guard | Metrics/features | Tool environment | Delftia/Klebsiella validation; 300-run guard | `ani.features.tsv`/metrics | No | All-vs-all ANI is skipped automatically in large mode above threshold unless forced. |
| Mash | Stable optional | Metrics | Tool environment | Delftia/Klebsiella validation | `mash_metrics.tsv` | No | Fast screening layer, not final taxonomy. |

## Optional Table-Analysis Validation

The repository also includes a reproducible optional-feature analysis validation under `validation/optional_feature_analysis/`, generated by:

```bash
python scripts/validate_optional_feature_analysis.py \
  --outdir validation/optional_feature_analysis \
  --force
```

This validation uses small local tables for AMR, VFDB, PlasmidFinder, MobileElementFinder, ISfinder-style BLAST, MOB-suite, prophage/geNomad, DefenseFinder, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder. It confirms that optional outputs can enter the same PanR2 analysis layer as AMR/VFDB:

```text
feature_files_checked=13
feature_rows=23
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

The run generated `all_features.tsv`, individual `*.features.tsv` files, feature matrices, co-occurrence/proximity outputs, metadata usability outputs, top-feature summaries, top findings, and HTML handoff pages. It validates the standardized table-to-analysis path, not every external runner/database installation.

## Current Validation Gaps

The next practical validation work should focus on targeted small real runs, not new modules:

1. **MobileElementFinder runner retry:** run on a small, controlled 5-10 genome subset and document parser behavior.
2. **Kaptive targeted Klebsiella subset:** run with a local Kaptive database path and verify capsule/O-locus feature export.
3. **geNomad runner smoke:** run with a local or freshly downloaded geNomad database path; verify prophage feature export.
4. **ECTyper targeted E. coli subset:** table-input first, runner mode later.

Large organism validations should use these modules only after small targeted smoke tests pass.

The optional feature-analysis result is documented at `validation/optional_feature_analysis/VALIDATION_RESULTS.md`.
The broad optional-module status for geNomad/prophage, DefenseFinder, and MobileElementFinder is documented at `validation/optional_runner_biological/BROAD_OPTIONAL_MODULES_STATUS.md`.
The current optional-runner smoke result is documented at `validation/optional_runner_smoke/OPTIONAL_RUNNER_SMOKE_RESULTS.md`.
The first real biological optional-runner result is documented at `validation/optional_runner_biological/KLEBSIELLA_2_OPTIONAL_RUNNER_RESULTS.md`.
The 100-genome MOB-suite/Kleborate biological validation and parallel-runner revalidation are documented at `validation/optional_runner_biological/KLEBSIELLA_100_MOBSUITE_KLEBORATE_RESULTS.md`.
