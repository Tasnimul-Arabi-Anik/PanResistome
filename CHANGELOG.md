# Changelog

## Unreleased

### Added
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
- `--local_samples` now ignores support directories such as `pipeline_versions`, preventing local validation runs from treating version folders as sample inputs.
- Empty sequence-QC header writing is now safe for column names containing percent signs.
- `envs/fetchm.yaml` now installs PanR2 from the GitHub source so PanResistome uses the current PanR2 reporting layer during integrated runs.
- Default thread count is now 8.
- CheckM2 process resources default to 8 CPUs and 8 GB RAM.
- FetchM process now writes to an internal work output folder so any `--outdir` value works.

### Validated
- `test.tsv` through FetchM, sequence QC, and CheckM2 with GTDB-Tk disabled.
- `test_small.tsv` through FetchM, sequence QC, CheckM2, ABRicate, PanR2, and result collection with `--qc_filter true`.
