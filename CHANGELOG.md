# Changelog

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
- Default thread count is now 8.
- CheckM2 process resources default to 8 CPUs and 8 GB RAM.
- FetchM process now writes to an internal work output folder so any `--outdir` value works.

### Validated
- `test.tsv` through FetchM, sequence QC, and CheckM2 with GTDB-Tk disabled.
- `test_small.tsv` through FetchM, sequence QC, CheckM2, ABRicate, PanR2, and result collection with `--qc_filter true`.
