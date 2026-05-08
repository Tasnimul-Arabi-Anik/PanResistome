# Klebsiella pneumoniae 100-Genome Validation Plan

This validation is intended for PanResistome `v0.3.0` development. It is not run in CI because it downloads public assemblies and executes heavy bioinformatics tools.

## Generate Input

Use the NCBI Assembly E-utilities helper to create a FetchM2/PanResistome-compatible input table:

```bash
python scripts/generate_ncbi_assembly_input.py \
  --organism "Klebsiella pneumoniae" \
  --limit 100 \
  --prefer-refseq \
  --out validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv
```

Record the generation date, query, and row count before committing any validation input.

## Recommended Validation Command

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_100 \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --panr2_native_feature_runners true \
  --threads 4 \
  --fetchm2_download_workers 2
```

## Expected Checks

After completion, run:

```bash
scripts/summarize_validation_run.py \
  --run-dir validation_runs/klebsiella_100 \
  --out-dir validation_runs/klebsiella_100
```

Inspect:

```text
validation_runs/klebsiella_100/<organism>/report/index.html
validation_runs/klebsiella_100/<organism>/panr2_inputs/report/panr2_handoff_index.html
validation_runs/klebsiella_100/<organism>/panr2_inputs/features/all_features.tsv
validation_runs/klebsiella_100/<organism>/panr2_inputs/cross_database/feature_proximity.tsv
validation_runs/klebsiella_100/<organism>/panr2_inputs/metadata_feature_analysis/top_findings.tsv
validation_runs/klebsiella_100/validation_summary.md
```

## Why Klebsiella

`Klebsiella pneumoniae` is expected to stress more of the comparative genomics layer than the Delftia validation:

- richer AMR diversity,
- plasmid replicon diversity,
- virulence features,
- MLST relevance,
- optional Kleborate/Kaptive relevance,
- richer public host/source/country/year metadata.

## Success Criteria

- FetchM2 downloads assemblies successfully.
- CheckM2, QUAST, ANI/skani, and Mash complete with GTDB-Tk disabled.
- PanResistome-native feature runners produce ABRicate, IntegronFinder, and MLST directories.
- PanR2 report and PanR2 handoff report pages are generated.
- Feature contract validation has zero unmatched, invalid, or duplicate rows unless documented.
- Same-contig/proximity outputs are present, even if they contain zero rows for the selected dataset.

For v0.3.0 development, run this validation after the standard native runners are split into finer-grained parallel Nextflow processes. The 45-genome `Delftia tsuruhatensis` native-runner validation should be re-run first to confirm the parallel path preserves the same feature-contract outputs.
