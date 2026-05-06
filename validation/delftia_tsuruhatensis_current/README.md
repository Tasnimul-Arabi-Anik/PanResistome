# Delftia tsuruhatensis Validation Input

Generated on: 2026-05-07

NCBI E-utilities query: `"Delftia tsuruhatensis"[Organism]`
Assembly records written: `45`

Files:
- `ncbi_dataset.tsv`: FetchM2/PanResistome-compatible Assembly TSV.
- `ncbi_assembly_esummary.json`: raw NCBI Assembly esummary JSON used to create the TSV.

Recommended validation command:

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_current \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --threads 8
```

The primary combined HTML output from a successful comprehensive run is:

```text
validation_runs/delftia_current/<organism>/report/index.html
```
