# Klebsiella pneumoniae Validation Input

Generated on: 2026-05-08

NCBI E-utilities query: `"Klebsiella pneumoniae"[Organism]`
Assembly records written: `100`
Candidate records requested: `500`
BioProject-diverse selection: `True`
RefSeq/GCF accessions sorted first: `true`

Files:
- `ncbi_dataset.tsv`: FetchM2/PanResistome-compatible Assembly TSV.
- `ncbi_assembly_esummary.json`: raw NCBI Assembly esummary JSON used to create the TSV.

Recommended validation command:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_pneumoniae_current \
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

The primary combined HTML output from a successful comprehensive run is:

```text
validation_runs/klebsiella_pneumoniae_current/<organism>/report/index.html
```
