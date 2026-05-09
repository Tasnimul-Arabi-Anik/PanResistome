# Klebsiella pneumoniae Validation Input

Generated on: 2026-05-09

NCBI E-utilities query: `"Klebsiella pneumoniae"[Organism]`
Assembly records written: `300`
Candidate records requested: `1200`
BioProject-diverse selection: `True`
RefSeq/GCF accessions sorted first: `true`

Files:
- `ncbi_dataset.tsv`: FetchM2/PanResistome-compatible Assembly TSV.
- `ncbi_assembly_esummary.json`: raw NCBI Assembly esummary JSON used to create the TSV.

Recommended desktop-safe large-mode validation command:

```bash
nextflow run main.nf -resume \
  --input validation/klebsiella_pneumoniae_300/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_300_large_checkm2_off_noamrfinder \
  -profile conda,mamba,desktop_parallel,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus false \
  --fetchm2_download_workers 2
```

This mode validates FetchM2 metadata/download, sequence QC, QUAST, Mash, native ABRicate/VFDB/PlasmidFinder, IntegronFinder, MLST, PanR2 large-mode report controls, and strict feature-contract export without the two longest desktop bottlenecks observed at 300 genomes: FastANI all-vs-all and AMRFinderPlus nucleotide `tblastn`.

For a full workstation/HPC benchmark, re-enable AMRFinderPlus and/or ANI intentionally:

```bash
nextflow run main.nf -resume \
  --input validation/klebsiella_pneumoniae_300/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_300_large_full_optional \
  -profile conda,mamba,workstation,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus true \
  --amrfinderplus_jobs 16 \
  --amrfinderplus_threads_per_sample 1 \
  --fetchm2_download_workers 2
```

Expect the AMRFinderPlus-enabled 300-genome run to take multiple hours on a desktop-class CPU.

The primary combined HTML output from a successful comprehensive run is:

```text
validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/<organism>/report/index.html
```

Validated result: `LARGE_MODE_CHECKM2_OFF_VALIDATION_RESULTS.md`.
