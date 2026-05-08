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
  --outdir validation_runs/delftia_fresh \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --threads 4 \
  --fetchm2_download_workers 2
```

The primary combined HTML output from a successful comprehensive run is:

```text
validation_runs/delftia_fresh/<organism>/report/index.html
```

Validated locally on 2026-05-08 with GTDB-Tk disabled, 4 threads, CheckM2 enabled with an existing local database path, AMRFinderPlus database auto-download enabled, and comprehensive PanR2 output generation. The run produced AMR, AMRFinderPlus, VFDB, PlasmidFinder, IntegronFinder, MLST, cross-database, temporal, QC, citation, software-version, and dashboard outputs.

Fresh-clone validation on 2026-05-08 also completed without a user-supplied CheckM2 database path. The pipeline automatically downloaded the CheckM2 database, updated AMRFinderPlus databases, verified ABRicate `ncbi/vfdb/plasmidfinder`, completed all 19 Nextflow processes, and generated a clean PanR2 handoff with 291 feature rows and zero unmatched, invalid, or duplicate feature rows. See `FRESH_CLONE_VALIDATION_RESULTS.md`.

Native feature-runner validation on 2026-05-08 completed all 20 Nextflow processes with `--panr2_native_feature_runners true`. PanResistome ran ABRicate, IntegronFinder, and MLST before PanR2 reporting, then exported 201 standardized feature rows with zero unmatched, invalid, or duplicate feature rows. MLST produced a valid header-only standardized table because this `Delftia` input had no supported PubMLST ST calls. See `NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md`.

The release-blocking fresh-clone validation path is documented in:

```text
docs/remote_user_validation.md
```
