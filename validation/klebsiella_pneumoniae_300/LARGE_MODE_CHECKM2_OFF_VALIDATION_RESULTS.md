# Klebsiella pneumoniae 300-Genome Large-Mode Validation

Date: 2026-05-10

Repository commit at run start: `e704791`

This validation tested the v0.4.0 large-dataset report path on a 300-record, BioProject-diverse `Klebsiella pneumoniae` input. The goal was to validate large report safeguards and native feature-runner scale on a desktop-safe profile, not to run every heavy optional module.

## Input

- Input file: `validation/klebsiella_pneumoniae_300/ncbi_dataset.tsv`
- Input records: 300
- Candidate Assembly records queried: 1200
- BioProject-diverse selection: enabled
- Preferred RefSeq/GCF accessions: enabled

## Command

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

## Result

Status: PASS

- Duration: 43m 12s
- CPU hours: 18.2
- Nextflow processes succeeded: 8
- Cached processes: 6
- Metadata rows after download: 299
- Downloaded FASTAs: 299
- Failed accession: `GCF_055382775.1`
- QC master status: 299 PASS
- CheckM2: skipped for this desktop-safe large-mode run
- GTDB-Tk: skipped
- ANI/FastANI: skipped after initial bottleneck observation
- AMRFinderPlus: skipped after initial bottleneck observation
- Dashboard/report: generated

## PanR2 Feature Contract

Schema validation passed.

- Feature rows: 36,638
- Databases seen: `amr`, `vfdb`, `plasmidfinder`, `integronfinder`, `mlst`
- Unmatched feature rows: 0
- Invalid feature rows: 0
- Duplicate feature rows: 0

Feature rows by table:

| Table | Feature rows |
| --- | ---: |
| `amr.features.tsv` | 3,487 |
| `vfdb.features.tsv` | 25,726 |
| `plasmidfinder.features.tsv` | 1,310 |
| `integronfinder.features.tsv` | 1,128 |
| `mlst.features.tsv` | 4,987 |
| `all_features.tsv` | 36,638 |

## Native Runner Audit

`panr2_inputs/manifest/native_runner_merge_audit.tsv` passed.

| Module | Expected raw tables | Observed raw tables | Samples processed | Feature rows | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| ABRicate | 897 | 897 | 299 | 30,523 | PASS |
| IntegronFinder | 299 | 299 | 299 | 1,128 | PASS |
| MLST | 1 | 1 | 299 | 2,346 | PASS |
| MobileElementFinder | 0 | 0 | 0 | 0 | SKIPPED |

## Large-Mode Report Controls

`panr2_inputs/manifest/report_controls.tsv` recorded:

| Setting | Value |
| --- | --- |
| `large_dataset` | `true` |
| `report_mode` | `compact` |
| `max_features_heatmap` | `150` |
| `max_features_network` | `150` |
| `max_metadata_columns` | `20` |
| `top_n_features_per_database` | `50` |
| `skip_heavy_interactive_plots` | `true` |
| `feature_rows` | `36638` |
| `unique_features` | `758` |

Complete feature tables were preserved under `panr2_inputs/features/`. Report-facing matrices and proximity/co-occurrence summaries were capped for readability. Complete proximity evidence was preserved as `panr2_inputs/cross_database/feature_proximity_all.tsv`.

## Runtime Notes

The completed no-AMRFinderPlus large-mode run showed these main process timings:

| Process | Runtime | Peak RSS |
| --- | ---: | ---: |
| `PANR2_FEATURE_RUNNERS` | 28.68m | 4.5 GiB |
| `PANR2_COMPREHENSIVE` | 8.28m | 7.1 GiB |
| `QUAST_QC` | 3.10m | 1.3 GiB |
| `MASH_PRESCREEN` | 1.93m | 0.04 GiB |
| `EXPORT_PANR2_INPUTS` | 32s | 0.44 GiB |

FetchM2 was cached in the final run; the cached trace recorded the original FetchM process at 22.28m.

## Bottleneck Findings

Two heavier options were intentionally excluded from the completed validation after direct observation:

1. FastANI all-vs-all became the first major bottleneck at 300 genomes. It was still processing query genomes one by one after several minutes and was stopped before completion. For 300+ genome validations, `--run_ani false` is the practical desktop default until ANI is chunked, prescreened, switched to a faster strategy, or run on a workstation/HPC node.

2. AMRFinderPlus nucleotide mode is CPU-bound through `tblastn`. In a preliminary 300-genome attempt with `--amrfinderplus_jobs 8 --amrfinderplus_threads_per_sample 1`, completed samples took roughly 5-8 minutes each. RAM use was safe, but wall time was multi-hour. This should be run as a separate workstation/HPC or overnight benchmark for 300+ genomes.

## Key Outputs

- Main report: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/Klebsiella_pneumoniae/report/index.html`
- Handoff report: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/Klebsiella_pneumoniae/panr2_inputs/report/panr2_handoff_index.html`
- Feature contract: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/Klebsiella_pneumoniae/panr2_inputs/features/all_features.tsv`
- Report controls: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/Klebsiella_pneumoniae/panr2_inputs/manifest/report_controls.tsv`
- Native runner audit: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/Klebsiella_pneumoniae/panr2_inputs/manifest/native_runner_merge_audit.tsv`
- Runtime summary: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/pipeline_runtime_summary.tsv`
- Compact validation summary: `validation_runs/klebsiella_300_large_checkm2_off_noamrfinder/validation_summary.md`

## Interpretation

This run validates the large-mode PanR2 report safeguards and native PanResistome feature-runner path on a 300-record dataset with CheckM2 disabled. It does not replace the 100-record comprehensive validation with CheckM2 and AMRFinderPlus enabled. Instead, it defines a practical desktop-scale mode and identifies the next scalability targets: ANI strategy, AMRFinderPlus large-run handling, and further splitting/chunking of native feature runners.
