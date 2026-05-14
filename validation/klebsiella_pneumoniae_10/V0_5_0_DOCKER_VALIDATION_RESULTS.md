# PanResistome v0.5.0 Compact Docker/GHCR Biological Validation

Date: 2026-05-14

Commit tested: `0f08ec1`

Input:

- Organism: `Klebsiella pneumoniae`
- Input records: 10
- Input file: `validation/klebsiella_pneumoniae_10/ncbi_dataset.tsv`
- Container image: `ghcr.io/tasnimul-arabi-anik/panresistome:experimental`
- Container engine: Docker through the Nextflow `docker` profile

Command:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_10/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_10_v050_docker \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus true \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 8 \
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_v050_kleb10_work
```

## Result

Status: PASS

Nextflow summary:

- Processes succeeded: 14/14
- Duration: 18m40s
- CPU hours: 2.7
- Runtime summary: `validation_runs/klebsiella_10_v050_docker/pipeline_runtime_summary.tsv`

Processes completed:

- `FETCHM_ENV_VERSIONS`
- `ABRICATE_ENV_VERSIONS`
- `PANR2_COMPREHENSIVE_ENV_VERSIONS`
- `AMRFINDERPLUS_ENV_VERSIONS`
- `MASH_ENV_VERSIONS`
- `FETCHM`
- `SEQUENCE_QC`
- `MASH_PRESCREEN`
- `COMBINED_QC`
- `AMRFINDERPLUS_ANALYSIS`
- `PANR2_FEATURE_RUNNERS`
- `PANR2_COMPREHENSIVE`
- `EXPORT_PANR2_INPUTS`
- `COLLECT_RESULTS`

## Feature Contract

PanR2 schema validation:

```text
feature_files_checked=6
feature_rows=1348
databases_seen=amr,amrfinderplus,mlst,plasmidfinder,vfdb
samples_seen=20
metadata_accessions=10
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature rows by table:

| Feature table | Rows excluding header |
| --- | ---: |
| `amr.features.tsv` | 132 |
| `amrfinderplus.features.tsv` | 144 |
| `integronfinder.features.tsv` | 0 |
| `mlst.features.tsv` | 170 |
| `plasmidfinder.features.tsv` | 42 |
| `vfdb.features.tsv` | 860 |
| `all_features.tsv` | 1348 |

The `integronfinder.features.tsv` file was header-only in this compact subset,
while the native runner itself completed successfully and wrote audited raw
outputs.

## Database and Runner Status

Database setup status highlights:

- FetchM2 metadata: PASS
- Sequence FASTA inputs: PASS, 10 FASTA files
- Mash: PASS
- PanR2: PASS
- ABRicate setup report: PASS
- ABRicate databases `ncbi`, `vfdb`, `plasmidfinder`: PASS
- IntegronFinder: PASS
- MLST: PASS
- AMRFinderPlus database and runner: PASS for 10 samples
- CheckM2, QUAST, ANI, GTDB-Tk: SKIPPED by validation design
- geNomad, MOB-suite, DefenseFinder, MobileElementFinder, ISfinder, Kaptive: SKIPPED/not requested

Native runner merge audit:

| Module | Mode | Expected raw tables | Observed raw tables | Samples processed | Status |
| --- | --- | ---: | ---: | ---: | --- |
| ABRicate | parallel | 30 | 30 | 10 | PASS |
| IntegronFinder | parallel | 10 | 10 | 10 | PASS |
| MLST | parallel | 1 | 1 | 10 | PASS |
| MobileElementFinder | parallel | 0 | 0 | 0 | SKIPPED |

## v0.5.0 Outputs Verified

The run generated the new v0.5.0 interpretation outputs:

- `panr2_inputs/manifest/reproducibility_manifest.json`
- `panr2_inputs/manifest/feature_contract.json`
- `panr2_inputs/metadata_feature_analysis/lineage_summary.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_feature_burden.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_metadata_confounding.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_adjusted_warnings.tsv`
- `panr2_inputs/metadata_feature_analysis/statistical_summary.tsv`
- `panr2_inputs/diversity/feature_richness_by_sample.tsv`
- `panr2_inputs/diversity/database_diversity_by_sample.tsv`
- `panr2_inputs/diversity/jaccard_distance_matrix.tsv`
- `panr2_inputs/diversity/core_accessory_rare_features.tsv`
- `panr2_inputs/diversity/pan_feature_accumulation.tsv`
- `panr2_inputs/report/lineage_context.html`
- `panr2_inputs/report/diversity_summary.html`
- `panr2_inputs/report/statistical_summary.html`

Output counts:

- `top_findings.tsv`: 50 findings plus header
- `lineage_summary.tsv`: 10 samples plus header
- `lineage_adjusted_warnings.tsv`: 50 findings plus header
- `core_accessory_rare_features.tsv`: 274 features plus header
- `pan_feature_accumulation.tsv`: 10 samples plus header
- Report HTML files in `panr2_inputs/report/`: 14

Example lineage context:

```text
assembly_accession  mlst_ST   ani_cluster  bioproject   lineage_data_status
GCF_041200225.2     ST_23                  PRJNA1143699 available
GCF_041200245.2     ST_23                  PRJNA224116  available
GCF_048279315.2     ST_14464               PRJNA1204189 available
```

The compact subset had MLST lineage context but no ANI clusters because ANI was
disabled by design.

## Large-Mode Report Controls

`report_controls.tsv` recorded:

```text
large_dataset=true
report_mode=compact
max_features_heatmap=150
max_features_network=150
max_metadata_columns=20
top_n_features_per_database=50
skip_heavy_interactive_plots=true
samples=10
feature_rows=1348
unique_features=274
databases=amr,amrfinderplus,mlst,plasmidfinder,vfdb
```

Complete feature tables remained available under
`panr2_inputs/features/`; large-mode caps affected report-facing summaries.

## Bottleneck Observed

AMRFinderPlus was the dominant runtime cost:

```text
AMRFINDERPLUS_ANALYSIS: 13.93m real time
PANR2_FEATURE_RUNNERS: 2.43m real time
PANR2_COMPREHENSIVE:   56.40s real time
EXPORT_PANR2_INPUTS:    3.00s real time
```

The AMRFinderPlus process actively ran `tblastn` jobs and was not stalled. For
larger validations, AMRFinderPlus should remain an intentionally enabled stage
with explicit resource planning.

## Interpretation

This compact Docker/GHCR biological validation confirms that the v0.5.0
reproducibility/contract manifests and the lineage, diversity, and statistical
summary outputs work on real Klebsiella feature-rich data while preserving clean
PanR2 feature-contract validation.

This is not a replacement for the 100-record and 300-record validations. It is a
bounded release-gate style validation for the new v0.5.0 interpretation layer.
