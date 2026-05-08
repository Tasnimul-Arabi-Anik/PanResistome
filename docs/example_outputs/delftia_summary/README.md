# Delftia tsuruhatensis Example Output Summary

This directory documents a lightweight example of the outputs produced by a comprehensive PanResistome + PanR2 validation run on 45 current `Delftia tsuruhatensis` assemblies from NCBI.

The full run output is intentionally not committed because it contains downloaded assemblies, CheckM2 intermediate files, large figures, and machine-specific work products. Instead, this summary records the validated output structure and key counts so users know what a successful run should produce.

## Validation Input

Input file:

```text
validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv
```

Input records:

```text
45 NCBI Assembly rows
```

## Validated Command

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
  --threads 4 \
  --fetchm2_download_workers 2
```

In the local validation, CheckM2 used an already-downloaded database path. On a fresh run, PanResistome can also download/build supported databases where the upstream tool allows it. GTDB-Tk remains off by default because it requires a large external database.

## Main HTML Output

The main combined dashboard is generated at:

```text
<outdir>/Delftia_tsuruhatensis/report/index.html
```

The dashboard links to:

- QC reports
- metadata completeness and bias summaries
- AMR analysis
- VFDB virulence analysis
- PlasmidFinder analysis
- IntegronFinder analysis
- MLST analysis
- cross-database association outputs
- temporal trend outputs
- citations
- software versions

## Output Structure

A successful comprehensive run produces folders similar to:

```text
Delftia_tsuruhatensis/
├── amrfinderplus/
├── ani/
├── checkm2/
├── cross_database/
├── integronfinder/
├── mash/
├── metadata_analysis/
├── metadata_output/
├── mlst/
├── ncbi/
├── panr2_inputs/
├── plasmidfinder/
├── qc/
├── quast/
├── report/
├── sequence/
├── sequence_filtered/
├── sequence_qc/
├── temporal/
├── tool_results/
└── vfdb/
```

## QC Summary

Validated QC status for the 45-genome run:

| QC field | PASS |
| --- | ---: |
| sequence QC | 45 |
| CheckM2 QC | 45 |
| combined QC | 45 |
| ANI species consistency | 45 |
| QC master status | 45 |

The main QC table is:

```text
qc/qc_master_report.csv
```

## FetchM2 Metadata Summary

FetchM2 enriched the metadata with standardized fields used by downstream PanR2 analysis:

| Metadata field | Non-missing |
| --- | ---: |
| Country | 38/45 |
| Collection_Year | 42/45 |
| Host_SD | 11/45 |
| Sample_Type_SD | 25/45 |
| Environment_Medium_SD | 10/45 |
| Isolation_Source_SD | 2/45 |

Metadata reports are generated under:

```text
metadata_analysis/
qc/metadata_completeness_report.csv
qc/metadata_bias_warning.txt
```

## PanR2 Feature Tables

The standardized PanR2 handoff is generated under:

```text
panr2_inputs/
```

Validated feature contract summary:

```text
feature_files_checked=6
feature_rows=291
databases_seen=amr,amrfinderplus,integronfinder,mlst,plasmidfinder,vfdb
metadata_accessions=45
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Validated feature tables:

| Table | Rows |
| --- | ---: |
| `panr2_inputs/features/all_features.tsv` | 291 |
| `panr2_inputs/features/amr.features.tsv` | 19 |
| `panr2_inputs/features/amrfinderplus.features.tsv` | 24 |
| `panr2_inputs/features/vfdb.features.tsv` | 40 |
| `panr2_inputs/features/plasmidfinder.features.tsv` | 2 |
| `panr2_inputs/features/integronfinder.features.tsv` | 116 |
| `panr2_inputs/features/mlst.features.tsv` | 90 |

## Important Analysis Outputs

Cross-database outputs:

```text
cross_database/analysis/cross_database_feature_matrix.csv
cross_database/analysis/cross_database_top_associations.csv
cross_database/analysis/amr_integron_associations.csv
cross_database/analysis/amr_plasmid_associations.csv
cross_database/analysis/amr_virulence_associations.csv
cross_database/figures/html_files/cross_database_feature_network.html
cross_database/figures/html_files/global_feature_association_heatmap.html
cross_database/figures/html_files/integrated_feature_presence_heatmap.html
```

Temporal outputs:

```text
temporal/analysis/temporal_feature_trends.csv
temporal/analysis/temporal_burden_trends.csv
temporal/figures/html_files/temporal_top_feature_trends.html
```

Report and reproducibility outputs:

```text
report/index.html
report/citations.md
report/citations.bib
report/software_versions.csv
panr2_inputs/manifest/schema_validation_summary.txt
panr2_inputs/manifest/feature_completeness_audit.tsv
panr2_inputs/manifest/module_status_summary.tsv
```

## Notes

AMRFinderPlus database auto-download was validated and all 45 AMRFinderPlus sample calls passed. MobileElementFinder remains opt-in because upstream parser behavior can vary. ISfinder is supported through an authorized local FASTA supplied by the user, not automatic redistribution.

