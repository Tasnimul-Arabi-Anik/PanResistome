# Klebsiella 100-Genome MOB-suite/Kleborate Optional Runner Validation

Date: 2026-05-12

Purpose: validate whether real opt-in runner outputs can produce PanR2-compatible analysis at 100-genome scale, similar to the standard AMR/VFDB feature-analysis path.

Input: 100 `Klebsiella pneumoniae` assemblies staged from the existing 300-genome validation cache under:

```text
validation_runs/optional_real_100_input/Klebsiella_pneumoniae
```

The staged input is intentionally not committed because it contains copied genome FASTA files.

## Command

```bash
nextflow run main.nf \
  -resume \
  -profile conda,mamba,desktop_parallel,large \
  --local_samples validation_runs/optional_real_100_input \
  --outdir validation_runs/optional_real_100_klebsiella_mobsuite_kleborate \
  --sequence_qc_engine python \
  --qc_filter true \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive false \
  --run_abricate false \
  --export_panr2_inputs true \
  --run_mobsuite true \
  --mobsuite_db validation_runs/mobsuite_db_test \
  --run_kleborate true \
  --run_kaptive false \
  --run_ectyper false \
  --run_genomad false \
  --threads 8 \
  --capture_versions false
```

The run intentionally disabled default AMR/VFDB/PlasmidFinder, CheckM2, ANI, Mash, AMRFinderPlus, QUAST, and GTDB-Tk. This isolates optional runner behavior and PanR2 feature-contract export.

## Result

Status: PASS

Completed processes:

```text
SEQUENCE_QC
COMBINED_QC
MOBSUITE_ANALYSIS
ORGANISM_SPECIFIC_TYPING
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

Runtime:

```text
Duration: 1h 44m 49s
CPU hours: 13.7
```

The MOB-suite module status reported:

```text
samples_input=100
samples_processed=100
samples_failed=0
raw_tables_created=100
status=PASS
```

The MOB-suite database audit reported:

```text
core_status=PASS
taxa_status=PASS
status=PASS
message=core MOB-suite database files are present; taxa.sqlite is present
```

## PanR2 Feature Contract Quality

Schema validation summary:

```text
feature_files_checked=2
feature_rows=17490
databases_seen=kleborate,mobsuite
samples_seen=100
metadata_accessions=100
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature completeness audit:

| Database | Feature rows | Unique features | Samples with features | Status |
| --- | ---: | ---: | ---: | --- |
| `mobsuite` | 16,332 | 489 | 100 | PASS |
| `kleborate` | 1,158 | 263 | 96 | PASS |

The four genomes without Kleborate feature rows produced no biological Kleborate calls after parsing; this did not create invalid or unmatched rows.

## MOB-suite Feature Quality

MOB-suite produced standardized plasmid/reconstruction-context features across all 100 genomes.

Top MOB-suite feature categories:

| Category | Rows |
| --- | ---: |
| molecule_type | 4,597 |
| mate_pair_formation | 1,764 |
| plasmid_cluster | 1,486 |
| plasmid_cluster_secondary | 1,396 |
| replicon | 1,302 |
| is5 | 1,086 |
| is3 | 808 |
| 23s | 752 |
| 16s | 751 |
| relaxase | 510 |
| mash_neighbor | 509 |
| observed_host_range | 189 |

Top MOB-suite features by sample prevalence included:

```text
molecule_type_chromosome
16s-rRNA
23s-rRNA
MPF_T
ISKpn1
molecule_type_plasmid
IncFIB
MPF_F
conjugative
Enterobacterales
IncFII
MOBF
```

Interpretation caution: MOB-suite outputs provide plasmid reconstruction/typing context, not direct proof of gene transfer or phenotype. Plasmid localization should still be interpreted with contig/reconstruction evidence and assembly quality.

## Kleborate Feature Quality

Kleborate produced standardized organism-specific typing, K/O-locus, virulence, resistance-score, and AMR marker features.

Top Kleborate categories:

| Category | Rows |
| --- | ---: |
| amr_marker | 129 |
| virulence_score | 96 |
| sequence_type | 96 |
| resistance_class_count | 96 |
| k_locus | 96 |
| k_type | 96 |
| resistance_score | 96 |
| o_locus | 96 |
| o_type | 96 |
| wzi | 86 |
| yersiniabactin | 64 |
| aerobactin | 40 |

Top Kleborate features included:

```text
SHV-11
resistance_score_2
OL2alpha.1 / OL2alpha.2
aerobactin_iuc_1
O2alpha
yersiniabactin_ybt_9
virulence_score_1
rmpadc_rmp_1
ST11
KL64 / K64
SHV-1
```

## PanR2 Analysis Outputs

The optional runner outputs generated the same downstream PanR2-style analysis layers used for standard AMR/VFDB analysis:

```text
panr2_inputs/features/mobsuite.features.tsv
panr2_inputs/features/kleborate.features.tsv
panr2_inputs/features/all_features.tsv
panr2_inputs/feature_matrices/all_features_presence_absence.tsv
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/database_cooccurrence_summary.tsv
panr2_inputs/cross_database/feature_proximity.tsv
panr2_inputs/metadata_feature_analysis/top_features_by_database.tsv
panr2_inputs/metadata_feature_analysis/top_findings.tsv
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/cross_database_interpretation.html
```

Output sizes:

```text
feature_cooccurrence.tsv rows: 11,175 data rows
all_features_presence_absence.tsv rows: 100 samples
top_findings.tsv rows: 50 findings
```

## Bottlenecks Learned

The first 100-genome validation showed that MOB-suite and Kleborate were biologically useful and PanR2-compatible, but the original runner implementation was too serial/batch-heavy for routine larger validation:

```text
MOB-suite started: 2026-05-12T02:54:42Z
MOB-suite completed: 2026-05-12T04:19:58Z
Approximate MOB-suite stage time: 1h 25m
```

Kleborate also currently runs as one batch process:

```text
Kleborate-only validation runtime: 27.82m
Combined-run organism typing stage: 17.00m
```

Recommended engineering improvement before routine 100+ genome optional-runner validation:

```text
Split MOB-suite per genome into independent Nextflow tasks.
Split Kleborate into per-genome or chunked Nextflow tasks.
Merge raw outputs and PanR2 feature tables after parallel execution.
```

## Parallel Optional-Runner Revalidation

Date: 2026-05-12

Purpose: revalidate the same 100-genome biological MOB-suite/Kleborate workflow after adding per-genome parallel dispatch for both optional runners.

Command:

```bash
nextflow run main.nf \
  -profile conda,mamba,desktop_parallel,large \
  --local_samples validation_runs/optional_real_100_input \
  --outdir validation_runs/optional_real_100_klebsiella_mobsuite_kleborate_parallel \
  --sequence_qc_engine python \
  --qc_filter true \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_panr2_comprehensive false \
  --run_abricate false \
  --export_panr2_inputs true \
  --run_mobsuite true \
  --mobsuite_db validation_runs/mobsuite_db_test \
  --mobsuite_jobs 8 \
  --mobsuite_threads_per_sample 1 \
  --run_kleborate true \
  --kleborate_jobs 8 \
  --run_kaptive false \
  --run_ectyper false \
  --run_genomad false \
  --threads 8 \
  --capture_versions false
```

Status: PASS

Completed processes:

```text
SEQUENCE_QC
COMBINED_QC
MOBSUITE_ANALYSIS
ORGANISM_SPECIFIC_TYPING
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

Runtime comparison:

| Run | Total duration | CPU hours | MOB-suite stage | Kleborate/typing stage |
| --- | ---: | ---: | ---: | ---: |
| Original serial/batch optional run | 1h 44m 49s | 13.7 | ~1h 25m | ~17m combined / ~28m Kleborate-only |
| Parallel optional-runner run | 1h 00m 13s | 7.9 | 55.63m | 3.78m |

The parallel run preserved the same PanR2-standardized feature contract result:

```text
feature_files_checked=2
feature_rows=17490
databases_seen=kleborate,mobsuite
samples_seen=100
metadata_accessions=100
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Per-sample runner status:

| Module | Samples input | Samples processed | Samples failed | Status |
| --- | ---: | ---: | ---: | --- |
| MOB-suite | 100 | 100 | 0 | PASS |
| Kleborate | 100 | 100 | 0 | PASS |

Feature completeness audit:

| Database | Feature rows | Unique features | Samples with features | Status |
| --- | ---: | ---: | ---: | --- |
| `mobsuite` | 16,332 | 489 | 100 | PASS |
| `kleborate` | 1,158 | 263 | 96 | PASS |

PanR2 analysis parity with AMR/VFDB-style analysis was confirmed. The optional outputs generated:

```text
panr2_inputs/features/mobsuite.features.tsv
panr2_inputs/features/kleborate.features.tsv
panr2_inputs/features/all_features.tsv
panr2_inputs/feature_matrices/mobsuite_presence_absence.tsv
panr2_inputs/feature_matrices/kleborate_presence_absence.tsv
panr2_inputs/feature_matrices/all_features_presence_absence.tsv
panr2_inputs/metadata_feature_analysis/feature_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/database_burden_by_sample.tsv
panr2_inputs/metadata_feature_analysis/database_burden_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/category_burden_by_sample.tsv
panr2_inputs/metadata_feature_analysis/category_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/top_features_by_database.tsv
panr2_inputs/metadata_feature_analysis/top_findings.tsv
panr2_inputs/metadata_feature_analysis/prevalence_tables/
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/database_cooccurrence_summary.tsv
panr2_inputs/cross_database/feature_proximity.tsv
panr2_inputs/report/panr2_handoff_index.html
panr2_inputs/report/top_findings.html
panr2_inputs/report/metadata_quality_and_bias.html
panr2_inputs/report/database_burden_by_metadata.html
panr2_inputs/report/cross_database_interpretation.html
```

Selected analysis output sizes:

```text
all_features_presence_absence.tsv: 100 samples plus header
mobsuite_presence_absence.tsv: 100 samples plus header
kleborate_presence_absence.tsv: 100 samples plus header
feature_cooccurrence.tsv: 11,175 data rows
top_findings.tsv: 50 findings
top_features_by_database.tsv: 50 MOB-suite + 50 Kleborate rows
prevalence_tables: 12 metadata-stratified tables
```

Top MOB-suite features by sample prevalence included:

```text
molecule_type_chromosome
16s-rRNA
23s-rRNA
MPF_T
ISKpn1
molecule_type_plasmid
IncFIB
MPF_F
conjugative
Enterobacterales
```

Top Kleborate features by sample prevalence included:

```text
SHV-11
resistance_score_2
OL2alpha.1 / OL2alpha.2
aerobactin_iuc_1
resistance_score_0
O2alpha
yersiniabactin_ybt_9
virulence_score_1
ST11
```

Interpretation: MOB-suite and Kleborate now behave like other PanR2 feature databases once their raw outputs are converted into standardized feature tables. PanR2 can summarize them by prevalence, burden, metadata association, co-occurrence, and report pages. The remaining limitation is scientific, not architectural: MOB-suite output should be interpreted as plasmid reconstruction/typing context, and Kleborate is organism-specific to Klebsiella-like analyses.

## Scope Limit

This 100-genome real validation covers:

```text
MOB-suite runner
Kleborate runner
PanR2 feature contract export
PanR2-style analysis/report outputs
```

It does not validate full 100-genome runner mode for:

```text
geNomad
MobileElementFinder
DefenseFinder
Kaptive
ECTyper
SerotypeFinder
SCCmecFinder
```

Those should be validated with targeted small real runs and required databases before scaling to 100 genomes.
