# Interpreting The Klebsiella Validation Output

The `Klebsiella pneumoniae` validation is intended to show what the comprehensive workflow produces on a biologically richer organism than the Delftia fresh-user validation.

## Start Here

Open:

```text
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/report/index.html
validation_runs/klebsiella_pneumoniae_parallel/Klebsiella_pneumoniae/panr2_inputs/report/panr2_handoff_index.html
```

The committed validation report is:

```text
validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md
```

## What Was Detected

The validation produced 12,838 standardized PanR2 feature rows from 99 analyzed genomes:

```text
AMR: 1234
AMRFinderPlus: 1350
VFDB: 7821
PlasmidFinder: 393
IntegronFinder: 418
MLST: 1622
```

This confirms that the feature space is large enough to test cross-database analysis and metadata-aware reporting.

## Key Tables

Feature tables:

```text
panr2_inputs/features/all_features.tsv
panr2_inputs/features/amr.features.tsv
panr2_inputs/features/amrfinderplus.features.tsv
panr2_inputs/features/vfdb.features.tsv
panr2_inputs/features/plasmidfinder.features.tsv
panr2_inputs/features/integronfinder.features.tsv
panr2_inputs/features/mlst.features.tsv
```

Audits:

```text
panr2_inputs/manifest/schema_validation_summary.txt
panr2_inputs/manifest/feature_completeness_audit.tsv
panr2_inputs/manifest/database_setup_status.tsv
panr2_inputs/manifest/native_runner_merge_audit.tsv
```

Metadata interpretation:

```text
panr2_inputs/metadata_feature_analysis/top_findings.tsv
panr2_inputs/metadata_feature_analysis/metadata_usability_summary.tsv
panr2_inputs/metadata_feature_analysis/metadata_column_eligibility.tsv
panr2_inputs/metadata_feature_analysis/feature_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/database_burden_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/bioproject_bias_report.tsv
```

Start with `metadata_usability_summary.tsv` before interpreting top findings. It lists missingness, cardinality, largest-group dominance, and whether a metadata column is recommended for comparative analysis.

`top_findings.tsv` includes interpretation-safety columns:

```text
supporting_samples
largest_bioproject
largest_bioproject_fraction
warning_flags
interpretation_label
```

Treat findings labeled `bioproject_bias_warning`, `low_sample_warning`, or `sparse_metadata_warning` as exploratory. A high `largest_bioproject_fraction` means one study may be driving the signal.

Cross-database interpretation:

```text
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/amr_mge_context.tsv
panr2_inputs/cross_database/amr_plasmid_context.tsv
panr2_inputs/cross_database/amrfinder_abricate_concordance.tsv
panr2_inputs/cross_database/feature_proximity.tsv
panr2_inputs/cross_database/feature_proximity_all.tsv
```

`amrfinder_abricate_concordance.tsv` classifies AMR calls as:

```text
called_by_both
abricate_only
amrfinderplus_only
same_symbol_different_samples
possible_class_match
```

`possible_class_match` means both tools reported an AMR feature in the same sample with a shared class label, but not the same normalized gene symbol.

`feature_proximity.tsv` and the AMR context tables include:

```text
evidence_level
interpretation_warning
```

Use `evidence_level` to distinguish same-contig coordinate context from weaker same-genome co-occurrence. Same-contig or within-10-kb evidence is stronger than sample-level co-occurrence, but it still does not prove transfer, expression, phenotype, or plasmid localization.

In large-dataset mode, `feature_proximity.tsv` is capped for report readability. Use `feature_proximity_all.tsv` when you need the complete proximity evidence table.

HTML interpretation pages:

```text
panr2_inputs/report/top_findings.html
panr2_inputs/report/metadata_quality_and_bias.html
panr2_inputs/report/bioproject_bias.html
panr2_inputs/report/cross_database_interpretation.html
panr2_inputs/report/amrfinder_abricate_concordance.html
```

## Interpretation Limits

- Sample-level co-occurrence means features occur in the same genome, not necessarily on the same DNA molecule.
- Same-contig evidence is stronger but still depends on assembly quality and contig fragmentation.
- Proximity evidence does not prove transfer, expression, phenotype, or plasmid localization.
- Metadata associations are exploratory unless group sizes, missingness, and BioProject/study bias are acceptable.
- AMRFinderPlus and ABRicate can disagree because they use different detection strategies and naming conventions.
- Concordance rows based on normalized symbols or class labels are screening summaries; inspect raw ABRicate and AMRFinderPlus tables before resolving disagreements.

## Practical Use

For manuscript-style interpretation, start with:

1. QC status and failed downloads.
2. Feature-contract validation.
3. Top AMR, VFDB, plasmid, integron, and MLST features.
4. Metadata usability, BioProject bias, and top findings.
5. AMRFinderPlus-vs-ABRicate concordance.
6. Same-contig/proximity context for AMR-plasmid, AMR-MGE, and AMR-integron questions.
