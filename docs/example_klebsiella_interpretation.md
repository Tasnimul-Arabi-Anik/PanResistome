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
panr2_inputs/metadata_feature_analysis/metadata_column_eligibility.tsv
panr2_inputs/metadata_feature_analysis/feature_metadata_associations.tsv
panr2_inputs/metadata_feature_analysis/database_burden_metadata_associations.tsv
```

Cross-database interpretation:

```text
panr2_inputs/cross_database/feature_cooccurrence.tsv
panr2_inputs/cross_database/amr_mge_context.tsv
panr2_inputs/cross_database/amr_plasmid_context.tsv
panr2_inputs/cross_database/amrfinder_abricate_concordance.tsv
panr2_inputs/cross_database/feature_proximity.tsv
```

## Interpretation Limits

- Sample-level co-occurrence means features occur in the same genome, not necessarily on the same DNA molecule.
- Same-contig evidence is stronger but still depends on assembly quality and contig fragmentation.
- Proximity evidence does not prove transfer, expression, phenotype, or plasmid localization.
- Metadata associations are exploratory unless group sizes, missingness, and BioProject/study bias are acceptable.
- AMRFinderPlus and ABRicate can disagree because they use different detection strategies and naming conventions.

## Practical Use

For manuscript-style interpretation, start with:

1. QC status and failed downloads.
2. Feature-contract validation.
3. Top AMR, VFDB, plasmid, integron, and MLST features.
4. Metadata column eligibility and top findings.
5. Same-contig/proximity context for AMR-plasmid, AMR-MGE, and AMR-integron questions.
