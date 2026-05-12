# Klebsiella 5-Genome MobileElementFinder Biological Validation

Date: 2026-05-12

Purpose: test whether the opt-in MobileElementFinder runner can generate real biological outputs and feed the same PanR2 feature-contract/reporting layer used for AMR/VFDB-style analyses.

Input:

```text
Source validation input: validation_runs/optional_real_100_input/Klebsiella_pneumoniae
Subset: 5 Klebsiella pneumoniae assemblies
Metadata: matching FetchM2 ncbi_clean.csv rows
```

This validation intentionally kept heavy unrelated modules off so the MobileElementFinder path could be tested quickly.

## Direct Runner Smoke

The cached PanR2 comprehensive environment with MobileElementFinder was used to run `panr2.runners.run_mobileelementfinder` directly on the 5 FASTAs.

Result:

```text
raw MobileElementFinder CSV files: 5/5
raw MobileElementFinder result text files: 5/5
raw MobileElementFinder MGE FASTA files: 5/5
PanR2-style MobileElementFinder result rows: 248 including header
summary rows: 6 including header
```

Then the PanR2 contract exporter was run on the generated tables.

Schema validation:

```text
feature_files_checked=1
feature_rows=245
databases_seen=mobileelementfinder
samples_seen=5
metadata_accessions=5
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Interpretation: MobileElementFinder produced biological MGE features, including insertion sequences and ICE calls, and those rows cleanly exported into `mobileelementfinder.features.tsv` and `all_features.tsv`.

## Nextflow Integration Validation

Command:

```bash
nextflow run main.nf \
  --local_samples validation_runs/broad_optional_small_input \
  --outdir validation_runs/mobileelementfinder_small_nextflow \
  -profile conda,mamba \
  --analysis_profile amr_basic \
  --sequence_qc_engine python \
  --qc_filter true \
  --stop_after_qc false \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode serial \
  --panr2_run_mobileelementfinder true \
  --panr2_mobileelementfinder_allow_failure false \
  --panr2_update_abricate_db false \
  --threads 2 \
  --capture_versions false
```

Result:

```text
Nextflow processes succeeded: 6/6
Duration: 3m 3s
CPU hours: 0.1
Samples processed: 5
```

Native runner audit:

```text
abricate: PASS, raw_tables=1, feature_rows=33, unique_features=29
mobileelementfinder: PASS, raw_tables=5, feature_rows=248, unique_features=95
integronfinder: SKIPPED
mlst: SKIPPED
```

Final PanR2 schema validation:

```text
feature_files_checked=2
feature_rows=753
databases_seen=amr,mobileelementfinder
samples_seen=5
metadata_accessions=5
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Final standardized feature files:

```text
panr2_inputs/features/amr.features.tsv
panr2_inputs/features/mobileelementfinder.features.tsv
panr2_inputs/features/all_features.tsv
```

PanR2 also generated MobileElementFinder-specific analysis/report outputs under:

```text
mobileelementfinder/analysis/
mobileelementfinder/figures/
panr2_inputs/cross_database/
panr2_inputs/metadata_feature_analysis/
panr2_inputs/report/
```

## Conclusion

MobileElementFinder is now biologically validated on a small real Klebsiella subset and can generate PanR2-compatible standardized features and report outputs.

Recommended support level after this validation:

```text
PanR2 table analysis: PASS
PanResistome opt-in runner: PASS on 5 real Klebsiella genomes
Default comprehensive mode: keep disabled
Reason to keep opt-in: runtime and upstream parser behavior should be tested on more species before default use
```

MobileElementFinder should remain opt-in via:

```bash
--panr2_run_mobileelementfinder true
```

For cautious runs, keep the default nonfatal behavior:

```bash
--panr2_mobileelementfinder_allow_failure true
```
