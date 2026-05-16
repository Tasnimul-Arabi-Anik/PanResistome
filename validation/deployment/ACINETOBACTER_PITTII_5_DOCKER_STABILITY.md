# Acinetobacter pittii 5-genome Docker stability validation

Date: 2026-05-16

Purpose: validate the hardened remote-user Docker/GHCR comprehensive route after the
CheckM2, AMRFinderPlus, geNomAD, DefenseFinder, and MobileElementFinder stability
updates.

## Command

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --organism_max_records 5 \
  --outdir validation_runs/acinetobacter_pittii_5_docker_stability \
  -profile docker \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 true \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad true \
  --genomad_jobs 1 \
  --genomad_threads_per_sample 1 \
  --genomad_splits 8 \
  --genomad_sensitivity 3.0 \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder false \
  --panr2_run_defensefinder false \
  --threads 8 \
  --checkm2_threads 4 \
  --fetchm2_download_workers 2
```

Post-run checker:

```bash
python scripts/check_comprehensive_validation_outputs.py \
  --run-dir validation_runs/acinetobacter_pittii_5_docker_stability \
  --require-checkm2 \
  --require-genomad \
  --expect-zero-schema-errors
```

## Result

Status: PASS

Nextflow summary:

```text
Succeeded: 23
Duration: 1h 59m 19s
CPU hours: 4.9
Results: validation_runs/acinetobacter_pittii_5_docker_stability
```

The validation checker reported:

```text
PASS validation_runs/acinetobacter_pittii_5_docker_stability/Acinetobacter_pittii
```

## Key outputs

Downloaded/analyzed genomes:

```text
input records: 5
CheckM2 rows: 5
QC PASS rows: 5
```

Feature-contract summary:

```text
feature_files_checked=7
feature_rows=630
databases_seen=amr,amrfinderplus,mlst,prophage,vfdb
samples_seen=10
metadata_accessions=5
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Selected feature tables:

```text
all_features.tsv: 630 feature rows
prophage.features.tsv: 27 feature rows
amrfinderplus.features.tsv: 28 feature rows
```

geNomAD module status:

```text
status=PASS
samples_input=5
samples_processed=5
samples_failed=0
raw_tables_created=5
feature_rows_created=27
unique_features_created=135
```

Native runner merge audit:

```text
abricate: PASS, observed_raw_tables=15, feature_rows=494
integronfinder: PASS, observed_raw_tables=5, feature_rows=0
mlst: PASS, samples_processed=5, feature_rows=38
mobileelementfinder: SKIPPED, not enabled for this profile
```

Report pages generated included:

```text
report/index.html
report/panr2_handoff_index.html
report/lineage_context.html
report/diversity_summary.html
report/statistical_summary.html
report/amrfinder_abricate_concordance.html
report/cross_database_interpretation.html
report/database_setup_and_contract.html
```

## What this validates

- Docker/GHCR comprehensive route works for a fresh biological taxon run.
- FetchM2 taxon input generation works with `--taxon "Acinetobacter pittii"`.
- CheckM2 auto database download, model-load smoke test, and QC prediction completed.
- QUAST, ANI, Mash, AMRFinderPlus, geNomAD, native PanR2 feature runners, PanR2
  comprehensive reporting, PanR2 handoff export, and result collection completed.
- geNomAD produced positive biological calls, not only header-only outputs.
- Feature-contract validation was clean.
- DefenseFinder remained disabled, matching the recommended stable route.
- MobileElementFinder remained disabled for this validation and is still opt-in.

## Notes

The first run paid database setup costs for CheckM2 and geNomAD. Subsequent runs
using the same output/database directory and `-resume` should avoid most of that
download/setup time.

`samples_seen=10` and `metadata_accessions=5` is expected for this contract
summary because feature rows include multiple compatible sample identifiers while
metadata accessions count the five downloaded assemblies. The unmatched count was
zero.
