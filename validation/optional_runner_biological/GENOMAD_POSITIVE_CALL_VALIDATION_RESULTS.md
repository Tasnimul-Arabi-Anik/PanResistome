# geNomAD Positive-Call Biological Validation

Date: 2026-05-13

Purpose: validate that the opt-in geNomAD runner can produce real biological viral/prophage calls, export clean PanR2-compatible `prophage.features.tsv`, and avoid the low-memory MMseqs failure observed in earlier complete-genome attempts.

## Result

Status: PASS with low-memory geNomAD settings.

Validated command pattern:

```bash
nextflow run main.nf \
  --input /tmp/panresistome_genomad_positive/klebsiella_1.tsv \
  --outdir validation_runs/genomad_positive_klebsiella_1_ghcr_final \
  -profile docker,large,genomad_host \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --container_run_options '-v /tmp/panresistome_genomad_db:/tmp/panresistome_genomad_db' \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --run_abricate false \
  --panr2_native_feature_runners false \
  --run_genomad true \
  --genomad_use_host_env true \
  --genomad_db /tmp/panresistome_genomad_db/genomad_db \
  --genomad_splits 8 \
  --genomad_sensitivity 3.0 \
  --threads 1 \
  --fetchm2_download_workers 1
```

## Validation Summary

```text
Nextflow processes succeeded: 11/11
runtime: 10m13s
CPU hours: 0.5
GENOMAD_PROPHAGE runtime: 7.63m
GENOMAD_PROPHAGE max RSS: 2.8 GB
GENOMAD_PROPHAGE max vmem: 11.4 GB
geNomAD database status: PASS, provided/cached DB
geNomAD samples input: 1
geNomAD samples processed: 1
geNomAD samples failed: 0
geNomAD collected region rows: 3
prophage.features.tsv rows: 3
all_features.tsv rows: 146
unmatched feature rows: 0
invalid feature rows: 0
duplicate feature rows: 0
```

The positive geNomAD calls were:

```text
plasmid_region:NZ_CP184712.1
viral_region:NZ_CP165751.1|provirus_4617910_4658925
viral_region:NZ_CP165751.1|provirus_2756504_2811596
```

## What Changed During Validation

The first complete-genome geNomAD attempts failed inside MMseqs `prefilter` with `SIGKILL`, even at `--threads 1`. Inspecting `genomad end-to-end --help` showed that geNomAD supports:

```text
--splits       split MMseqs searches to reduce memory usage
--sensitivity  marker-search sensitivity
```

PanResistome now exposes these as:

```text
--genomad_splits
--genomad_sensitivity
```

Using `--genomad_splits 8 --genomad_sensitivity 3.0 --threads 1` allowed the same focused biological validation to complete.

The generic optional-table collector was also tightened for geNomAD so that PanR2 consumes geNomAD summary outputs, not low-level MMseqs/gene/intermediate tables. The exported geNomAD feature table now reports region-level rows with `viral_region` and `plasmid_region` feature categories.

## Interpretation

This validates the geNomAD runner path for a small real bacterial genome with positive viral/prophage calls through Docker/GHCR and a mounted cached geNomAD database.

Recommended geNomAD command additions for memory-constrained desktops:

```text
--threads 1
--genomad_splits 8
--genomad_sensitivity 3.0
```

geNomAD should remain opt-in. The next validation step is a 2-10 genome run using the same split settings, followed by a larger validation only on a machine with enough memory and disk.

## Remaining Limits

This does not prove:

```text
fresh geNomAD database auto-download on every network
large 100-genome geNomAD scalability
Apptainer execution
DefenseFinder runner execution
```

Those remain separate validation targets.
