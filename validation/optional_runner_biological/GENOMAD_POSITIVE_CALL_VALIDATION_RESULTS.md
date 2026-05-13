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

## Follow-up 2-Genome Scale Check

Date: 2026-05-13

Purpose: confirm that the positive-call route works beyond one assembly before attempting larger geNomAD validation.

Command pattern:

```bash
nextflow run main.nf \
  --input /tmp/panresistome_genomad_positive/klebsiella_2.tsv \
  --outdir validation_runs/genomad_positive_klebsiella_2_ghcr \
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

Result:

```text
Nextflow processes succeeded: 11/11
runtime: 19m37s
CPU hours: 1.3
GENOMAD_PROPHAGE runtime: 14.52m
GENOMAD_PROPHAGE max RSS: 2.8 GB
GENOMAD_PROPHAGE max vmem: 11.4 GB
geNomAD database status: PASS, provided/cached DB
geNomAD samples input: 2
geNomAD samples processed: 2
geNomAD samples failed: 0
geNomAD collected region rows: 6
prophage.features.tsv rows: 6
all_features.tsv rows: 292
unmatched feature rows: 0
invalid feature rows: 0
duplicate feature rows: 0
```

Positive calls:

```text
GCF_041200225.2: 2 viral/prophage regions, 1 plasmid-like region
GCF_041200245.2: 2 viral/prophage regions, 1 plasmid-like region
```

Interpretation: the Docker/GHCR geNomAD path is now validated for a small positive 2-genome biological run with mounted cached DB and memory-safe MMseqs splitting. The next scale step should be 5-10 genomes on a machine where a 1-2 hour optional-module validation is acceptable.

## Follow-up 2-Genome Parallel geNomAD Check

Date: 2026-05-13

Purpose: test whether bounded per-genome geNomAD parallelism can reduce wall time without changing the PanR2 feature contract.

Command pattern:

```bash
nextflow run main.nf \
  --input /tmp/panresistome_genomad_positive/klebsiella_2.tsv \
  --outdir validation_runs/genomad_parallel_klebsiella_2_ghcr \
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
  --genomad_jobs 2 \
  --genomad_threads_per_sample 1 \
  --threads 2 \
  --fetchm2_download_workers 1
```

Fresh parallel result:

```text
Nextflow processes succeeded: 11/11
runtime: 14m45s
CPU hours: 0.9
GENOMAD_PROPHAGE runtime: 11.28m
GENOMAD_PROPHAGE max RSS: 5.6 GB
GENOMAD_PROPHAGE max vmem: 21.8 GB
geNomAD samples input: 2
geNomAD samples processed: 2
geNomAD samples failed: 0
geNomAD collected region rows: 6
prophage.features.tsv rows: 6
all_features.tsv rows: 292
unmatched feature rows: 0
invalid feature rows: 0
duplicate feature rows: 0
```

Per-sample audit:

```text
GCF_041200225.2_ASM4120022v2_genomic: PASS, 677 seconds
GCF_041200245.2_ASM4120024v2_genomic: PASS, 676 seconds
```

Resume behavior:

```text
--genomad_reuse_existing true
genomad_sample_status.tsv records PASS_REUSED for existing non-empty per-sample summaries
GENOMAD_PROPHAGE resume runtime: 0.10s
```

Interpretation: parallelism worked, but it is not free. Two concurrent one-thread geNomAD jobs reduced total wall time from 19m37s to 14m45s and reduced `GENOMAD_PROPHAGE` from 14.52m to 11.28m, while roughly doubling peak RSS from 2.8 GB to 5.6 GB and peak vmem from 11.4 GB to 21.8 GB. The practical recommendation is:

```text
16 GB desktop: --genomad_jobs 1 or 2, --genomad_threads_per_sample 1
32 GB workstation: --genomad_jobs 2-4, --genomad_threads_per_sample 1
HPC/container node: increase --genomad_jobs gradually after checking trace memory
```

PanResistome now writes `prophage/genomad_sample_status.tsv` so users can see which per-genome geNomAD jobs passed, failed, or were reused.

## Follow-up 5-Genome geNomAD Scale Check

Date: 2026-05-13

Purpose: validate that the bounded parallel geNomAD route can scale beyond two genomes while keeping the PanR2 feature contract clean.

Command pattern:

```bash
nextflow run main.nf \
  --input /tmp/panresistome_genomad_positive/klebsiella_5.tsv \
  --outdir validation_runs/genomad_parallel_klebsiella_5_ghcr \
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
  --genomad_jobs 2 \
  --genomad_threads_per_sample 1 \
  --threads 2 \
  --fetchm2_download_workers 1
```

Result:

```text
Nextflow processes succeeded: 11/11
runtime: 29m39s
CPU hours: 1.7
GENOMAD_PROPHAGE runtime: 24.32m
GENOMAD_PROPHAGE max RSS: 5.6 GB
GENOMAD_PROPHAGE max vmem: 23.9 GB
geNomAD samples input: 5
geNomAD samples processed: 5
geNomAD samples failed: 0
geNomAD collected region rows: 23
prophage.features.tsv rows: 23
all_features.tsv rows: 601
unmatched feature rows: 0
invalid feature rows: 0
duplicate feature rows: 0
```

Per-sample audit:

```text
GCF_041200225.2_ASM4120022v2_genomic: PASS, 514 seconds
GCF_041200245.2_ASM4120024v2_genomic: PASS, 514 seconds
GCF_050247555.1_ASM5024755v1_genomic: PASS, 495 seconds
GCF_048279315.2_ASM4827931v2_genomic: PASS, 503 seconds
GCF_050269205.1_ASM5026920v1_genomic: PASS, 449 seconds
```

Interpretation: this validates the 5-genome Docker/GHCR geNomAD scale route with a mounted cached database, two concurrent one-thread jobs, low-memory MMseqs splitting, and clean PanR2 handoff output. The memory profile stayed similar to the two-genome parallel run because concurrency remained capped at two jobs. Larger geNomAD runs should continue increasing genome count before increasing `--genomad_jobs`.
