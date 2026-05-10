# Klebsiella Optional Runner Biological Validation

Date: 2026-05-10

Purpose: run a small real biological optional-runner validation while keeping unrelated heavy/default modules disabled.

Input: two complete `Klebsiella pneumoniae` assemblies selected from the 300-record validation run:

```text
GCA_041085125.2_ASM4108512v2
GCA_041283555.2_ASM4128355v2
```

The input subset was staged locally under `validation_runs/optional_biological_input/` and is intentionally not committed because it contains copied genome FASTA files.

## Command

The run reused cached local environments for Kleborate and MOB-suite and disabled Conda environment creation:

```bash
PATH=$PWD/work/conda/env-dec1009326d0898c-b1741fff61178898a0bf90bfb494bd6f/bin:\
$PWD/work/conda/env-d29089f7ecf3a5bb-babb9a4070eff291b15019f5ec492fe1/bin:$PATH \
nextflow run main.nf \
  -profile test \
  --local_samples validation_runs/optional_biological_input \
  --outdir validation_runs/optional_runner_biological_klebsiella2_v4 \
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
  --run_kleborate true \
  --run_kaptive false \
  --run_ectyper false \
  --run_genomad false \
  --threads 2 \
  --capture_versions false
```

## Result

Status: PARTIAL PASS

Completed processes:

```text
SEQUENCE_QC
COMBINED_QC
MOBSUITE_ANALYSIS
ORGANISM_SPECIFIC_TYPING
EXPORT_PANR2_INPUTS
COLLECT_RESULTS
```

PanR2 feature outputs:

| Feature table | Rows | Status |
| --- | ---: | --- |
| `kleborate.features.tsv` | 25 | PASS |
| `mobsuite.features.tsv` | 0 | WARNING_EMPTY |

Schema validation:

```text
feature_files_checked=2
feature_rows=25
databases_seen=kleborate
samples_seen=2
metadata_accessions=2
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

## Kleborate Biological Calls

Kleborate v3.2.4 ran successfully with `--preset kpsc`. The PanR2 contract exporter converted real Kleborate output into standardized features including:

```text
ST147
ST592
KL64 / K64
KL57 / K57
OL2α.1 / O2α
OL3γ / O3γ
yersiniabactin_ybt_9
aerobactin_iuc_1_(truncated)
salmochelin_iro_1_(truncated)
rmpadc_rmp_1
wzi_wzi64
wzi_wzi206
SHV-11
DHA-1
SHV-26
```

Feature-contract validation was clean: zero unmatched, invalid, or duplicate feature rows.

## MOB-suite Result

MOB-suite was attempted on the same two assemblies, but the cached environment is broken:

```text
ModuleNotFoundError: No module named 'mob_suite.mob_recon'
```

The workflow now preserves per-sample MOB-suite logs under:

```text
mobsuite/raw/<sample>/mob_recon.stderr
```

and records the failure in:

```text
mobsuite/module_status.tsv
panr2_inputs/manifest/module_status_summary.tsv
```

The module status was:

```text
status=WARNING_EMPTY
samples_input=2
samples_processed=0
samples_failed=2
feature_rows_created=0
```

This means MOB-suite orchestration and failure reporting are working, but the local MOB-suite installation must be repaired before biological MOB-suite validation can pass.

## MOB-suite Environment Repair Follow-up

After the initial biological attempt, a clean disposable environment was tested. Installing `bioconda::mob_suite=3.1.9` exposed a working `mob_recon` module during installation, but the Conda post-link script ran `mob_init` and stalled during database initialization. This explains why interrupted or partially rolled-back cached environments can leave `mob_recon` entry points without the expected Python modules.

The repository MOB-suite environment was changed to install the heavy dependencies with Conda and `mob-suite==3.1.9` with pip. A disposable environment created from `envs/mobsuite.yaml` passed:

```text
import mob_suite.mob_recon
mob_recon 3.1.9
```

The two-genome run was repeated with that repaired environment. MOB-suite launched successfully and found `blastn`, `makeblastdb`, and `tblastn`, but it could not initialize/download its database in the current restricted environment:

```text
MOB-databases need to be initialized
ERROR: Something went wrong with database download or unpacking
```

The current status is therefore:

```text
MOB-suite executable/environment: fixed
PanResistome orchestration and logging: working
PanR2 header-only empty export: working
Biological MOB-suite feature validation: still blocked until a preinitialized MOB-suite database is supplied with --mobsuite_db or runtime database download succeeds
```

PanResistome now supports:

```text
--mobsuite_db /path/to/mob_suite/databases
```

which is passed to `mob_recon --database_directory`.

## Fixes Added From This Validation

1. Kleborate runner now uses `--preset kpsc`, required by Kleborate v3.
2. Kleborate output is collected from its output directory instead of treating `-o` as a single TSV file.
3. PanR2 contract export now converts real Kleborate output into feature rows for ST, virulence/resistance scores, K/O loci, siderophore loci, wzi, and AMR markers.
4. Kleborate sample mapping now prefers the `strain` column so assembly accessions match metadata.
5. MOB-suite now writes per-sample stdout/stderr and a `module_status.tsv` instead of hiding runtime failures behind an empty table.
6. MOB-suite environment creation now avoids the fragile Bioconda `mob_init` post-link by installing `mob-suite==3.1.9` with pip inside a Conda-managed environment.
7. MOB-suite runner mode now accepts `--mobsuite_db` for a preinitialized MOB-suite database directory.

## Next Biological Optional-Runner Targets

1. Supply a preinitialized MOB-suite database with `--mobsuite_db`, then rerun this same two-genome validation.
2. Add a Kaptive database path and run Kaptive on the same Klebsiella subset.
3. Add geNomad DB path and run geNomad on the same subset.
4. Use an E. coli two-genome subset for ECTyper validation.
