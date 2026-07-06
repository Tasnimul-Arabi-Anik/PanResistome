# Future Agent Runbook

This runbook is for future AI agents or maintainers working on PanResistome.
It summarizes what the pipeline does, how to run it, and where interpretation
can go wrong.

## What PanResistome Does

PanResistome analyzes bacterial genome assemblies at population scale. The
normal comprehensive route can:

- fetch metadata and assemblies through FetchM2/NCBI inputs;
- run sequence QC, optional CheckM2, optional QUAST, optional ANI/skani, and
  optional Mash;
- annotate AMR, virulence, plasmid replicons, integrons, MLST, and optional
  modules such as AMRFinderPlus, geNomAD, MOB-suite, MobileElementFinder,
  Kaptive, Kleborate, and ISfinder-compatible BLAST;
- export standardized `panr2_inputs/features/*.features.tsv` files;
- generate `basic/`, `important/`, or complete result bundles;
- build `important/results.html`, a portable dashboard report with figures,
  interpretation tables, warnings, and downloads.

The architecture is:

```text
PanResistome = heavy execution, database setup, QC, provenance, feature export
PanR2        = standardized comparative analysis and report generation
```

## Recommended Output Mode

For large real runs, use:

```text
--output_mode important
--large_dataset true
--report_mode compact
```

This keeps the user-facing output manageable while preserving complete TSVs and
the full `panr2_inputs/` handoff where needed.

## Shared Database Runs

Do not let every run download large databases into its own output directory.
On shared workstations, use a stable database root such as:

```text
/mnt/storage/db
```

Use explicit parameters:

```bash
--checkm2_db /mnt/storage/db/checkm2/v1.1.0/extracted/CheckM2_database/uniref100.KO.1.dmnd
--checkm2_db_dir /mnt/storage/db/checkm2/v1.1.0
--gtdbtk_data_path /mnt/storage/db/gtdbtk/r232/extracted/release232
--genomad_db /mnt/storage/db/genomad/<genomad_db_version>
--genomad_db_dir /mnt/storage/db/genomad
--mobsuite_db_dir /mnt/storage/db/mobsuite
--db /mnt/storage/db/abricate/legacy_abricate_env/current
--container_run_options "-v /mnt/storage/db:/mnt/storage/db"
```

`--container_run_options` is required for Docker/Apptainer-style runs when the
database path is outside the repository or output directory.

## Dulab WGS Workspace Pattern

On the Dulab workstation, use:

```text
Pipeline clone:  ~/Work/Bioinformatics/wgs/04_workflows/PanResistome
Runs:            ~/Work/Bioinformatics/wgs/06_runs/<project>/
Results:         ~/Work/Bioinformatics/wgs/07_results/<project>/
Reports:         ~/Work/Bioinformatics/wgs/08_reports/<project>/
Shared DB root:  /mnt/storage/db
```

The reusable shared-database wrapper is:

```bash
templates/run_panresistome_important_shared_db.sh
```

Copy or symlink it into the WGS workspace `09_scripts/` directory and edit only
project-specific input/output values.

## Large Acinetobacter-Scale Run Pattern

For a large `Acinetobacter pittii` important report on a workstation with ample
RAM, start with GTDB-Tk and CheckM2 off unless a QC/taxonomy audit specifically
requires them:

```bash
nextflow run main.nf \
  --taxon "Acinetobacter pittii" \
  --outdir /path/to/results_acinetobacter_pittii_all_important \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --container_run_options "-v /mnt/storage/db:/mnt/storage/db" \
  --analysis_profile comprehensive \
  --output_mode important \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast true \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder true \
  --threads 32 \
  --amrfinderplus_jobs 24 \
  --amrfinderplus_threads_per_sample 1 \
  --fetchm2_download_workers 4
```

Enable geNomAD only when runtime and memory are acceptable:

```bash
--run_genomad true \
--genomad_db /mnt/storage/db/genomad/<genomad_db_version> \
--genomad_jobs 1 \
--genomad_threads_per_sample 1 \
--genomad_splits 8 \
--genomad_sensitivity 3.0
```

## Important Cautions

- Important report summaries are exploratory unless a section says otherwise.
- Geographic patterns reflect the analyzed dataset, not global prevalence.
- Temporal trends can be invalid if placeholder dates are not cleaned.
- Co-occurrence does not prove physical linkage; proximity/same-contig evidence
  is stronger but still not proof of transfer or expression.
- Lineage context does not replace phylogenetic analysis.
- Notable-genome scores are research prioritization only, not clinical risk.
- Small metadata groups, BioProject dominance, and clone/ST dominance can drive
  apparent associations.
- GTDB-Tk and CheckM2 are useful but resource-heavy; skipping them must be
  explicitly stated in interpretation.

## Regenerating Only The Important Report

If heavy annotation has already completed, regenerate report-facing outputs from
an existing sample directory:

```bash
python scripts/export_panr2_inputs.py \
  --sample-dir /path/to/<outdir>/<sample_dir> \
  --large-dataset \
  --report-mode compact \
  --output-mode important \
  --max-features-network 150 \
  --max-features-heatmap 150 \
  --skip-heavy-interactive-plots
```

This is the preferred route while polishing visualizations and write-ups because
it avoids rerunning expensive annotation.

## Report QA

After regeneration:

```bash
python scripts/check_important_report_outputs.py /path/to/<outdir>/<sample_dir>
python scripts/check_important_report_visual_layout.py /path/to/<outdir>/<sample_dir> --browser skip
```

For visual changes, also run browser screenshot QA when possible:

```bash
python scripts/check_important_report_visual_layout.py \
  /path/to/<outdir>/<sample_dir> \
  --browser require \
  --screenshots-dir /tmp/panresistome_visual_qa \
  --section-screenshots
```

## What Not To Change Casually

- Do not delete `panr2_inputs/` outputs.
- Do not remove complete TSVs to make the report smaller.
- Do not add new default databases without validation and documentation.
- Do not hide optional-module failures; write module status and warnings.
- Do not reclassify warning-heavy figures as publication-ready only because
  PNG/SVG/PDF/data files exist.
- Do not use raw warning-row counts as the headline warning number without
  grouping into interpretable warning categories.
