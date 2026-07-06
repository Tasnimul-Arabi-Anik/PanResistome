#!/usr/bin/env bash
set -euo pipefail

# Reusable Docker/GHCR launcher for large PanResistome important-report runs
# with shared databases. Defaults match the Dulab WGS workstation layout, but
# all paths can be overridden through environment variables.

PIPELINE_DIR="${PIPELINE_DIR:-$HOME/Work/Bioinformatics/wgs/04_workflows/PanResistome}"
DB_ROOT="${PANRESISTOME_DB_ROOT:-/mnt/storage/db}"
OUTDIR="${OUTDIR:-$HOME/Work/Bioinformatics/wgs/07_results/panresistome_important}"
WORKDIR="${WORKDIR:-$HOME/Work/Bioinformatics/wgs/work/panresistome}"
TAXON="${TAXON:-Acinetobacter pittii}"
THREADS="${THREADS:-32}"
FETCH_WORKERS="${FETCH_WORKERS:-4}"
AMRFINDER_JOBS="${AMRFINDER_JOBS:-24}"
CONTAINER_IMAGE="${CONTAINER_IMAGE:-ghcr.io/tasnimul-arabi-anik/panresistome:experimental}"

if [[ ! -d "$PIPELINE_DIR" ]]; then
  echo "ERROR: PIPELINE_DIR does not exist: $PIPELINE_DIR" >&2
  exit 1
fi
if [[ ! -d "$DB_ROOT" ]]; then
  echo "ERROR: PANRESISTOME_DB_ROOT does not exist: $DB_ROOT" >&2
  exit 1
fi

CHECKM2_DB="${CHECKM2_DB:-}"
if [[ -z "$CHECKM2_DB" ]]; then
  CHECKM2_DB="$(find "$DB_ROOT/checkm2" -path '*CheckM2_database/*.dmnd' -type f 2>/dev/null | sort | tail -n 1 || true)"
fi

GTDBTK_DATA_PATH="${GTDBTK_DATA_PATH:-}"
if [[ -z "$GTDBTK_DATA_PATH" ]]; then
  GTDBTK_DATA_PATH="$(find "$DB_ROOT/gtdbtk" -maxdepth 5 -type d -name 'release*' 2>/dev/null | sort | tail -n 1 || true)"
fi

GENOMAD_DB="${GENOMAD_DB:-}"
if [[ -z "$GENOMAD_DB" ]]; then
  GENOMAD_DB="$(find "$DB_ROOT/genomad" -maxdepth 1 -type d -name 'genomad_db*' 2>/dev/null | sort | tail -n 1 || true)"
fi

ABRICATE_DB="${ABRICATE_DB:-$DB_ROOT/abricate/legacy_abricate_env/current}"
MOBSUITE_DB_DIR="${MOBSUITE_DB_DIR:-$DB_ROOT/mobsuite}"

mkdir -p "$OUTDIR" "$WORKDIR" "$MOBSUITE_DB_DIR"

cd "$PIPELINE_DIR"

echo "PanResistome shared-database run"
echo "  pipeline:     $PIPELINE_DIR"
echo "  db root:      $DB_ROOT"
echo "  output:       $OUTDIR"
echo "  work:         $WORKDIR"
echo "  taxon:        $TAXON"
echo "  CheckM2 DB:   ${CHECKM2_DB:-not found; pass CHECKM2_DB or disable CheckM2}"
echo "  GTDB-Tk data: ${GTDBTK_DATA_PATH:-not found; GTDB-Tk should stay disabled}"
echo "  geNomAD DB:   ${GENOMAD_DB:-not found; geNomAD should stay disabled or auto-download elsewhere}"
echo "  ABRicate DB:  ${ABRICATE_DB:-not found}"

nextflow run main.nf \
  -work-dir "$WORKDIR" \
  --taxon "$TAXON" \
  --outdir "$OUTDIR" \
  -profile docker,large \
  --container_image "$CONTAINER_IMAGE" \
  --container_run_options "-v $DB_ROOT:$DB_ROOT" \
  --analysis_profile comprehensive \
  --output_mode important \
  --large_dataset true \
  --report_mode compact \
  --qc_filter true \
  --run_gtdbtk false \
  --gtdbtk_data_path "$GTDBTK_DATA_PATH" \
  --run_checkm2 false \
  --checkm2_db "$CHECKM2_DB" \
  --checkm2_db_dir "$DB_ROOT/checkm2" \
  --run_quast true \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus true \
  --run_genomad false \
  --genomad_db "$GENOMAD_DB" \
  --genomad_db_dir "$DB_ROOT/genomad" \
  --mobsuite_db_dir "$MOBSUITE_DB_DIR" \
  --db "$ABRICATE_DB" \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --panr2_run_mobileelementfinder true \
  --panr2_run_defensefinder false \
  --panr2_update_abricate_db false \
  --threads "$THREADS" \
  --fetchm2_download_workers "$FETCH_WORKERS" \
  --amrfinderplus_jobs "$AMRFINDER_JOBS" \
  --amrfinderplus_threads_per_sample 1 \
  "$@"
