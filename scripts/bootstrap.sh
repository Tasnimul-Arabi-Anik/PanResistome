#!/usr/bin/env bash
set -euo pipefail

CHECKM2_DB=""
CHECKM2_DB_DIR=""
ABRICATE_DB="./db"
GTDBTK_DATA=""
DOWNLOAD_CHECKM2_DB=false
RUN_ENV_CHECK=false

usage() {
    cat <<'EOF'
PanResistome bootstrap/preflight

Usage:
  scripts/bootstrap.sh [options]

Options:
  --checkm2-db PATH          Existing CheckM2 DIAMOND database file
  --download-checkm2-db DIR  Download and extract the CheckM2 database under DIR
  --abricate-db DIR          ABRicate database directory [default: ./db]
  --gtdbtk-data DIR          Optional GTDB-Tk reference data directory
  --check-envs               Ask Nextflow/Conda to resolve pipeline environments
  -h, --help                 Show this help

Example:
  scripts/bootstrap.sh --download-checkm2-db "$HOME/databases" --abricate-db ./db
EOF
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

warn() {
    echo "WARN: $*" >&2
}

have() {
    command -v "$1" >/dev/null 2>&1
}

abs_path() {
    local path="$1"
    if [[ -d "$path" ]]; then
        (cd "$path" && pwd)
    else
        local dir
        dir="$(dirname "$path")"
        local base
        base="$(basename "$path")"
        (cd "$dir" && printf "%s/%s\n" "$(pwd)" "$base")
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --checkm2-db)
            CHECKM2_DB="${2:-}"
            shift 2
            ;;
        --download-checkm2-db)
            CHECKM2_DB_DIR="${2:-}"
            DOWNLOAD_CHECKM2_DB=true
            shift 2
            ;;
        --abricate-db)
            ABRICATE_DB="${2:-}"
            shift 2
            ;;
        --gtdbtk-data)
            GTDBTK_DATA="${2:-}"
            shift 2
            ;;
        --check-envs)
            RUN_ENV_CHECK=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            fail "Unknown option: $1"
            ;;
    esac
done

[[ -f main.nf ]] || fail "Run this script from the PanResistome repository root."

echo "== PanResistome preflight =="

have git || fail "git is not installed or not on PATH."
have nextflow || fail "nextflow is not installed or not on PATH."
have conda || fail "conda/mamba compatible conda command is not installed or not on PATH."
PROFILE="conda"
if have mamba; then
    PROFILE="conda,mamba"
    echo "OK: mamba found; recommended Nextflow profile: $PROFILE"
else
    warn "mamba was not found. Plain conda is supported, but optional heavy-tool environments can solve slowly; install mamba/Miniforge for the recommended profile."
fi

if ! java -version >/dev/null 2>&1; then
    fail "Java is not installed or not on PATH; Nextflow requires Java."
fi

echo "OK: required commands found"

for env_file in envs/fetchm.yaml envs/checkm2.yaml envs/abricate.yaml envs/amrfinderplus.yaml envs/ani.yaml envs/quast.yaml envs/mash.yaml envs/panr2_comprehensive.yaml envs/mobsuite.yaml envs/genomad.yaml envs/organism_typing.yaml; do
    [[ -f "$env_file" ]] || fail "Missing environment file: $env_file"
done

if [[ ! -d "$ABRICATE_DB" ]]; then
    fail "ABRicate database directory not found: $ABRICATE_DB"
fi
if ! find "$ABRICATE_DB" -mindepth 2 -type f -print -quit | grep -q .; then
    fail "ABRicate database directory exists but does not contain database files: $ABRICATE_DB"
fi
ABRICATE_DB="$(abs_path "$ABRICATE_DB")"
echo "OK: ABRicate database directory: $ABRICATE_DB"

if [[ "$DOWNLOAD_CHECKM2_DB" == true ]]; then
    [[ -n "$CHECKM2_DB_DIR" ]] || fail "--download-checkm2-db requires a directory"
    mkdir -p "$CHECKM2_DB_DIR"
    CHECKM2_DB_DIR="$(abs_path "$CHECKM2_DB_DIR")"
    archive="$CHECKM2_DB_DIR/checkm2_database.tar.gz"
    CHECKM2_DB="$CHECKM2_DB_DIR/CheckM2_database/uniref100.KO.1.dmnd"

    if [[ ! -s "$CHECKM2_DB" ]]; then
        have curl || fail "curl is required for --download-checkm2-db."
        echo "Downloading CheckM2 database to $archive"
        curl -L -C - -o "$archive" "https://zenodo.org/api/records/14897628/files/checkm2_database.tar.gz/content"
        echo "Extracting CheckM2 database under $CHECKM2_DB_DIR"
        tar -xzf "$archive" -C "$CHECKM2_DB_DIR"
    fi
fi

if [[ -n "$CHECKM2_DB" ]]; then
    [[ -s "$CHECKM2_DB" ]] || fail "CheckM2 database file not found or empty: $CHECKM2_DB"
    CHECKM2_DB="$(abs_path "$CHECKM2_DB")"
    echo "OK: CheckM2 database: $CHECKM2_DB"
else
    warn "No CheckM2 database provided. The pipeline now attempts automatic CheckM2 database download by default; use --checkm2_auto_download_db false to disable it."
fi

if [[ -n "$GTDBTK_DATA" ]]; then
    [[ -d "$GTDBTK_DATA" ]] || fail "GTDB-Tk data directory not found: $GTDBTK_DATA"
    GTDBTK_DATA="$(abs_path "$GTDBTK_DATA")"
    echo "OK: GTDB-Tk data directory: $GTDBTK_DATA"
else
    echo "OK: GTDB-Tk data not required unless --run_gtdbtk true"
fi

echo "Checking Nextflow script syntax"
nextflow run main.nf --help --run_gtdbtk false --stop_after_qc true >/dev/null
echo "OK: Nextflow script syntax"

if [[ "$RUN_ENV_CHECK" == true ]]; then
    echo "Resolving Conda environments through Nextflow"
    nextflow run main.nf --help -profile "$PROFILE" --run_gtdbtk false --stop_after_qc true >/dev/null
    echo "OK: $PROFILE profile resolved"
fi

echo
echo "Preflight complete."
echo
echo "Recommended QC validation command:"
printf "nextflow run main.nf --input test_small.tsv --outdir results_small -profile %s --stop_after_qc true --run_gtdbtk false --threads 8 --db %q" "$PROFILE" "$ABRICATE_DB"
if [[ -n "$CHECKM2_DB" ]]; then
    printf " --checkm2_db %q" "$CHECKM2_DB"
fi
echo
echo
echo "Recommended full command after QC validation:"
printf "nextflow run main.nf --input test_small.tsv --outdir results_small -profile %s --run_gtdbtk false --qc_filter true --threads 8 --db %q" "$PROFILE" "$ABRICATE_DB"
if [[ -n "$CHECKM2_DB" ]]; then
    printf " --checkm2_db %q" "$CHECKM2_DB"
fi
echo
echo
echo "Optional comparative-genomics command:"
printf "nextflow run main.nf --input test_small.tsv --outdir results_comparative -profile %s --run_gtdbtk false --run_quast true --run_ani true --run_mash true --qc_filter true --threads 8 --db %q" "$PROFILE" "$ABRICATE_DB"
if [[ -n "$CHECKM2_DB" ]]; then
    printf " --checkm2_db %q" "$CHECKM2_DB"
fi
echo
echo
echo "Recommended standard comprehensive validation command:"
printf "nextflow run main.nf --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv --outdir validation_runs/delftia_fresh -profile %s --analysis_profile comprehensive --qc_filter true --run_gtdbtk false --run_quast true --run_ani true --run_mash true --run_amrfinderplus true --threads 4 --fetchm2_download_workers 2 --db %q" "$PROFILE" "$ABRICATE_DB"
if [[ -n "$CHECKM2_DB" ]]; then
    printf " --checkm2_db %q" "$CHECKM2_DB"
fi
echo
echo "See docs/remote_user_validation.md for expected outputs and release-passing criteria."
