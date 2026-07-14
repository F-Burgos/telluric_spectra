#!/usr/bin/env bash

set -u
set -o pipefail

warn() {
  echo "[WARNING] $*" >&2
}

info() {
  echo "[INFO] $*"
}

usage() {
  cat <<'EOF'
Usage:
  run_pipeline.sh [options]

Runs the telluric pipeline in order:
  1) phase1/raw_table.py
  2) phase1/cleantable.py
  3) phase2/smoke_object_group_telluric_v2.py

Options:
  --science PATH          Science root path (default: /mnt/disco_datos/data/HARPS/science)
  --calib PATH            Calibration root path (default: /mnt/disco_datos/data/HARPS/calib)
  --phase1-dir PATH       Phase1 directory (default: <script_dir>/phase1)
  --phase2-dir PATH       Phase2 directory (default: <script_dir>/phase2)
  --work-dir PATH         Stage 1 Parquet directory (default: <phase1_dir>)
  --template-config PATH  Template config for smoke script (default: <phase2_dir>/smoke_config.ini)
  --telluric-script PATH  Path to telluric_spectra.py (default: <phase2_dir>/telluric_spectra.py)
  --tables-path PATH      Path for tables output (default: <phase2_dir>/spectra_results)
  --output PATH           Output base for object groups (default: <current_dir>/spectra)
  --star NAME             Process only one star (default: all)
  --products LIST         FITS products scanned in phase 1
                          (default: e2ds_A,ccf_A,s1d_A,bis_A)
  --radius-near ARCSEC    Radius for cleantable.py (default: 50.0)
  --plot-orders LIST      Optional telluric preview order(s) to plot; examples:
                          63, 10,20, 50-55, all. Default: no plots
  --python BIN            Python executable (default: python3)
  --no-fresh              Do not delete previous phase1 parquet outputs before running
  -h, --help              Show this help message
EOF
}

run_step() {
  local step_name="$1"
  shift

  info "Starting ${step_name}"
  echo "[CMD] $*"

  "$@"
  local status=$?

  if [ "$status" -ne 0 ]; then
    warn "${step_name} failed with exit code ${status}."
    return "$status"
  fi

  info "Completed ${step_name}"
  return 0
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$(pwd -P)"

SCIENCE_PATH="/mnt/disco_datos/data/HARPS/science"
CALIB_PATH="/mnt/disco_datos/data/HARPS/calib"
PHASE1_DIR="${SCRIPT_DIR}/phase1"
PHASE2_DIR="${SCRIPT_DIR}/phase2"
WORK_DIR=""
TEMPLATE_CONFIG=""
TELLURIC_SCRIPT=""
TABLES_PATH=""
OUTPUT_BASE=""
STAR="all"
PRODUCTS="e2ds_A,ccf_A,s1d_A,bis_A"
RADIUS_NEAR="50.0"
PLOT_ORDERS=""
PYTHON_BIN="python3"
FRESH_RUN="1"

while [ "$#" -gt 0 ]; do
  case "$1" in
    --science)
      SCIENCE_PATH="$2"
      shift 2
      ;;
    --calib)
      CALIB_PATH="$2"
      shift 2
      ;;
    --phase1-dir)
      PHASE1_DIR="$2"
      shift 2
      ;;
    --phase2-dir)
      PHASE2_DIR="$2"
      shift 2
      ;;
    --work-dir)
      WORK_DIR="$2"
      shift 2
      ;;
    --template-config)
      TEMPLATE_CONFIG="$2"
      shift 2
      ;;
    --telluric-script)
      TELLURIC_SCRIPT="$2"
      shift 2
      ;;
    --tables-path)
      TABLES_PATH="$2"
      shift 2
      ;;
    --output)
      OUTPUT_BASE="$2"
      shift 2
      ;;
    --star)
      STAR="$2"
      shift 2
      ;;
    --products)
      PRODUCTS="$2"
      shift 2
      ;;
    --radius-near)
      RADIUS_NEAR="$2"
      shift 2
      ;;
    --plot-orders)
      PLOT_ORDERS="$2"
      shift 2
      ;;
    --python)
      PYTHON_BIN="$2"
      shift 2
      ;;
    --no-fresh)
      FRESH_RUN="0"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      warn "Unknown argument: $1"
      usage
      exit 2
      ;;
  esac
done

if [ -z "${TEMPLATE_CONFIG}" ]; then
  TEMPLATE_CONFIG="${PHASE2_DIR}/smoke_config.ini"
fi
if [ -z "${TELLURIC_SCRIPT}" ]; then
  TELLURIC_SCRIPT="${PHASE2_DIR}/telluric_spectra.py"
fi
if [ -z "${TABLES_PATH}" ]; then
  TABLES_PATH="${PHASE2_DIR}/spectra_results"
fi
if [ -z "${OUTPUT_BASE}" ]; then
  OUTPUT_BASE="${RUN_DIR}/spectra"
fi
if [ -z "${WORK_DIR}" ]; then
  WORK_DIR="${PHASE1_DIR}"
fi

RAW_SCRIPT="${PHASE1_DIR}/raw_table.py"
CLEAN_SCRIPT="${PHASE1_DIR}/cleantable.py"
SMOKE_SCRIPT="${PHASE2_DIR}/smoke_object_group_telluric_v2.py"

RAW_PATH="${WORK_DIR}/metadata_raw.parquet"
CLEAN_PATH="${WORK_DIR}/metadata.parquet"
FINAL_PATH="${WORK_DIR}/metadata_final.parquet"
REMAINING_PATH="${WORK_DIR}/metadata_remaining.parquet"

for required_file in "$RAW_SCRIPT" "$CLEAN_SCRIPT" "$SMOKE_SCRIPT" "$TEMPLATE_CONFIG" "$TELLURIC_SCRIPT"; do
  if [ ! -f "$required_file" ]; then
    warn "Required file not found: $required_file"
    exit 3
  fi
done

for required_dir in "$SCIENCE_PATH" "$CALIB_PATH" "$PHASE1_DIR" "$PHASE2_DIR"; do
  if [ ! -d "$required_dir" ]; then
    warn "Required directory not found: $required_dir"
    exit 4
  fi
done

mkdir -p "$WORK_DIR" "$OUTPUT_BASE" "$TABLES_PATH"
mkdir -p "$WORK_DIR/.matplotlib"
export MPLCONFIGDIR="$WORK_DIR/.matplotlib"

if [ "$FRESH_RUN" = "1" ]; then
  info "Fresh run enabled: removing previous phase1 parquet outputs"
  rm -f "$RAW_PATH" "$CLEAN_PATH" "$FINAL_PATH" "$REMAINING_PATH"
fi

if [ "$FRESH_RUN" = "1" ] || [ ! -f "$FINAL_PATH" ]; then
  run_step "Step 1/3 - raw_table.py" \
    "$PYTHON_BIN" "$RAW_SCRIPT" \
    --root_dir "$SCIENCE_PATH" \
    --output_dir "$WORK_DIR" \
    --output_name "metadata_raw.parquet" \
    --products "$PRODUCTS"
  status=$?
  if [ "$status" -ne 0 ]; then
    warn "Pipeline aborted at Step 1/3."
    exit "$status"
  fi

  run_step "Step 2/3 - cleantable.py" \
    "$PYTHON_BIN" "$CLEAN_SCRIPT" \
    --base_dir "$WORK_DIR" \
    --raw_path "$RAW_PATH" \
    --clean_path "$CLEAN_PATH" \
    --final_path "$FINAL_PATH" \
    --remaining_path "$REMAINING_PATH" \
    --radius_near "$RADIUS_NEAR"
  status=$?
  if [ "$status" -ne 0 ]; then
    warn "Pipeline aborted at Step 2/3."
    exit "$status"
  fi
else
  info "Reusing existing Stage 1 result: $FINAL_PATH"
fi

run_step "Step 3/3 - smoke_object_group_telluric_v2.py" \
  "$PYTHON_BIN" "$SMOKE_SCRIPT" \
  --parquet "$FINAL_PATH" \
  --template_config "$TEMPLATE_CONFIG" \
  --script "$TELLURIC_SCRIPT" \
  --output "$OUTPUT_BASE" \
  --calib "$CALIB_PATH" \
  --science_root "$SCIENCE_PATH" \
  --tables_path "$TABLES_PATH" \
  --star "$STAR" \
  --plot_orders "$PLOT_ORDERS" \
  --python "$PYTHON_BIN"
status=$?
if [ "$status" -ne 0 ]; then
  warn "Pipeline aborted at Step 3/3."
  exit "$status"
fi

info "Pipeline finished successfully."
exit 0
