#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  setup_data_links.sh --science PATH --calib PATH [--force]

Creates the local data docking folder expected by the pipeline:

  Data/science -> PATH
  Data/calib   -> PATH

The script creates symbolic links only. It does not copy or move the underlying
HARPS data.

Options:
  --science PATH  Real HARPS science root
  --calib PATH    Real HARPS calibration root
  --force         Replace existing Data/science or Data/calib entries
  -h, --help      Show this help message
EOF
}

warn() {
  echo "[WARNING] $*" >&2
}

info() {
  echo "[INFO] $*"
}

SCIENCE_PATH=""
CALIB_PATH=""
FORCE="0"

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
    --force)
      FORCE="1"
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

if [ -z "$SCIENCE_PATH" ] || [ -z "$CALIB_PATH" ]; then
  warn "Both --science and --calib are required."
  usage
  exit 2
fi

if [ ! -d "$SCIENCE_PATH" ]; then
  warn "Science directory not found: $SCIENCE_PATH"
  exit 4
fi

if [ ! -d "$CALIB_PATH" ]; then
  warn "Calibration directory not found: $CALIB_PATH"
  exit 4
fi

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${PROJECT_DIR}/Data"
mkdir -p "$DATA_DIR"

link_path() {
  local target="$1"
  local link="$2"
  local label="$3"

  target="$(realpath "$target")"

  if [ -e "$link" ] || [ -L "$link" ]; then
    if [ "$FORCE" != "1" ]; then
      warn "$link already exists. Use --force to replace it."
      exit 5
    fi
    rm -rf "$link"
  fi

  ln -s "$target" "$link"
  info "${label}: ${link} -> ${target}"
}

link_path "$SCIENCE_PATH" "${DATA_DIR}/science" "Science"
link_path "$CALIB_PATH" "${DATA_DIR}/calib" "Calibration"

info "Data links are ready. The pipeline defaults now use Data/science and Data/calib."
