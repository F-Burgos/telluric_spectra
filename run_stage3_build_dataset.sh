#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${PROJECT_DIR}/phase3/build_telluric_dataset.py"

if [ ! -f "$SCRIPT" ]; then
  echo "Stage 3 script not found: $SCRIPT" >&2
  exit 3
fi

exec python "$SCRIPT" "$@"
