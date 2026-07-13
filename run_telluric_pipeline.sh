#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${PROJECT_DIR}/p_apoyo/p_apoyo/run_pipeline.sh"

if [ ! -x "$PIPELINE" ]; then
  echo "Canonical pipeline is not executable: $PIPELINE" >&2
  exit 3
fi

exec "$PIPELINE" "$@"
