#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <sample_result_dir> [plot options...]" >&2
    exit 1
fi

SAMPLE_DIR="$1"
shift

exec "$PYTHON" "$SCRIPT_DIR/PlotHCMExonDelMultiTrack.py" "$SAMPLE_DIR" "$@"
