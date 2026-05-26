#!/bin/bash
# process_directory.sh
#
# TEMPLATE for batch processing. Open this file in a text editor and set
# BOLD_DIR and OUTPUT_DIR to your real paths before running. Then run from
# a terminal with:
#     bash examples/process_directory.sh
#
# This script loops over every '*_desc-preproc_bold.nii.gz' in BOLD_DIR
# and calls reparcellate.py on each, writing the output TSVs to OUTPUT_DIR.

set -euo pipefail

# === Edit these two paths ===
BOLD_DIR="/path/to/bolds"
OUTPUT_DIR="/path/to/outputs"
# ============================

# Safety check: refuse to run with placeholder paths.
if [ "$BOLD_DIR" = "/path/to/bolds" ] || [ "$OUTPUT_DIR" = "/path/to/outputs" ]; then
    echo "ERROR: This is a template script. Please open process_directory.sh"
    echo "       in a text editor and set BOLD_DIR and OUTPUT_DIR to your"
    echo "       real paths before running." >&2
    exit 1
fi

REPARC_SCRIPT="$(dirname "$0")/../reparcellate.py"

if [ ! -f "$REPARC_SCRIPT" ]; then
    echo "ERROR: Cannot find reparcellate.py at $REPARC_SCRIPT" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
bolds=("$BOLD_DIR"/*_desc-preproc_bold.nii.gz)
n=${#bolds[@]}
if [ "$n" -eq 0 ]; then
    echo "No *_desc-preproc_bold.nii.gz files found in $BOLD_DIR" >&2
    exit 1
fi
echo "Found $n BOLD file(s) to process."

i=0
for bold in "${bolds[@]}"; do
    i=$((i + 1))
    echo ""
    echo "===== [$i/$n] $(basename "$bold") ====="
    python "$REPARC_SCRIPT" --bold "$bold" --out "$OUTPUT_DIR"
done

echo ""
echo "All $n subject(s) processed. Output: $OUTPUT_DIR"
