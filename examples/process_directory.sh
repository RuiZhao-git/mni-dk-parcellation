#!/bin/bash
# process_directory.sh
# Example: process every BOLD NIfTI in a directory through the parcellation
# pipeline. Adjust BOLD_DIR and OUTPUT_DIR to your paths.

set -euo pipefail

# === Edit these two paths ===
BOLD_DIR="/path/to/bolds"
OUTPUT_DIR="/path/to/outputs"
# ============================

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
