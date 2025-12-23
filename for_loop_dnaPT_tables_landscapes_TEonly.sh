#!/bin/bash

BASE_DIR="/DATABIG/sara.sebestova/SRAs/B_atticus/cleaned_fastq/reproductive_tissue"
SCRIPT="$BASE_DIR/dnaPT_tables_landscapes_TEonly.sh"

for FG_DIR in "$BASE_DIR"/*_FG/; do
    echo "Processing: $FG_DIR"
    bash "$SCRIPT" -I "$FG_DIR"
done
