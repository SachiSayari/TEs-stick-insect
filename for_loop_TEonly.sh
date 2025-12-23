#!/bin/bash

BASE_DIR="/DATABIG/sara.sebestova/SRAs/B_grandii/cleaned_fastq/male_reproductive_tissue"
SCRIPT="$BASE_DIR/dnaPT_charts_TEonly.sh"

for folder in "$BASE_DIR"/*/; do
    echo "Processing folder: $folder"
    
    # cd into the folder and run the script
    (cd "$folder" && bash "$SCRIPT" -I "$folder")
done

