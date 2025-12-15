#!/usr/bin/env bash

set -euo pipefail

### USER CONFIGURATION ###

# Path to Singularity image

if [[ -z "${CONTAINERS}" ]]; then
    echo "ERROR: \$CONTAINERS environment variable is not set." >&2
    exit 1
fi
IMG="$CONTAINERS/dnapipete.img"

# Species directory containing fastq files
CLEANED_FAST_DIR="/DATABIG/sara.sebestova/SRAs/B_grandii/cleaned_fastq"

# Repeat library directory
RP_LIB_DIR="/DATABIG/sara.sebestova/SRAs/Bacillus_lib"
RP_LIB="cleaned_classified_library.fa"

# dnaPipeTE parameters
GENOME_SIZE=2000000000
COV=0.1
SAMPLE_NUMBER=2
RM_T=0.25
CPU=6

### PRE-RUN CHECKS ###

# Check if directories exist
if [[ ! -d "$CLEANED_FAST_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $CLEANED_FAST_DIR" >&2
    exit 1
fi

if [[ ! -d "$RP_LIB_DIR" ]]; then
    echo "ERROR: Repeat Library directory not found: $RP_LIB_DIR" >&2
    exit 1
fi

if ! command -v singularity &> /dev/null; then
    echo "ERROR: singularity command not found." >&2
    exit 1
fi

### RUN PIPELINES ###

# Change to the directory where FASTQ files are located
cd "$CLEANED_FAST_DIR"

# Get all R1 FASTQ files
fastqs=( *clean.R1.fastq )

# Check if any FASTQ files were found
if [[ ${#fastqs[@]} -eq 0 ]]; then
    echo "ERROR: No FASTQ files found matching '*clean.R1.fastq' in $CLEANED_FAST_DIR" >&2
    exit 1
fi

echo "Found ${#fastqs[@]} FASTQ files. Time to say FINGERS CROSSED."

# --- Function to run dnaPipeTE for a single sample ---
run_job() {
    fq="$1"
    # Derive sample name from the FASTQ filename
    sample=$(basename "$fq" .clean.R1.fastq)
    # Since we are in $CLEANED_FAST_DIR, this path is relative
    outdir="newRPlib_${sample}"

    # SKIP CHECK: Check if the output directory exists and is not empty.
    if [ -d "$outdir" ] && [ -n "$(ls -A "$outdir" 2>/dev/null)" ]; then
        echo ">>> Skipping $sample (output directory exists and is not empty)"
        return
    fi

    echo ">>> Starting a sample: $sample"
   
    # Run dnaPipeTE inside the Singularity container
    singularity exec \
        --bind "$CLEANED_FAST_DIR":/data/fastq,"$RP_LIB_DIR":/data/RP_lib \
        "$IMG" \
        bash -lc "
    cd /opt/dnaPipeTE/

    # Explicitly use container-bundled Java 8 (required by Trinity)
    # Without specifying the script will fail on which Java version should be used
    export JAVA_HOME=/opt/dnaPipeTE/bin/OpenJDK-1.8.0.141-x86_64-bin
    export PATH=\$JAVA_HOME/bin:\$PATH

    # Set locale to avoid Python/Perl Unicode warnings 
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

            python3 dnaPipeTE.py \
                -input /data/fastq/$fq \
                -output /data/fastq/$outdir \
                -genome_size $GENOME_SIZE \
                -genome_coverage $COV \
                -sample_number $SAMPLE_NUMBER \
                -RM_lib /data/RP_lib/$RP_LIB \
                -RM_t $RM_T \
                -cpu $CPU
        "

    # Check the exit status of the singularity command
    if [ $? -eq 0 ]; then
        echo ">>> Successfully Finished: $sample"
    else
        echo ">>> WARNING: $sample FAILED with exit code $?" >&2
    fi
}
# --- End of run_job function ---

### SEQUENTIAL EXECUTION ###

# Loop through each FASTQ file and run the job function
for fq in "${fastqs[@]}"; do
    run_job "$fq"
    echo "--- Moving to the next sample ---"
done

echo ">>> All jobs finished! FINGERS UNCROSSED."
