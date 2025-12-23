#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

# CONTAINERS (safe with set -u)
: "${CONTAINERS:?ERROR: \$CONTAINERS environment variable is not set.}"

IMG="$CONTAINERS/dnapipete.img"

# fastq files
CLEANED_FAST_DIR="/DATABIG/sara.sebestova/SRAs/B_rossius/cleaned_fastq"

# library
RP_LIB_DIR="/DATABIG/sara.sebestova/SRAs/Bacillus_lib"
RP_LIB="cleaned_classified_library.fa"

# parameters
GENOME_SIZE=2000000000
COV=0.5
SAMPLE_NUMBER=2
RM_T=0.25
CPU=4

[[ -f "$IMG" ]] || { echo "ERROR: Singularity image not found: $IMG" >&2; exit 1; }
[[ -d "$CLEANED_FAST_DIR" ]] || { echo "ERROR: FASTQ directory not found: $CLEANED_FAST_DIR" >&2; exit 1; }
[[ -d "$RP_LIB_DIR" ]] || { echo "ERROR: Repeat Library directory not found: $RP_LIB_DIR" >&2; exit 1; }
[[ -f "$RP_LIB_DIR/$RP_LIB" ]] || { echo "ERROR: Repeat library file not found: $RP_LIB_DIR/$RP_LIB" >&2; exit 1; }

command -v singularity >/dev/null 2>&1 || { echo "ERROR: singularity command not found." >&2; exit 1; }

cd "$CLEANED_FAST_DIR"

# gzipped FASTQs
fastqs=( *clean.R1.fastq *clean.R1.fastq.gz )

if (( ${#fastqs[@]} == 0 )); then
  echo "ERROR: No FASTQ files found matching '*clean.R1.fastq' or '*clean.R1.fastq.gz' in $CLEANED_FAST_DIR" >&2
  exit 1
fi

echo "Found ${#fastqs[@]} FASTQ files."

run_job() {
  local fq="$1"

  # sample name 
  local base="${fq%.gz}"
  local sample="${base%.clean.R1.fastq}"
  local outdir="newRPlib_05_${sample}"

  # Skip if output exists and has content
  if [[ -d "$outdir" ]] && find "$outdir" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    echo ">>> Skipping $sample (output directory exists and is not empty): $outdir"
    return 0
  fi

  echo ">>> Starting sample: $sample"

  if singularity exec \
      --bind "${CLEANED_FAST_DIR}:/data/fastq,${RP_LIB_DIR}:/data/RP_lib,/DATABIG/sara.sebestova/tmp:/DATABIG/sara.sebestova/tmp" \
      "$IMG" \
      bash -lc "
        set -euo pipefail
        cd /opt/dnaPipeTE/

	export TMPDIR=/DATABIG/sara.sebestova/tmp
	export TMP=/DATABIG/sara.sebestova/tmp
	export TEMP=/DATABIG/sara.sebestova/tmp
        export JAVA_HOME=/opt/dnaPipeTE/bin/OpenJDK-1.8.0.141-x86_64-bin
        export PATH=\$JAVA_HOME/bin:\$PATH
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
  then
    echo ">>> Successfully finished: $sample"
  else
    local rc=$?
    echo ">>> WARNING: $sample FAILED (exit code $rc)" >&2
    return "$rc"
  fi
}

for fq in "${fastqs[@]}"; do
  run_job "$fq" || true
  echo "--- Moving to the next sample ---"
done

echo ">>> All jobs finished."
 
