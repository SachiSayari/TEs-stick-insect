#!/bin/bash

# -----------------------------------------------------------
# dnaPT_percentages.sh
# Extract percentages from dnaPipeTE output
# -----------------------------------------------------------

function usage() {
    echo "Usage: $0 -I <dnaPipeTE_output_folder>"
    exit 1
}

# --------------------------
# Parse arguments
# --------------------------
while (( "$#" )); do
  case "$1" in
    -I|--input-dir)
      DSA=$2
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown argument: $1"
      usage
      ;;
  esac
done

if [[ -z "$DSA" ]]; then
    echo "Error: missing input directory"
    usage
fi

echo ">>> Processing dnaPipeTE folder:"
echo "$DSA"
echo ""

# --------------------------
# Run R to compute percentages
# --------------------------
Rscript - <<EOF
library(tidyverse)

# -------- Read dnaPipeTE table --------
rca <- read.table("$DSA/reads_per_component_and_annotation",
                  fill = TRUE, na.strings=c("","NA"))

# Split annotation field
rca <- separate(rca, V6, c("subclass","superfamily"), sep="/", fill="right")
rca\$superfamily[is.na(rca\$superfamily)] <- "Unknown"
rca\$subclass[is.na(rca\$subclass)] <- "Unknown"
rca\$superfamily <- paste(rca\$subclass, rca\$superfamily, sep="/")
rca\$superfamily[rca\$superfamily == "Unknown/Unknown"] <- "Unknown"

# Rename columns
names(rca) <- c("reads","bp","contig","q_length","target_name","t_subc","t_supf","hit_cov")

# -------- Get total mapped bp from Counts.txt --------
reads.b <- as.numeric(system("tail -n 2 $DSA/Counts.txt | head -n 1 | cut -f 2", intern=TRUE))

# -------- Compute % genome for each contig --------
rca <- rca[order(rca\$bp, decreasing=TRUE),]
rca\$perc_g <- rca\$bp / reads.b * 100

write.table(
    rca[,c("contig","t_subc","bp","perc_g")],
    file="$DSA/contig_percentages.tsv",
    sep="\t", row.names=FALSE, quote=FALSE
)

# -------- Compute % per repeat subclass --------
rca_pie <- rca %>%
    group_by(t_subc) %>%
    summarise(perc_g = sum(perc_g)) %>%
    mutate(perc_r = perc_g / sum(perc_g) * 100)

write.table(
    rca_pie,
    file="$DSA/repeat_class_percentages.tsv",
    sep="\t", row.names=FALSE, quote=FALSE
)

cat("Generated tables:\n")
cat("   $DSA/contig_percentages.tsv\n")
cat("   $DSA/repeat_class_percentages.tsv\n")

EOF

echo ">>> Done."
