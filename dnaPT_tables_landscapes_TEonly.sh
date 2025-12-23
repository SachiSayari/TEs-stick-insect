#!/bin/bash
#
# Changelog
# V.0 | 02.08.22 - first version (dnaPT_landscapes.sh)
# V.1 | 23.12.25 - TE-only scope (DNA, LINE, SINE, LTR, RC) + TE-normalized landscape
# V.2 | 23.12.25 - add percentage table output
#
set -euo pipefail

function usage() {
  cat << 'HEREDOC'
   **************************************
   >>>   dnaPT_landscapes_TEonly.sh   <<<
   **************************************
   TE landscape: histogram of blastn divergence between raw reads (TE copies) and
   their consensus sequences assembled in Trinity.fasta.
   Scope (hard filter): DNA, LINE, SINE, LTR, RC
   Dependencies:
   - R + packages: ggplot2, tidyr, dplyr
   Usage: ./dnaPT_landscapes_TEonly.sh -I <dataset_directory> [options]
HEREDOC
}

if [[ $# -eq 0 ]]; then
  echo "Error: No mandatory argument given."
  usage
  exit 1
fi

PARAMS=""
while (( "$#" )); do
  case "$1" in
    -I|--input-dir) DSA="$2"; shift 2 ;;
    -p|--prefix) PREF="$2"; shift 2 ;;
    -y|--ylim) YMAX="$2"; shift 2 ;;
    -o|--output) OUTF="$2"; shift 2 ;;
    -S|--superfamily) SF=TRUE; shift ;;
    -h|--help) usage; exit 0 ;;
    *) shift ;;
  esac
done

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
OUTF="${OUTF:-$DSA}"
PREF="${PREF:-dnaPipeTE}"
YMAX="${YMAX:-FALSE}"
SF="${SF:-FALSE}"

mkdir -p "$OUTF"

export DSA OUTF PREF YMAX SF DIR

Rscript - <<'RSCRIPT'
################################################################################
# packages
################################################################################
packages <- c("ggplot2", "tidyr", "dplyr")
invisible(lapply(packages, library, character.only = TRUE))

################################################################################
# variables
################################################################################
DSA <- Sys.getenv("DSA")
OUTF <- Sys.getenv("OUTF")
PREF <- Sys.getenv("PREF")
YMAX <- Sys.getenv("YMAX")
SF_CHOICE <- as.logical(Sys.getenv("SF"))
DIR <- Sys.getenv("DIR")

te_classes <- c("DNA", "LINE", "SINE", "LTR", "RC")

################################################################################
# load data
################################################################################
blast3 <- file.path(DSA, "Annotation", "sorted_blast3")
onehit <- file.path(DSA, "Annotation", "one_RM_hit_per_Trinity_contigs")

cmd <- paste0(
  "join -a1 -12 -21 -o 1.3,2.4,2.5 ",
  shQuote(blast3), " ", shQuote(onehit),
  " | awk '{if (NF == 1) {print $1\"\\tUnknown_repeat\\tUnknown\"} else {print $0}}'",
  " | sed 's/ /\\t/g'"
)

land <- read.table(
  sep = "\t",
  text = system(cmd, intern = TRUE),
  fill = TRUE,
  stringsAsFactors = FALSE
)

land <- tidyr::separate(land, V3, into = c("Sub_class", "SF"), sep = "/", fill = "right")
names(land) <- c("div", "TE_family", "TE_subclass", "TE_SF")

land$div <- as.numeric(land$div)
land$TE_subclass[land$TE_subclass == "MITE"] <- "DNA"

land <- land[land$TE_subclass %in% te_classes, , drop = FALSE]
reads_te <- nrow(land)

################################################################################
# percentage table
################################################################################
land$div_bin <- floor(100 - land$div)

if (SF_CHOICE) {
  perc_table <- land %>%
    group_by(div_bin, TE_SF) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percent_TE = count / reads_te * 100,
           TE_class = TE_SF)

  tab_out <- paste0(PREF, "_TEonly_landscapes_superfamily_percentages.tsv")

} else {
  perc_table <- land %>%
    group_by(div_bin, TE_subclass) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percent_TE = count / reads_te * 100,
           TE_class = TE_subclass)

  tab_out <- paste0(PREF, "_TEonly_landscapes_subclass_percentages.tsv")
}

write.table(
  perc_table,
  file = file.path(OUTF, tab_out),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

################################################################################
# plot
################################################################################
if (SF_CHOICE) {
  lscapes <- ggplot(land, aes(100 - div, fill = TE_SF)) +
    geom_histogram(aes(y = ..count.. / reads_te * 100), binwidth = 1)
  out <- paste0(PREF, "_TEonly_landscapes_superfamily.pdf")
} else {
  lscapes <- ggplot(land, aes(100 - div, fill = TE_subclass)) +
    geom_histogram(aes(y = ..count.. / reads_te * 100), binwidth = 1)
  out <- paste0(PREF, "_TEonly_landscapes_subclass.pdf")
}

lscapes <- lscapes +
  ylab("TE %") +
  xlab("blastn divergence (read vs dnaPipeTE contig)") +
  theme_classic()

ggsave(
  file = out,
  plot = lscapes,
  device = "pdf",
  path = OUTF,
  width = 2600,
  height = 2000,
  units = "px"
)

RSCRIPT

echo "All done. Plot + percentage table written to $OUTF"
