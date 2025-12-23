#!/bin/bash
#
# dnaPipeTE charts - TE-only (DNA, LINE, SINE, LTR, RC)
# Base: Clément Goubert - goubert.clement@gmail.com
# Modified: TE-only scope + robust typing + optional Unknown exclusion for genome denominator

set -euo pipefail

function usage() {
   cat << 'HEREDOC'

   **************************************
   >>>  dnaPT_charts_TEonly.sh  <<<
   **************************************

   Usage:
     ./dnaPT_charts_TEonly.sh -I <dnaPipeTE_output_dir> [options]

   Mandatory:
     -I, --input-dir              dnaPipeTE output directory

   Options:
     -p, --prefix                 prefix for output files (default: dnaPipeTE)
     -t, --percent_threshold      contigs under threshold (%) are grouped (default: 0.001)
                                  NOTE: threshold is applied on TE-total denominator
     -o, --output                 output folder (default: input directory)
     -y, --ylim                   max value for y axis (0-100); applies to TE% axis
     -U, --exclude-unknown        exclude contigs with t_subc == "Unknown" (source table)
                                  and exclude Unknown bp from genome denominator (genome pie)
     -h, --help                   print this message and exit

   Scope:
     TE-only = {DNA, LINE, SINE, LTR, RC}
     Everything else is excluded before plotting.

   Outputs:
     If -U is OFF:
       <prefix>_TEonly_charts.pdf
       <prefix>_TEonly_barplot_table.tsv
       <prefix>_TEonly_pie_relative_table.tsv
       <prefix>_TEonly_pie_genome_table.tsv

     If -U is ON:
       <prefix>_TEonly_charts_noUnknown.pdf
       <prefix>_TEonly_barplot_table_noUnknown.tsv
       <prefix>_TEonly_pie_relative_table_noUnknown.tsv
       <prefix>_TEonly_pie_genome_table_noUnknown.tsv

HEREDOC
}

if [[ $# -eq 0 ]]; then
  echo "Error: No arguments given."
  usage
  exit 1
fi

PARAMS=""
while (( "$#" )); do
  case "$1" in
    -I|--input-dir)
      if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
        DSA="$2"
        shift 2
      else
        echo "Error: missing dnaPipeTE directory after $1" >&2
        usage
        exit 1
      fi
      ;;
    -p|--prefix)
      if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
        PREF="$2"
        shift 2
      else
        shift
      fi
      ;;
    -y|--ylim)
      if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
        YMAX="$2"
        shift 2
      else
        shift
      fi
      ;;
    -o|--output)
      if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
        OUTF="$2"
        shift 2
      else
        echo "Error: missing output folder after $1" >&2
        usage
        exit 1
      fi
      ;;
    -t|--percent_threshold)
      if [[ -n "${2:-}" && "${2:0:1}" != "-" ]]; then
        TPERC="$2"
        shift 2
      else
        shift
      fi
      ;;
    -U|--exclude-unknown)
      EXCLUDE_UNKNOWN="TRUE"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*|--*=)
      echo "Error: Unsupported argument $1" >&2
      usage
      exit 1
      ;;
    *)
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done
eval set -- "$PARAMS"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

OUTF="${OUTF:-$DSA}"
PREF="${PREF:-dnaPipeTE}"
YMAX="${YMAX:-FALSE}"
TPERC="${TPERC:-0.001}"
EXCLUDE_UNKNOWN="${EXCLUDE_UNKNOWN:-FALSE}"

echo "input dataset:      $DSA"
echo "output folder:      $OUTF"
echo "output prefix:      $PREF"
echo "custom ylim:        $YMAX"
echo "exclude Unknown:    $EXCLUDE_UNKNOWN"
echo "scope:              TE-only (DNA, LINE, SINE, LTR, RC)"

mkdir -p "$OUTF"

export DSA OUTF PREF YMAX TPERC EXCLUDE_UNKNOWN DIR

Rscript - <<'RSCRIPT'
################################################################################
# packages
################################################################################
packages <- c("ggplot2", "tidyverse", "cowplot")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  message("ERROR: missing R packages: ", paste(packages[!installed_packages], collapse = ", "))
  quit(status = 1)
}
invisible(lapply(packages, library, character.only = TRUE))

################################################################################
# parameters from env
################################################################################
DSA <- Sys.getenv("DSA")
OUTF <- Sys.getenv("OUTF")
PREF <- Sys.getenv("PREF")
YMAX <- Sys.getenv("YMAX")
TPERC <- as.numeric(Sys.getenv("TPERC"))
EXCLUDE_UNKNOWN <- as.logical(Sys.getenv("EXCLUDE_UNKNOWN"))
DIR <- Sys.getenv("DIR")
ymax <- YMAX

te_classes <- c("DNA", "LINE", "SINE", "LTR", "RC")

################################################################################
# read input table
################################################################################
anno_file <- file.path(DSA, "reads_per_component_and_annotation")
if (!file.exists(anno_file)) {
  message("ERROR: cannot find ", anno_file)
  quit(status = 1)
}

rca <- read.table(
  anno_file,
  fill = TRUE,
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

rca <- tidyr::separate(rca, V6, into = c("subclass", "superfamily"), sep = "/", fill = "right")
rca$superfamily[is.na(rca$superfamily)] <- "Unknown"
rca$subclass[is.na(rca$subclass)] <- "Unknown"

rca$superfamily <- paste(rca$subclass, rca$superfamily, sep = "/")
rca$superfamily[rca$superfamily == "Unknown/Unknown"] <- "Unknown"

if (ncol(rca) < 8) {
  rca$V7 <- NA
}
names(rca) <- c("reads", "bp", "contig", "q_length", "target_name", "t_subc", "t_supf", "hit_cov")

rca$reads <- suppressWarnings(as.numeric(rca$reads))
rca$bp    <- suppressWarnings(as.numeric(rca$bp))

################################################################################
# genome proxy + unknown bp (for genome denominator if -U)
################################################################################
reads.b <- as.numeric(system(paste0("tail -n 2 ", shQuote(file.path(DSA, "Counts.txt")), " | head -n 1 | cut -f 2"), intern = TRUE))
unknown_bp_total <- sum(rca[rca$t_subc == "Unknown", "bp"], na.rm = TRUE)

denom_genome <- reads.b
if (EXCLUDE_UNKNOWN) {
  denom_genome <- reads.b - unknown_bp_total
  if (is.na(denom_genome) || denom_genome <= 0) {
    message("WARNING: genome denominator without Unknown is <= 0; using total genome denominator.")
    denom_genome <- reads.b
  }
}

################################################################################
# scope selection: (optional) remove Unknown class + TE-only filter
################################################################################
rca_use <- rca
if (EXCLUDE_UNKNOWN) {
  rca_use <- rca_use[rca_use$t_subc != "Unknown", , drop = FALSE]
}

# TE-only
rca_use <- rca_use[rca_use$t_subc %in% te_classes, , drop = FALSE]

if (nrow(rca_use) == 0) {
  message("ERROR: After TE-only filtering (DNA/LINE/SINE/LTR/RC), no rows remain.")
  quit(status = 1)
}

rca_use <- rca_use[order(rca_use$bp, decreasing = TRUE), , drop = FALSE]

denom_te_bp <- sum(rca_use$bp, na.rm = TRUE)
if (is.na(denom_te_bp) || denom_te_bp <= 0) {
  message("ERROR: TE denominator <= 0 (sum of bp for TE-only rows).")
  quit(status = 1)
}

################################################################################
# thresholding on TE denominator
################################################################################
bp_cutoff <- TPERC * denom_te_bp / 100
leftover_label <- paste0("TE_under_", TPERC, "%")

rca_t <- rca_use[rca_use$bp >= bp_cutoff, , drop = FALSE]

leftover_df <- data.frame(
  reads       = sum(rca_use[rca_use$bp < bp_cutoff, "reads"], na.rm = TRUE),
  bp          = sum(rca_use[rca_use$bp < bp_cutoff, "bp"],    na.rm = TRUE),
  contig      = NA_character_,
  q_length    = NA_real_,
  target_name = NA_character_,
  t_subc      = leftover_label,
  t_supf      = leftover_label,
  hit_cov     = NA_real_,
  stringsAsFactors = FALSE
)

rca_t <- dplyr::bind_rows(rca_t, leftover_df)
rca_t$reads <- as.numeric(rca_t$reads)
rca_t$bp    <- as.numeric(rca_t$bp)

################################################################################
# TE-normalized plotting percentages
################################################################################
rca_t$perc_plot <- rca_t$bp / denom_te_bp * 100
rca_t$x2 <- cumsum(rca_t$perc_plot)
rca_t$x1 <- rca_t$x2 - rca_t$perc_plot

rca_t$t_subc <- factor(rca_t$t_subc, levels = unique(rca_t$t_subc))

suffix <- if (EXCLUDE_UNKNOWN) "_noUnknown" else ""

write.table(
  rca_t,
  file = file.path(OUTF, paste0(PREF, "_TEonly_barplot_table", suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

################################################################################
# colors
################################################################################
cols_file <- file.path(DIR, "colors.land")
if (!file.exists(cols_file)) {
  message("ERROR: cannot find colors file: ", cols_file)
  quit(status = 1)
}
cols <- read.table(cols_file, sep = "\t", stringsAsFactors = FALSE)
cols[nrow(cols) + 1, ] <- c(leftover_label, "grey10")

classes <- levels(rca_t$t_subc)
col.bars <- rep("grey50", length(classes))

print("Recognized classes (TE-only):")
for (i in seq_along(classes)) {
  print(paste0("^", classes[i], "$"))
  hit <- cols$V2[grep(pattern = paste0("^", classes[i], "$"), x = cols$V1)]
  if (length(hit) >= 1 && !is.na(hit[1])) {
    col.bars[i] <- hit[1]
  }
}

################################################################################
# barplot (TE %)
################################################################################
barplot <- ggplot(rca_t) +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = perc_plot, fill = t_subc)) +
  scale_fill_manual(values = col.bars, name = "TE class") +
  xlab("cumulative TE %") +
  ylab("TE % (composition)") +
  {if (ymax != FALSE) ylim(0, as.numeric(ymax))} +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

################################################################################
# pie charts
################################################################################
# Relative TE composition (sums to 100 by construction)
rca_pie <- rca_t %>%
  group_by(t_subc) %>%
  summarise(perc_g = sum(perc_plot, na.rm = TRUE), .groups = "drop") %>%
  mutate(subc = as.character(t_subc)) %>%
  select(perc_g, subc)

rca_pie$perc_r <- rca_pie$perc_g / sum(rca_pie$perc_g) * 100

write.table(
  rca_pie,
  file = file.path(OUTF, paste0(PREF, "_TEonly_pie_relative_table", suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

df2 <- rca_pie %>%
  mutate(
    csum = rev(cumsum(rev(perc_r))),
    pos  = perc_r / 2 + lead(csum, 1),
    pos  = if_else(is.na(pos), perc_r / 2, pos)
  )
df2$subc <- factor(df2$subc, levels = df2$subc, labels = paste0(round(df2$perc_r, 2), "% ", df2$subc))

pie_rel <- ggplot(df2, aes(x = "", y = perc_r, fill = subc)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  scale_fill_manual(values = col.bars, name = "TE class (composition)") +
  theme_void()

# Genome pie: TE fraction of genome (+ non-TE remainder)
rca_pie_g <- rca_t %>%
  group_by(t_subc) %>%
  summarise(bp = sum(bp, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    perc_g = bp / denom_genome * 100,
    subc   = as.character(t_subc)
  ) %>%
  select(perc_g, subc)

non_te <- 100 - sum(rca_pie_g$perc_g, na.rm = TRUE)
rca_pie_g <- rbind(rca_pie_g, data.frame(perc_g = non_te, subc = "non-TE"))
rca_pie_g$perc_g <- as.numeric(rca_pie_g$perc_g)

# Add a color for non-TE
col.bars_g <- c(col.bars, "grey90")

write.table(
  rca_pie_g,
  file = file.path(OUTF, paste0(PREF, "_TEonly_pie_genome_table", suffix, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

dfg <- rca_pie_g %>%
  mutate(
    csum = rev(cumsum(rev(perc_g))),
    pos  = perc_g / 2 + lead(csum, 1),
    pos  = if_else(is.na(pos), perc_g / 2, pos)
  )
dfg$subc <- factor(dfg$subc, levels = dfg$subc, labels = paste0(round(dfg$perc_g, 2), "% ", dfg$subc))

pie_g <- ggplot(dfg, aes(x = "", y = perc_g, fill = subc)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  scale_fill_manual(values = col.bars_g, name = if (EXCLUDE_UNKNOWN) "TE class (genome %, Unknown excluded)" else "TE class (genome %)") +
  theme_void()

################################################################################
# output PDF
################################################################################
top <- plot_grid(pie_rel, pie_g, ncol = 2)
charts <- plot_grid(top, barplot, ncol = 1)

outfile <- if (EXCLUDE_UNKNOWN) paste0(PREF, "_TEonly_charts_noUnknown.pdf") else paste0(PREF, "_TEonly_charts.pdf")

save_plot(
  file = outfile,
  charts,
  path = OUTF,
  base_width = 10,
  base_height = 8
)
RSCRIPT

echo "All done, results in $OUTF/"

