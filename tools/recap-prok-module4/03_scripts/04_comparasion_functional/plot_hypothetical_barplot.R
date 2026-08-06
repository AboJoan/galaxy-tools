#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_hypothetical_barplot.R <in.tsv> <out.pdf> <out.tiff>")
}

in_tsv  <- args[1]
out_pdf <- args[2]
out_tif <- args[3]

tab <- read.delim(in_tsv, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Keep only "Hypothetical proteins" rows
tab <- tab[tab$Metric == "Hypothetical proteins", , drop = FALSE]
if (nrow(tab) == 0) stop("No rows with Metric == 'Hypothetical proteins' found.")

# Convert "123 (45.67%)" -> numeric percent
extract_pct <- function(x) {
  m <- regmatches(x, regexec("\\(([-0-9.]+)%\\)", x))
  sapply(m, function(v) if (length(v) >= 2) as.numeric(v[2]) else NA_real_)
}

# Column names (3rd and 4th columns are the two compared datasets)
col_re <- colnames(tab)[3]  # e.g., "Re-annotated"
col_rf <- colnames(tab)[4]  # e.g., "NCBI Reference genome"

tab$Reann_pct <- extract_pct(tab[[col_re]])
tab$Ref_pct   <- extract_pct(tab[[col_rf]])

if (all(is.na(tab$Reann_pct)) || all(is.na(tab$Ref_pct))) {
  stop("Could not parse percentages from table columns 3/4. Expected format like '123 (45.67%)'.")
}

# ---- Option C: show Species + ID together (nice labels) ----
# We prefer "Species (ID)" if possible, otherwise keep whatever is already there.
# If your compare script already writes "Species (Re-annotated_01)" this does nothing harmful.
make_species_id_label <- function(x) {
  # If already contains parentheses, assume it already has the ID
  ifelse(grepl("\\(.+\\)$", x), x, x)
}
tab$XLabel <- make_species_id_label(tab$Species)

# Order by Re-annotation percent (nice look)
ord <- order(tab$Reann_pct, decreasing = TRUE, na.last = TRUE)
tab <- tab[ord, , drop = FALSE]

# Build matrix for base R barplot:
# Put Reference first (left bar), Re-annotated second (right bar)
m <- rbind(tab$Ref_pct, tab$Reann_pct)
colnames(m) <- tab$XLabel

# Colors (NCBI = red, Re-annotated = teal)
col_ref  <- "#4C72B0"
col_rean <- "#DD8452"
bar_cols <- c(col_ref, col_rean)

plot_one <- function(device_open) {
  device_open()
  
  par(mar = c(12, 5, 4, 1), xpd = NA)
  
  bp <- barplot(
    m,
    beside = TRUE,
    las = 2,
    ylim = c(0, 100),
    ylab = "Hypothetical proteins (%)",
    main = "Hypothetical proteins: Re-annotated vs Reference",
    col = bar_cols,
    border = NA
  )
  
  # Legend (match bar order)
  legend(
    "top",
    inset = c(0, -0.10),
    horiz = TRUE,
    bty = "n",
    fill = bar_cols,
    legend = c(col_rf, col_re)
  )
  
  # Add percent labels above bars
  vals <- as.vector(m)
  text(
    x = as.vector(bp),
    y = vals,
    labels = ifelse(is.na(vals), "", sprintf("%.1f%%", vals)),
    pos = 3,
    cex = 0.85,
    offset = 0.35
  )
  
  dev.off()
}

# PDF
plot_one(function() pdf(out_pdf, width = 14, height = 6))

# TIFF
plot_one(function() tiff(out_tif, width = 14, height = 6, units = "in", res = 300, compression = "lzw"))