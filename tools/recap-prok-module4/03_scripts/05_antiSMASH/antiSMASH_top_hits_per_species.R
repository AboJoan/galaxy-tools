#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
})

bgc_tsv <- Sys.getenv("BGC_TSV")
outdir  <- Sys.getenv("OUTDIR")

if (bgc_tsv == "" || outdir == "") {
  stop("Missing ENV vars: BGC_TSV / OUTDIR")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read_tsv(bgc_tsv, show_col_types = FALSE)

stopifnot(all(c("Genomes", "Canonical_Class", "BGC_count") %in% colnames(df)))

class_order <- c(
  "Terpene",
  "RiPP",
  "PKS",
  "Phosphonate",
  "Deazapurine",
  "NRPS",
  "Auto-inducer",
  "Other"
)

df <- df %>%
  mutate(
    Genomes = as.character(Genomes),
    Canonical_Class = as.character(Canonical_Class),
    BGC_count = as.numeric(BGC_count)
  )

missing_classes <- setdiff(unique(df$Canonical_Class), class_order)
class_order_final <- c(class_order, sort(missing_classes))

genome_order <- df %>%
  group_by(Genomes) %>%
  summarise(total = sum(BGC_count, na.rm = TRUE), .groups = "drop") %>%
  arrange(total, Genomes) %>%
  pull(Genomes)

df <- df %>%
  mutate(
    Genomes = factor(Genomes, levels = genome_order),
    Canonical_Class = factor(Canonical_Class, levels = class_order_final)
  )

p <- ggplot(df, aes(x = Genomes, y = BGC_count, fill = Canonical_Class)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of BGCs",
    fill = "BGC class",
    title = "BGC classes per re-annotated genome"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

png_out <- file.path(outdir, "BGC_classes_per_genome.png")
pdf_out <- file.path(outdir, "BGC_classes_per_genome.pdf")

ggsave(png_out, p, width = 10, height = 6, dpi = 300)
ggsave(pdf_out, p, width = 10, height = 6)

cat("Saved:", png_out, "\n")
cat("Saved:", pdf_out, "\n")
