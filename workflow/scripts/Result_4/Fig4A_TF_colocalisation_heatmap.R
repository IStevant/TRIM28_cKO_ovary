source(".Rprofile")

###########################################
# Libraries
###########################################

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("rtracklayer")
  library("ggplot2")
  library("cowplot")
})

###########################################
# Logging
###########################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

###########################################
# Load data
###########################################

TF_list <- list(
  FOXL2 = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2 = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2  = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX = rtracklayer::import(snakemake@input[["RUNX"]])
)

TF_names <- names(TF_list)

###########################################
# Functions
###########################################

compute_overlap_percentage <- function(
  query,
  subject
) {
  overlap <- GenomicRanges::countOverlaps(
    query,
    subject,
    ignore.strand = TRUE
  ) > 0

  100 * sum(overlap) / length(query)
}

###########################################
# Pairwise overlap matrix
###########################################

log_message("Computing pairwise TF overlap matrix")

overlap_matrix <- matrix(
  NA_real_,
  nrow = length(TF_names),
  ncol = length(TF_names),
  dimnames = list(
    TF_names,
    TF_names
  )
)

for (TF1 in TF_names) {

  for (TF2 in TF_names) {

    overlap_matrix[TF1, TF2] <- compute_overlap_percentage(
      query = TF_list[[TF1]],
      subject = TF_list[[TF2]]
    )

  }

}

###########################################
# Export table
###########################################

write.table(
  round(overlap_matrix, 2),
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

###########################################
# Prepare plotting table
###########################################

plot_table <- expand.grid(
  TF_row = TF_names,
  TF_col = TF_names,
  stringsAsFactors = FALSE
)

plot_table$overlap_pct <- as.vector(
  overlap_matrix
)

overlap_counts <- matrix(
  NA_integer_,
  nrow = length(TF_names),
  ncol = length(TF_names),
  dimnames = list(TF_names, TF_names)
)

for (TF1 in TF_names) {
  for (TF2 in TF_names) {

    overlap_counts[TF1, TF2] <- sum(
      GenomicRanges::countOverlaps(
        TF_list[[TF1]],
        TF_list[[TF2]],
        ignore.strand = TRUE
      ) > 0
    )

  }
}

plot_table$overlap_n <- as.vector(
  overlap_counts
)

plot_table$label <- paste0(
  round(plot_table$overlap_pct),
  "%\n(",
  scales::comma(plot_table$overlap_n),
  ")"
)

plot_table$TF_row <- factor(
  plot_table$TF_row,
  levels = rev(TF_names)
)

plot_table$TF_col <- factor(
  plot_table$TF_col,
  levels = TF_names
)

###########################################
# Heatmap
###########################################

heatmap_plot <- ggplot(
  plot_table,
  aes(
    x = TF_col,
    y = TF_row,
    fill = overlap_pct
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = label),
    size = 4
  ) +
  scale_fill_gradient(
    low = "#FFF9EC",
    high = "#7D1128",
    limits = c(0, 100)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "% overlap",
    title = "Pairwise overlap between ovarian transcription factors"
  ) +
  theme_light() +
  theme(
    plot.title = element_text(
      size = 12,
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.text.y = element_text(size = 12)
  )

###########################################
# Save figure
###########################################

save_plot(
  snakemake@output[["pdf"]],
  heatmap_plot,
  base_width = 14,
  base_height = 12,
  units = "cm",
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  heatmap_plot,
  base_width = 14,
  base_height = 12,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")