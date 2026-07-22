source(".Rprofile")

###########################################
# Logging functions
###########################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}

log_message("Starting H3K9me3 and TRIM28 overlap analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("cowplot")
  library("gVenn")
  library("grid")
  library("gridExtra")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

distance_to_H3K9me3 <- as.numeric(snakemake@params[["distance_to_H3K9me3"]])

if (is.na(distance_to_H3K9me3)) {
  log_error("Parameter 'distance_to_H3K9me3' must be numeric")
}

###########################################
# Load peak data
###########################################

log_message("Importing peak files")

H3K9me3_peaks <- rtracklayer::import(snakemake@input[["H3K9me3_peaks"]])
TRIM28_peaks <- rtracklayer::import(snakemake@input[["TRIM28_peaks"]])

log_message("Number of H3K9me3 peaks: ", length(H3K9me3_peaks))
log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

log_message("Extending TRIM28 peaks by ", distance_to_H3K9me3, " bp")

TRIM28_peaks_extended <- TRIM28_peaks + distance_to_H3K9me3

###########################################
# Functions
###########################################

make_region_id <- function(regions) {
  paste0(
    GenomicRanges::seqnames(regions),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

make_region_table <- function(regions) {
  data.frame(
    chromosome = as.character(GenomicRanges::seqnames(regions)),
    start = GenomicRanges::start(regions),
    end = GenomicRanges::end(regions),
    region = make_region_id(regions),
    stringsAsFactors = FALSE
  )
}

write_region_table <- function(regions, output_file) {
  write.table(
    make_region_table(regions),
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

make_gvenn_plot <- function(sets_list) {
  overlaps <- gVenn::computeOverlaps(sets_list)

  gVenn::plotVenn(
    overlaps,
    fills = list(
      fill = c("#74cfb1", "#ebb655"),
      alpha = 0.45
    ),
    edges = list(
      col = c("#74cfb1", "#ebb655"),
      lwd = 2
    ),
    quantities = list(fontsize = 10)
  )
}

make_summary_table <- function(summary_stats) {
  table_data <- data.frame(
    Category = c(
      "H3K9me3 total",
      "TRIM28 total",
      "H3K9me3 overlapping TRIM28",
      "TRIM28 overlapping H3K9me3",
      "H3K9me3 not overlapping TRIM28",
      "TRIM28 not overlapping H3K9me3"
    ),
    Count = c(
      summary_stats$n_h3k9me3_total,
      summary_stats$n_trim28_total,
      summary_stats$n_h3k9me3_overlap,
      summary_stats$n_trim28_overlap,
      summary_stats$n_h3k9me3_not_overlap,
      summary_stats$n_trim28_not_overlap
    ),
    Percent = c(
      "100%",
      "100%",
      paste0(summary_stats$pct_h3k9me3_overlap, "%"),
      paste0(summary_stats$pct_trim28_overlap, "%"),
      paste0(summary_stats$pct_h3k9me3_not_overlap, "%"),
      paste0(summary_stats$pct_trim28_not_overlap, "%")
    )
  )

  gridExtra::tableGrob(
    table_data,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 9,
      core = list(
        fg_params = list(hjust = 0, x = 0.05)
      ),
      colhead = list(
        fg_params = list(fontface = "bold")
      )
    )
  )
}

make_main_title <- function(distance_to_H3K9me3) {
  distance_kb <- format(
    distance_to_H3K9me3 / 1000,
    trim = TRUE,
    scientific = FALSE
  )

  paste0(
    "Overlap of H3K9me3 and TRIM28 (+/- ",
    distance_kb,
    " kb)"
  )
}

run_overlap_analysis <- function(
  H3K9me3_peaks,
  TRIM28_peaks_extended,
  common_h3k9me3_file,
  common_trim28_file,
  distance_to_H3K9me3
) {
  log_message("Finding overlaps between H3K9me3 and extended TRIM28 peaks")

  hits <- GenomicRanges::findOverlaps(
    H3K9me3_peaks,
    TRIM28_peaks_extended,
    ignore.strand = TRUE
  )

  H3K9me3_hits <- unique(S4Vectors::queryHits(hits))
  TRIM28_hits <- unique(S4Vectors::subjectHits(hits))

  common_h3k9me3 <- H3K9me3_peaks[H3K9me3_hits]
  common_trim28_extended <- TRIM28_peaks_extended[TRIM28_hits]
  common_trim28 <- common_trim28_extended - distance_to_H3K9me3

  n_h3k9me3_total <- length(H3K9me3_peaks)
  n_trim28_total <- length(TRIM28_peaks_extended)

  n_h3k9me3_overlap <- length(common_h3k9me3)
  n_trim28_overlap <- length(common_trim28)

  n_h3k9me3_not_overlap <- n_h3k9me3_total - n_h3k9me3_overlap
  n_trim28_not_overlap <- n_trim28_total - n_trim28_overlap

  pct_h3k9me3_overlap <- round(
    100 * n_h3k9me3_overlap / n_h3k9me3_total,
    1
  )

  pct_trim28_overlap <- round(
    100 * n_trim28_overlap / n_trim28_total,
    1
  )

  pct_h3k9me3_not_overlap <- round(
    100 * n_h3k9me3_not_overlap / n_h3k9me3_total,
    1
  )

  pct_trim28_not_overlap <- round(
    100 * n_trim28_not_overlap / n_trim28_total,
    1
  )

  summary_stats <- list(
    n_h3k9me3_total = n_h3k9me3_total,
    n_trim28_total = n_trim28_total,
    n_h3k9me3_overlap = n_h3k9me3_overlap,
    n_trim28_overlap = n_trim28_overlap,
    n_h3k9me3_not_overlap = n_h3k9me3_not_overlap,
    n_trim28_not_overlap = n_trim28_not_overlap,
    pct_h3k9me3_overlap = pct_h3k9me3_overlap,
    pct_trim28_overlap = pct_trim28_overlap,
    pct_h3k9me3_not_overlap = pct_h3k9me3_not_overlap,
    pct_trim28_not_overlap = pct_trim28_not_overlap
  )

  log_message("H3K9me3 total peaks: ", n_h3k9me3_total)
  log_message("TRIM28 total peaks: ", n_trim28_total)
  log_message("H3K9me3 peaks overlapping TRIM28: ", n_h3k9me3_overlap)
  log_message("TRIM28 peaks overlapping H3K9me3: ", n_trim28_overlap)
  log_message("H3K9me3 peaks not overlapping TRIM28: ", n_h3k9me3_not_overlap)
  log_message("TRIM28 peaks not overlapping H3K9me3: ", n_trim28_not_overlap)

  log_message("Writing common H3K9me3 peaks")

  write_region_table(
    common_h3k9me3,
    common_h3k9me3_file
  )

  log_message("Writing common TRIM28 peaks")

  write_region_table(
    common_trim28,
    common_trim28_file
  )

  log_message("Generating gVenn plot and summary table")

  venn_plot <- make_gvenn_plot(
    list(
      H3K9me3 = H3K9me3_peaks,
      TRIM28 = TRIM28_peaks_extended
    )
  )

  summary_table <- make_summary_table(summary_stats)

  combined_plot <- cowplot::plot_grid(
    venn_plot,
    summary_table,
    ncol = 2,
    rel_widths = c(1.2, 1)
  )

  main_title <- make_main_title(distance_to_H3K9me3)

  cowplot::plot_grid(
    grid::textGrob(
      main_title,
      gp = grid::gpar(
        fontsize = 16,
        fontface = "bold"
      )
    ),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.15, 1)
  )
}

###########################################
# Run analysis
###########################################

venn_plot <- run_overlap_analysis(
  H3K9me3_peaks = H3K9me3_peaks,
  TRIM28_peaks_extended = TRIM28_peaks_extended,
  common_h3k9me3_file = snakemake@output[["common_h3k9me3"]],
  common_trim28_file = snakemake@output[["common_trim28"]],
  distance_to_H3K9me3 = distance_to_H3K9me3
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = venn_plot,
  base_width = 18,
  base_height = 7,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = venn_plot,
  base_width = 18,
  base_height = 7,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")