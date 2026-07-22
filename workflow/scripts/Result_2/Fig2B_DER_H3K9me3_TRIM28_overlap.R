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

log_message("Starting H3K9me3 cluster overlap analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("txdbmaker")
  library("ChIPseeker")
  library("cowplot")
  library("ggplot2")
  library("scales")
  library("grid")
})

###########################################
# ChIPseeker options
###########################################

log_message("Setting ChIPseeker options")

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

promoter <- as.numeric(snakemake@params[["promoter"]])
distance_to_H3K9me3 <- as.numeric(snakemake@params[["distance_to_H3K9me3"]])

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}

if (is.na(distance_to_H3K9me3)) {
  log_error("Parameter 'distance_to_H3K9me3' must be numeric")
}

###########################################
# Load input data
###########################################

log_message("Loading H3K9me3 peak cluster table")

H3K9me3_peaks <- read.csv(
  snakemake@input[["H3K9me3_DER"]],
  header = TRUE,
  row.names = 1
)

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Number of H3K9me3 peaks in table: ", nrow(H3K9me3_peaks))
log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

log_message("Extending TRIM28 peaks by ", distance_to_H3K9me3, " bp")

TRIM28_peaks_extended <- TRIM28_peaks + distance_to_H3K9me3

###########################################
# Prepare output directory
###########################################

log_message("Preparing output directory")

dir.create(
  snakemake@output[["common_H3K9me3"]],
  showWarnings = FALSE,
  recursive = TRUE
)

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

log_message("Preparing gene ID to symbol mapping")

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

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

make_overlap_pie_plot <- function(
  n_total,
  n_overlap
) {
  n_not_overlap <- n_total - n_overlap

  plot_data <- data.frame(
    group = c(
      "H3K9me3",
      "H3K9me3 + TRIM28"
    ),
    count = c(
      n_not_overlap,
      n_overlap
    )
  )

  plot_data$percent <- if (sum(plot_data$count) == 0) {
    0
  } else {
    100 * plot_data$count / sum(plot_data$count)
  }

  plot_data$label <- paste0(
    round(plot_data$percent, 1),
    "%\n(",
    scales::comma(plot_data$count),
    ")"
  )

  plot_data$group <- factor(
    plot_data$group,
    levels = c(
      "H3K9me3",
      "H3K9me3 + TRIM28"
    )
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = "",
      y = percent,
      fill = group
    )
  ) +
    ggplot2::geom_bar(
      width = 1,
      stat = "identity",
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::coord_polar(
      theta = "y",
      start = 0
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 4.5,
      colour = "black",
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "H3K9me3" = "#c2e48f",
        "H3K9me3 + TRIM28" = "#629b0f"
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      )
    )
}

annotate_common_peaks <- function(
  common_H3K9me3,
  TxDb,
  gene2symbol
) {
  if (length(common_H3K9me3) == 0) {
    return(data.frame(
      region = character(),
      annotation = character(),
      nearest_gene = character(),
      distance_to_TSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      common_H3K9me3,
      tssRegion = c(-promoter, 0),
      TxDb = TxDb,
      level = "gene",
      overlap = "all"
    )
  )

  annotation$geneId <- gene2symbol[
    annotation$geneId,
    "gene_name"
  ]

  data.frame(
    region = make_region_id(common_H3K9me3),
    annotation = annotation$annotation,
    nearest_gene = annotation$geneId,
    distance_to_TSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

write_common_peak_table <- function(
  common_H3K9me3,
  annotation,
  output_file
) {
  output_table <- data.frame(
    chromosome = as.character(
      GenomicRanges::seqnames(common_H3K9me3)
    ),
    start = GenomicRanges::start(common_H3K9me3),
    end = GenomicRanges::end(common_H3K9me3),
    region = make_region_id(common_H3K9me3),
    stringsAsFactors = FALSE
  )

  output_table <- merge(
    output_table,
    annotation,
    by = "region",
    all.x = TRUE,
    sort = FALSE
  )

  write.table(
    output_table,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

get_overlap <- function(
  H3K9me3_regions,
  TRIM28_regions,
  output_file,
  TxDb,
  gene2symbol
) {
  hits <- GenomicRanges::findOverlaps(
    H3K9me3_regions,
    TRIM28_regions,
    ignore.strand = TRUE
  )

  query_hits <- unique(
    S4Vectors::queryHits(hits)
  )

  common_H3K9me3 <- H3K9me3_regions[
    query_hits
  ]

  n_total <- length(H3K9me3_regions)
  n_overlap <- length(common_H3K9me3)

  log_message("Total H3K9me3 peaks in cluster: ", n_total)
  log_message("H3K9me3 peaks overlapping TRIM28: ", n_overlap)
  log_message("H3K9me3 peaks not overlapping TRIM28: ", n_total - n_overlap)

  annotation <- annotate_common_peaks(
    common_H3K9me3 = common_H3K9me3,
    TxDb = TxDb,
    gene2symbol = gene2symbol
  )

  write_common_peak_table(
    common_H3K9me3 = common_H3K9me3,
    annotation = annotation,
    output_file = output_file
  )

  make_overlap_pie_plot(
    n_total = n_total,
    n_overlap = n_overlap
  )
}

###########################################
# Draw plots
###########################################

log_message("Preparing cluster-level overlap plots")

clusters <- unique(H3K9me3_peaks$x)

plots <- lapply(clusters, function(cluster_id) {
  log_message("Processing cluster: ", cluster_id)

  peak_regions <- rownames(
    H3K9me3_peaks[
      H3K9me3_peaks$x == cluster_id,
      ,
      drop = FALSE
    ]
  )

  output_file <- paste0(
    snakemake@output[["common_H3K9me3"]],
    "/DER_H3K9me3_TRIM28_",
    cluster_id,
    ".csv"
  )

  H3K9me3_regions <- GenomicRanges::GRanges(
    peak_regions
  )

  get_overlap(
    H3K9me3_regions = H3K9me3_regions,
    TRIM28_regions = TRIM28_peaks_extended,
    output_file = output_file,
    TxDb = TxDb,
    gene2symbol = gene2symbol
  )
})

names(plots) <- clusters

###########################################
# Combine plots with shared legend
###########################################

legend <- cowplot::get_legend(
  plots[[1]] +
    ggplot2::theme(
      legend.position = "bottom"
    )
)

plots_no_legend <- lapply(
  plots,
  function(plot) {
    plot +
      ggplot2::theme(
        legend.position = "none"
      )
  }
)

plot_panel <- cowplot::plot_grid(
  plotlist = plots_no_legend,
  labels = names(plots),
  ncol = length(plots),
  align = "h"
)

figure <- cowplot::plot_grid(
  plot_panel,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 15,
  base_height = 8,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 15,
  base_height = 8,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")