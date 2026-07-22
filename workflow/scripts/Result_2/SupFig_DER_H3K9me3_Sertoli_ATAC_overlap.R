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

log_message("Starting H3K9me3 and ATAC cell-specific overlap analysis script")

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

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

###########################################
# Load input data
###########################################

log_message("Loading H3K9me3 peak cluster table")

H3K9me3_peaks <- read.csv(
  snakemake@input[["H3K9me3_DER"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"x" %in% colnames(H3K9me3_peaks)) {
  log_error("Column 'x' not found in H3K9me3 peak table")
}

log_message("Importing ATAC Sertoli-specific peaks")

ATAC_Sertoli_peaks <- rtracklayer::import(
  snakemake@input[["ATAC_Sertoli_peaks"]]
)

log_message("Importing ATAC Granulosa-biased/specific peaks")

ATAC_Granulosa_peaks <- rtracklayer::import(
  snakemake@input[["ATAC_Granulosa_peaks"]]
)

###########################################
# Prepare output directory
###########################################

dir.create(
  snakemake@output[["common_H3K9me3"]],
  showWarnings = FALSE,
  recursive = TRUE
)

###########################################
# Load genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

txdb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

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

get_overlap_status <- function(query_regions, subject_regions) {
  hits <- GenomicRanges::findOverlaps(
    query_regions,
    subject_regions,
    ignore.strand = TRUE
  )

  seq_along(query_regions) %in% S4Vectors::queryHits(hits)
}

make_overlap_category <- function(overlap_sertoli, overlap_granulosa) {
  category <- rep("H3K9me3 only", length(overlap_sertoli))

  category[overlap_sertoli] <- "H3K9me3 + ATAC Sertoli-specific"
  category[!overlap_sertoli & overlap_granulosa] <- "H3K9me3 + ATAC Granulosa-biased/specific"

  factor(
    category,
    levels = c(
      "H3K9me3 only",
      "H3K9me3 + ATAC Sertoli-specific",
      "H3K9me3 + ATAC Granulosa-biased/specific"
    )
  )
}

make_overlap_pie_plot <- function(overlap_category, show_legend = TRUE) {
  plot_data <- as.data.frame(
    table(overlap_category),
    stringsAsFactors = FALSE
  )

  colnames(plot_data) <- c("group", "count")

  plot_data <- plot_data[
    plot_data$count > 0,
    ,
    drop = FALSE
  ]

  plot_data$percent <- 100 * plot_data$count / sum(plot_data$count)

  plot_data$label <- paste0(
    round(plot_data$percent, 1),
    "%\n(",
    scales::comma(plot_data$count),
    ")"
  )

  plot_data$group <- factor(
    plot_data$group,
    levels = c(
      "H3K9me3 only",
      "H3K9me3 + ATAC Sertoli-specific",
      "H3K9me3 + ATAC Granulosa-biased/specific"
    )
  )

  legend_position <- if (show_legend) "bottom" else "none"

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
      size = 4.2,
      colour = "black",
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "H3K9me3 only" = "#98d344",
        "H3K9me3 + ATAC Sertoli-specific" = "#316695",
        "H3K9me3 + ATAC Granulosa-biased/specific" = "#b741b7"
      ),
      drop = FALSE
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = 10)
    )
}

annotate_peaks <- function(H3K9me3_regions, txdb, gene2symbol) {
  if (length(H3K9me3_regions) == 0) {
    return(data.frame(
      region = character(),
      annotation = character(),
      nearest_gene = character(),
      distanceToTSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      H3K9me3_regions,
      tssRegion = c(-2000, 0),
      TxDb = txdb,
      level = "gene",
      overlap = "all"
    )
  )

  annotation$geneId <- gene2symbol[
    annotation$geneId,
    "gene_name"
  ]

  data.frame(
    region = make_region_id(H3K9me3_regions),
    annotation = annotation$annotation,
    nearest_gene = annotation$geneId,
    distanceToTSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

write_overlap_table <- function(
  H3K9me3_regions,
  overlap_sertoli,
  overlap_granulosa,
  overlap_category,
  annotation,
  output_file
) {
  output_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(H3K9me3_regions)),
    start = GenomicRanges::start(H3K9me3_regions),
    end = GenomicRanges::end(H3K9me3_regions),
    region = make_region_id(H3K9me3_regions),
    overlap_ATAC_Sertoli = ifelse(overlap_sertoli, "Yes", "No"),
    overlap_ATAC_Granulosa = ifelse(overlap_granulosa, "Yes", "No"),
    ATAC_category = as.character(overlap_category),
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
  ATAC_Sertoli_regions,
  ATAC_Granulosa_regions,
  output_file,
  txdb,
  gene2symbol
) {
  overlap_sertoli <- get_overlap_status(
    H3K9me3_regions,
    ATAC_Sertoli_regions
  )

  overlap_granulosa <- get_overlap_status(
    H3K9me3_regions,
    ATAC_Granulosa_regions
  )

  overlap_category <- make_overlap_category(
    overlap_sertoli,
    overlap_granulosa
  )

  annotation <- annotate_peaks(
    H3K9me3_regions,
    txdb,
    gene2symbol
  )

  write_overlap_table(
    H3K9me3_regions,
    overlap_sertoli,
    overlap_granulosa,
    overlap_category,
    annotation,
    output_file
  )

  make_overlap_pie_plot(
    overlap_category,
    show_legend = TRUE
  )
}

###########################################
# Draw plots
###########################################

clusters <- unique(H3K9me3_peaks$x)

results <- lapply(clusters, function(cluster_id) {
  peak_regions <- rownames(
    H3K9me3_peaks[
      H3K9me3_peaks$x == cluster_id,
      ,
      drop = FALSE
    ]
  )

  H3K9me3_regions <- GenomicRanges::GRanges(
    peak_regions
  )

  output_file <- file.path(
    snakemake@output[["common_H3K9me3"]],
    paste0("DER_H3K9me3_ATAC_cell_specific_", cluster_id, ".tsv")
  )

  get_overlap(
    H3K9me3_regions,
    ATAC_Sertoli_peaks,
    ATAC_Granulosa_peaks,
    output_file,
    txdb,
    gene2symbol
  )
})

names(results) <- clusters

shared_legend <- cowplot::get_legend(
  results[[1]] +
    ggplot2::theme(legend.position = "bottom")
)

plots <- lapply(
  results,
  function(plot) {
    plot +
      ggplot2::theme(legend.position = "none")
  }
)

plot_grid_only <- cowplot::plot_grid(
  plotlist = plots,
  labels = names(plots),
  ncol = length(plots),
  align = "h"
)

figure <- cowplot::plot_grid(
  plot_grid_only,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

###########################################
# Save output files
###########################################

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 9,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 9,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")