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

log_message("Starting ATAC / Sertoli / Granulosa overlap analysis script")

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

log_message("Loading ATAC peak cluster table")

ATAC_peaks <- read.csv(
  snakemake@input[["ATAC_peaks"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"cluster" %in% colnames(ATAC_peaks)) {
  log_error("Column 'cluster' not found in ATAC peak table")
}

log_message("Importing ATAC-Sertoli peaks")

ATAC_Sertoli_peaks <- rtracklayer::import(
  snakemake@input[["ATAC_Sertoli_peaks"]]
)

log_message("Importing ATAC-Granulosa biased/specific peaks")

ATAC_Granulosa_peaks <- rtracklayer::import(
  snakemake@input[["ATAC_Granulosa_peaks"]]
)

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
# Load output file
###########################################

common_ATAC_file <- snakemake@output[["common_ATAC"]]

if (is.null(common_ATAC_file)) {
  log_error("Output 'common_ATAC' must be defined as a single output file")
}

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
  category <- rep("ATAC only", length(overlap_sertoli))

  category[overlap_sertoli] <- "ATAC + Sertoli-specific"
  category[!overlap_sertoli & overlap_granulosa] <- "ATAC + Granulosa-biased/specific"

  factor(
    category,
    levels = c(
      "ATAC only",
      "ATAC + Sertoli-specific",
      "ATAC + Granulosa-biased/specific"
    )
  )
}

make_overlap_pie_plot <- function(overlap_category, title, show_legend = TRUE) {
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
      "ATAC only",
      "ATAC + Sertoli-specific",
      "ATAC + Granulosa-biased/specific"
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
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "ATAC only" = "#d267d9",
        "ATAC + Sertoli-specific" = "#316695",
        "ATAC + Granulosa-biased/specific" = "#b741b7"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 13
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = 10)
    )
}

annotate_ATAC_peaks <- function(ATAC_regions, txdb, gene2symbol) {
  if (length(ATAC_regions) == 0) {
    return(data.frame(
      region = character(),
      annotation = character(),
      nearest.gene = character(),
      distanceToTSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      ATAC_regions,
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
    region = make_region_id(ATAC_regions),
    annotation = annotation$annotation,
    nearest.gene = annotation$geneId,
    distanceToTSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

make_overlap_table <- function(
  ATAC_regions,
  overlap_sertoli,
  overlap_granulosa,
  overlap_category,
  cluster_id,
  txdb,
  gene2symbol
) {
  annotation <- annotate_ATAC_peaks(
    ATAC_regions,
    txdb,
    gene2symbol
  )

  output_table <- data.frame(
    cluster = cluster_id,
    chromosome = as.character(GenomicRanges::seqnames(ATAC_regions)),
    start = GenomicRanges::start(ATAC_regions),
    end = GenomicRanges::end(ATAC_regions),
    region = make_region_id(ATAC_regions),
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

  output_table[, c(
    "cluster",
    "chromosome",
    "start",
    "end",
    "region",
    "overlap_ATAC_Sertoli",
    "overlap_ATAC_Granulosa",
    "ATAC_category",
    "annotation",
    "nearest.gene",
    "distanceToTSS"
  )]
}

get_overlap <- function(
  ATAC_regions,
  ATAC_Sertoli_regions,
  ATAC_Granulosa_regions,
  cluster_id,
  txdb,
  gene2symbol,
  title
) {
  overlap_sertoli <- get_overlap_status(
    ATAC_regions,
    ATAC_Sertoli_regions
  )

  overlap_granulosa <- get_overlap_status(
    ATAC_regions,
    ATAC_Granulosa_regions
  )

  overlap_category <- make_overlap_category(
    overlap_sertoli,
    overlap_granulosa
  )

  overlap_table <- make_overlap_table(
    ATAC_regions,
    overlap_sertoli,
    overlap_granulosa,
    overlap_category,
    cluster_id,
    txdb,
    gene2symbol
  )

  plot <- make_overlap_pie_plot(
    overlap_category,
    title,
    show_legend = TRUE
  )

  list(
    plot = plot,
    table = overlap_table
  )
}

###########################################
# Run overlap analysis
###########################################

clusters <- unique(ATAC_peaks$cluster)

results <- lapply(clusters, function(cluster_id) {
  peak_regions <- rownames(
    ATAC_peaks[
      ATAC_peaks$cluster == cluster_id,
      ,
      drop = FALSE
    ]
  )

  ATAC_regions <- GenomicRanges::GRanges(
    peak_regions
  )

  get_overlap(
    ATAC_regions,
    ATAC_Sertoli_peaks,
    ATAC_Granulosa_peaks,
    cluster_id,
    txdb,
    gene2symbol,
    paste0("Cluster ", cluster_id)
  )
})

names(results) <- clusters

plots_with_legend <- lapply(
  results,
  function(result) result$plot
)

shared_legend <- cowplot::get_legend(
  plots_with_legend[[1]] +
    ggplot2::theme(legend.position = "bottom")
)

plots <- lapply(
  plots_with_legend,
  function(plot) {
    plot +
      ggplot2::theme(legend.position = "none")
  }
)

overlap_tables <- lapply(
  results,
  function(result) result$table
)

final_overlap_table <- do.call(
  rbind,
  overlap_tables
)

###########################################
# Write output table
###########################################

write.table(
  final_overlap_table,
  file = snakemake@output[["common_ATAC"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Draw plots
###########################################

plots_grid <- cowplot::plot_grid(
  plotlist = plots,
  labels = names(plots),
  ncol = length(plots),
  align = "h"
)

figure <- cowplot::plot_grid(
  plots_grid,
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