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

log_message("Starting DAR and TRIM28 overlap analysis script")

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

log_message("Loading DAR peak cluster table")

DAR_peaks <- read.csv(
  snakemake@input[["DAR_peaks"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"cluster" %in% colnames(DAR_peaks)) {
  log_error("Column 'cluster' not found in DAR peak table")
}

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

###########################################
# Prepare output directory
###########################################

log_message("Preparing output directory")

dir.create(
  snakemake@output[["common_ATAC"]],
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

log_message("Building TxDb object")

txdb <- txdbmaker::makeTxDbFromGFF(
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
      "ATAC",
      "ATAC + TRIM28"
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
      "ATAC",
      "ATAC + TRIM28"
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
        "ATAC" = "#d267d9",
        "ATAC + TRIM28" = "#8b228d"
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 11)
    )
}

annotate_ATAC_peaks <- function(
  ATAC_regions,
  TRIM28_status,
  txdb,
  gene2symbol
) {
  if (length(ATAC_regions) == 0) {
    return(data.frame(
      chromosome = character(),
      start = integer(),
      end = integer(),
      region = character(),
      TRIM28 = character(),
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
    chromosome = as.character(GenomicRanges::seqnames(ATAC_regions)),
    start = GenomicRanges::start(ATAC_regions),
    end = GenomicRanges::end(ATAC_regions),
    region = make_region_id(ATAC_regions),
    TRIM28 = TRIM28_status,
    annotation = annotation$annotation,
    nearest.gene = annotation$geneId,
    distanceToTSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

write_ATAC_table <- function(
  common_ATAC,
  not_common_ATAC,
  output_file,
  txdb,
  gene2symbol
) {
  common_table <- annotate_ATAC_peaks(
    ATAC_regions = common_ATAC,
    TRIM28_status = rep("Yes", length(common_ATAC)),
    txdb = txdb,
    gene2symbol = gene2symbol
  )

  not_common_table <- annotate_ATAC_peaks(
    ATAC_regions = not_common_ATAC,
    TRIM28_status = rep("No", length(not_common_ATAC)),
    txdb = txdb,
    gene2symbol = gene2symbol
  )

  output_table <- rbind(
    common_table,
    not_common_table
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
  DAR_regions,
  TRIM28_regions,
  output_file,
  txdb,
  gene2symbol
) {
  hits <- GenomicRanges::findOverlaps(
    DAR_regions,
    TRIM28_regions,
    ignore.strand = TRUE
  )

  query_hits <- unique(
    S4Vectors::queryHits(hits)
  )

  common_ATAC <- DAR_regions[
    query_hits
  ]

  not_common_ATAC <- DAR_regions[
    setdiff(
      seq_along(DAR_regions),
      query_hits
    )
  ]

  n_total <- length(DAR_regions)
  n_overlap <- length(common_ATAC)

  log_message("Total DAR/ATAC peaks in cluster: ", n_total)
  log_message("DAR/ATAC peaks overlapping TRIM28: ", n_overlap)
  log_message("DAR/ATAC peaks not overlapping TRIM28: ", n_total - n_overlap)

  write_ATAC_table(
    common_ATAC = common_ATAC,
    not_common_ATAC = not_common_ATAC,
    output_file = output_file,
    txdb = txdb,
    gene2symbol = gene2symbol
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

clusters <- unique(DAR_peaks$cluster)

plots <- lapply(clusters, function(cluster_id) {
  log_message("Processing cluster: ", cluster_id)

  peak_regions <- rownames(
    DAR_peaks[
      DAR_peaks$cluster == cluster_id,
      ,
      drop = FALSE
    ]
  )

  DAR_regions <- GenomicRanges::GRanges(
    peak_regions
  )

  output_file <- file.path(
    snakemake@output[["common_ATAC"]],
    paste0("DAR_TRIM28_", cluster_id, ".tsv")
  )

  get_overlap(
    DAR_regions = DAR_regions,
    TRIM28_regions = TRIM28_peaks,
    output_file = output_file,
    txdb = txdb,
    gene2symbol = gene2symbol
  )
})

names(plots) <- clusters

figure <- cowplot::plot_grid(
  plotlist = plots,
  labels = names(plots),
  ncol = length(plots),
  align = "h"
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