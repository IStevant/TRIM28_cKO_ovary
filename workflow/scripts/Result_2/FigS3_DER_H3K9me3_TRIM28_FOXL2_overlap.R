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

log_message("Starting H3K9me3 / TRIM28 / FOXL2 overlap analysis script")

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
  library("ComplexHeatmap")
  library("gVenn")
  library("cowplot")
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

H3K9me3_table <- read.csv(
  snakemake@input[["H3K9me3_DER"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"x" %in% colnames(H3K9me3_table)) {
  log_error("Column 'x' not found in H3K9me3 peak table")
}

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28"]]
)

log_message("Importing FOXL2 peaks")

FOXL2_peaks <- rtracklayer::import(
  snakemake@input[["FOXL2"]]
)

log_message("Number of H3K9me3 peaks: ", nrow(H3K9me3_table))
log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))
log_message("Number of FOXL2 peaks: ", length(FOXL2_peaks))

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

make_overlap_plots <- function(
  sets_list
) {
  log_message("Computing overlaps")

  overlaps <- gVenn::computeOverlaps(
    sets_list
  )

  log_message("Generating UpSet plot")

  plot_upset <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      gVenn::plotUpSet(overlaps)
    )
  )

  log_message("Generating Venn plot")

  plot_venn <- gVenn::plotVenn(
    overlaps,
    fills = list(
      fill = c(
        H3K9me3 = "#98d344",
        TRIM28 = "#fbc62f",
        FOXL2 = "#b741b7"
      ),
      alpha = 0.45
    ),
    edges = list(
      col = c(
        H3K9me3 = "#98d344",
        TRIM28 = "#fbc62f",
        FOXL2 = "#b741b7"
      ),
      lwd = 2
    )
  )

  cowplot::plot_grid(
    plot_upset,
    plot_venn,
    ncol = 1,
    rel_heights = c(0.65, 0.35)
  )
}

annotate_H3K9me3_peaks <- function(
  H3K9me3_regions,
  cluster_id
) {
  log_message("Annotating H3K9me3 peaks for cluster: ", cluster_id)

  if (length(H3K9me3_regions) == 0) {
    return(data.frame(
      cluster = character(),
      region = character(),
      chromosome = character(),
      start = integer(),
      end = integer(),
      width = integer(),
      annotation = character(),
      nearest_gene = character(),
      distance_to_TSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  peak_anno <- as.data.frame(
    ChIPseeker::annotatePeak(
      H3K9me3_regions,
      genomicAnnotationPriority = c(
        "Promoter",
        "5UTR",
        "Exon",
        "Intron",
        "3UTR",
        "Downstream",
        "Intergenic"
      ),
      tssRegion = c(-promoter, 0),
      TxDb = TxDb,
      level = "gene",
      overlap = "all"
    )
  )

  peak_anno$geneId <- gene2symbol[
    peak_anno$geneId,
    "gene_name"
  ]

  peak_anno$region <- make_region_id(
    H3K9me3_regions
  )

  peak_anno <- peak_anno[, c(
    "region",
    "seqnames",
    "start",
    "end",
    "width",
    "annotation",
    "geneId",
    "distanceToTSS"
  )]

  colnames(peak_anno) <- c(
    "region",
    "chromosome",
    "start",
    "end",
    "width",
    "annotation",
    "nearest_gene",
    "distance_to_TSS"
  )

  peak_anno$cluster <- cluster_id

  peak_anno <- peak_anno[, c(
    "cluster",
    "region",
    "chromosome",
    "start",
    "end",
    "width",
    "annotation",
    "nearest_gene",
    "distance_to_TSS"
  )]

  peak_anno
}

make_overlap_table <- function(
  H3K9me3_regions,
  TRIM28_regions,
  FOXL2_regions,
  cluster_id
) {
  log_message("Computing overlap table for cluster: ", cluster_id)

  overlap_table <- data.frame(
    cluster = cluster_id,
    region = make_region_id(H3K9me3_regions),
    overlap_TRIM28 = ifelse(
      GenomicRanges::countOverlaps(
        H3K9me3_regions,
        TRIM28_regions,
        ignore.strand = TRUE
      ) > 0,
      "YES",
      "NO"
    ),
    overlap_FOXL2 = ifelse(
      GenomicRanges::countOverlaps(
        H3K9me3_regions,
        FOXL2_regions,
        ignore.strand = TRUE
      ) > 0,
      "YES",
      "NO"
    ),
    stringsAsFactors = FALSE
  )

  peak_anno <- annotate_H3K9me3_peaks(
    H3K9me3_regions = H3K9me3_regions,
    cluster_id = cluster_id
  )

  final_table <- merge(
    overlap_table,
    peak_anno,
    by = c("cluster", "region"),
    all.x = TRUE,
    sort = FALSE
  )

  final_table
}

get_cluster_results <- function(
  cluster_id
) {
  log_message("Processing cluster: ", cluster_id)

  peak_regions <- rownames(
    H3K9me3_table[
      H3K9me3_table$x == cluster_id,
      ,
      drop = FALSE
    ]
  )

  H3K9me3_peaks <- GenomicRanges::GRanges(
    peak_regions
  )

  log_message(
    "Number of H3K9me3 peaks in cluster ",
    cluster_id,
    ": ",
    length(H3K9me3_peaks)
  )

  log_message("Extending H3K9me3 peaks by ", distance_to_H3K9me3, " bp")

  H3K9me3_peaks <- H3K9me3_peaks + distance_to_H3K9me3

  peak_list <- list(
    H3K9me3 = H3K9me3_peaks,
    TRIM28 = TRIM28_peaks,
    FOXL2 = FOXL2_peaks
  )

  plot <- make_overlap_plots(
    sets_list = peak_list
  )

  table <- make_overlap_table(
    H3K9me3_regions = H3K9me3_peaks,
    TRIM28_regions = TRIM28_peaks,
    FOXL2_regions = FOXL2_peaks,
    cluster_id = cluster_id
  )

  list(
    plot = plot,
    table = table
  )
}

###########################################
# Run analysis
###########################################

log_message("Preparing cluster-level overlap plots and tables")

clusters <- unique(
  H3K9me3_table$x
)

results <- lapply(
  clusters,
  get_cluster_results
)

names(results) <- clusters

plots <- lapply(
  results,
  function(result) {
    result$plot
  }
)

names(plots) <- clusters

tables <- lapply(
  results,
  function(result) {
    result$table
  }
)

final_table <- do.call(
  rbind,
  tables
)

log_message("Final output table rows: ", nrow(final_table))

###########################################
# Draw final figure
###########################################

figure <- cowplot::plot_grid(
  plotlist = plots,
  labels = names(plots),
  ncol = 2,
  align = "hv"
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 28,
  base_height = 9,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 28,
  base_height = 9,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Writing output table")

write.table(
  final_table,
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

log_message("Analysis completed successfully")