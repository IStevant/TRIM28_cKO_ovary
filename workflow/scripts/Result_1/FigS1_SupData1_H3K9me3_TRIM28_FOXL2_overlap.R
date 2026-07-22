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

log_message("Starting overlap analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("cowplot")
  library("gVenn")
  library("grid")
})

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
# Load peak data
###########################################

log_message("Importing peak files")

H3K9me3 <- rtracklayer::import(snakemake@input[["H3K9me3"]])
TRIM28  <- rtracklayer::import(snakemake@input[["TRIM28"]])
FOXL2   <- rtracklayer::import(snakemake@input[["FOXL2"]])

log_message("Extending H3K9me3 peaks by ", distance_to_H3K9me3, " bp")

H3K9me3 <- H3K9me3 + distance_to_H3K9me3

log_message("Number of H3K9me3 peaks: ", length(H3K9me3))
log_message("Number of TRIM28 peaks: ", length(TRIM28))
log_message("Number of FOXL2 peaks: ", length(FOXL2))

peak_GRlist <- list(
  H3K9me3 = H3K9me3,
  TRIM28 = TRIM28,
  FOXL2 = FOXL2
)

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(snakemake@input[["genome"]])

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(snakemake@input[["genome"]])

log_message("Preparing gene ID to symbol mapping")

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

###########################################
# Functions
###########################################

make_overlap_plots <- function(sets_list) {
  log_message("Computing overlaps")

  overlaps <- computeOverlaps(sets_list)

  log_message("Generating UpSet plot")

  plot_upset <- grid::grid.grabExpr(
    ComplexHeatmap::draw(plotUpSet(overlaps))
  )

  log_message("Generating Venn plot")

  plot_venn <- plotVenn(
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

  list(plot_upset, plot_venn)
}

make_region_id <- function(seqnames, start, end) {
  paste0(seqnames, ":", start, "-", end)
}

###########################################
# Plot overlaps
###########################################

log_message("Creating overlap plots")

plots <- make_overlap_plots(peak_GRlist)

figure <- cowplot::plot_grid(
  plotlist = plots,
  ncol = 2
)

###########################################
# Annotate peaks
###########################################

log_message("Annotating H3K9me3 peaks")

peak_anno <- as.data.frame(
  ChIPseeker::annotatePeak(
    H3K9me3,
    genomicAnnotationPriority = c(
      "Promoter", "5UTR", "Exon", "Intron",
      "3UTR", "Downstream", "Intergenic"
    ),
    tssRegion = c(-promoter, 0),
    TxDb = TxDb,
    level = "gene",
    overlap = "all"
  )
)

log_message("Annotated peaks: ", nrow(peak_anno))

peak_anno$geneId <- gene2symbol[peak_anno$geneId, "gene_name"]

peak_anno$region <- make_region_id(
  peak_anno$seqnames,
  peak_anno$start,
  peak_anno$end
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

###########################################
# Create overlap table
###########################################

log_message("Computing peak overlaps")

overlap_table <- data.frame(
  region = make_region_id(
    seqnames(H3K9me3),
    start(H3K9me3),
    end(H3K9me3)
  ),
  overlap_TRIM28 = ifelse(
    countOverlaps(H3K9me3, TRIM28) > 0,
    "YES",
    "NO"
  ),
  overlap_FOXL2 = ifelse(
    countOverlaps(H3K9me3, FOXL2) > 0,
    "YES",
    "NO"
  )
)

log_message("Merging annotation and overlap tables")

final_table <- merge(
  overlap_table,
  peak_anno,
  by = "region",
  all.x = TRUE,
  sort = FALSE
)

log_message("Final table rows: ", nrow(final_table))

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 35,
  base_height = 9,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 35,
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