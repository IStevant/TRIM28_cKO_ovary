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

log_message("Starting single-region genomic track plotting script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("IRanges")
  library("rtracklayer")
  library("txdbmaker")
  library("Gviz")
  library("cowplot")
  library("grid")
})

options(ucscChromosomeNames = FALSE)

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

H3K9me3_bigwig <- snakemake@input[["H3K9me3_bigwig"]]
H3K9me3_domains_file <- snakemake@input[["H3K9me3_domains"]]

TRIM28_peaks_file <- snakemake@input[["TRIM28_peaks"]]
FOXL2_peaks_file <- snakemake@input[["FOXL2_peaks"]]

repeat_masker_file <- snakemake@input[["repeatMasker"]]
CRE_file <- snakemake@input[["CRE"]]

genome_file <- snakemake@input[["genome"]]

region <- snakemake@params[["region"]]

###########################################
# Functions
###########################################

parse_region <- function(region) {
  if (inherits(region, "GRanges")) {
    return(region)
  }

  region <- as.character(region)

  if (!grepl("^.+:[0-9]+-[0-9]+$", region)) {
    log_error(
      "Region must be formatted as 'chr:start-end'. Current value: ",
      region
    )
  }

  chromosome <- sub(":.*", "", region)
  coordinates <- sub(".*:", "", region)

  start_position <- as.numeric(
    sub("-.*", "", coordinates)
  )

  end_position <- as.numeric(
    sub(".*-", "", coordinates)
  )

  GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(
      start = start_position,
      end = end_position
    )
  )
}

make_bigwig_track <- function(
  bigwig_file,
  plot_region
) {
  log_message("Importing bigWig signal")

  coverage <- rtracklayer::import(
    bigwig_file,
    which = plot_region
  )

  if (length(coverage) == 0) {
    log_message("No signal found in selected region")

    coverage <- plot_region
    coverage$score <- 0
  }

  max_score <- round(
    max(coverage$score, na.rm = TRUE),
    digits = 2
  )

  if (!is.finite(max_score) || max_score == 0) {
    max_score <- 1
  }

  Gviz::DataTrack(
    range = coverage,
    ucscChromosomeNames = FALSE,
    type = "hist",
    baseline = 0,
    lwd.baseline = 1,
    chromosome = as.character(
      GenomicRanges::seqnames(plot_region)
    ),
    name = "H3K9me3 signal",
    col.baseline = "#333333",
    col.histogram = 0,
    fill.histogram = "#98d344",
    ylim = c(0, max_score),
    yTicksAt = c(0, max_score),
    rotation.title = 0,
    lwd = 0,
    # alpha = 0.8,
    sizes = 0.4
  )
}

make_annotation_track <- function(
  regions,
  plot_region,
  track_name,
  colour
) {
  regions <- IRanges::subsetByOverlaps(
    regions,
    plot_region,
    ignore.strand = TRUE
  )

  if (length(regions) == 0) {
    regions <- plot_region
    regions$empty <- TRUE
  }

  Gviz::AnnotationTrack(
    range = regions,
    chromosome = as.character(
      GenomicRanges::seqnames(plot_region)
    ),
    start = GenomicRanges::start(regions),
    end = GenomicRanges::end(regions),
    strand = as.character(
      GenomicRanges::strand(regions)
    ),
    name = track_name,
    stacking = "dense",
    col.line = NULL,
    col = 0,
    fill = colour,
    sizes = 0.4,
    rotation.title = 0
  )
}

make_repeatmasker_track <- function(
  repeat_masker_file,
  plot_region
) {
  log_message("Loading RepeatMasker annotation")

  repeats <- read.csv(
    repeat_masker_file,
    sep = "\t",
    header = TRUE
  )

  repeats <- GenomicRanges::makeGRangesFromDataFrame(
    repeats,
    keep.extra.columns = TRUE
  )

  repeats <- repeats[
    grep("ERV", repeats$repFamily, ignore.case = TRUE)
  ]

  repeats <- IRanges::subsetByOverlaps(
    repeats,
    plot_region,
    ignore.strand = TRUE
  )

  if (length(repeats) == 0) {
    repeats <- plot_region
    repeats$repName <- "No_ERV"
  }

  Gviz::AnnotationTrack(
    range = repeats,
    chromosome = as.character(
      GenomicRanges::seqnames(plot_region)
    ),
    start = GenomicRanges::start(repeats),
    end = GenomicRanges::end(repeats),
    strand = "*",
    name = "ERVs",
    stacking = "dense",
    col.line = NULL,
    col = 0,
    fill = "#d23a5b",
    sizes = 0.4,
    rotation.title = 0
  )
}

make_gene_track <- function(
  genome_file,
  plot_region
) {
  log_message("Preparing genome annotation track")

  txdb <- txdbmaker::makeTxDbFromGFF(
    genome_file
  )

  genome_gtf <- rtracklayer::import(
    genome_file
  )


  gene2symbol <- GenomicRanges::mcols(
    genome_gtf
  )[, c("gene_id", "gene_name")]

  gene2symbol <- unique(gene2symbol)

  gene2symbol$gene_name <- paste0(
    gene2symbol$gene_name,
    " "
  )

  rownames(gene2symbol) <- gene2symbol$gene_id

  gene_track <- Gviz::GeneRegionTrack(
    txdb,
    collapseTranscripts = "meta",
    name = "Reference genome",
    chromosome = as.character(
      GenomicRanges::seqnames(plot_region)
    ),
    start = GenomicRanges::start(plot_region),
    end = GenomicRanges::end(plot_region),
    fontface.group = "italic",
    col = 0,
    col.line = NULL,
    fill = "#333333",
    fontcolor.group = "#333333",
    fontsize.group = 16,
    sizes = 0.4,
    rotation.title = 0,
    thinBoxFeature = "UTR"
  )

  ranges(gene_track)$symbol <- gene2symbol[
    ranges(gene_track)$gene,
    "gene_name"
  ]

  gene_track
}

make_gtrack_plot <- function(
  tracks,
  plot_region,
  title
) {
  grid::grid.grabExpr({
    Gviz::plotTracks(
      tracks,
      chromosome = as.character(
        GenomicRanges::seqnames(plot_region)
      ),
      from = GenomicRanges::start(plot_region),
      to = GenomicRanges::end(plot_region),
      transcriptAnnotation = "symbol",
      background.title = "transparent",
      col.border.title = "transparent",
      col.title = "#333333",
      col.axis = "#333333",
      alpha.title = 1,
      alpha.axis = 1,
      sizes = c(
        1,
        0.35,
        0.8,
        0.35,
        0.35,
        0.35,
        0.35
      ),
      cex.axis = 0.7,
      cex.title = 0.8,
      cex.id = 0.5,
      title.width = 1.5,
      main = title,
      cex.main = 1,
      col.main = "#333333"
    )
  })
}

###########################################
# Load input data
###########################################

log_message("Parsing genomic region")

plot_region <- parse_region(region)

log_message(
  "Region to plot: ",
  GenomicRanges::seqnames(plot_region),
  ":",
  GenomicRanges::start(plot_region),
  "-",
  GenomicRanges::end(plot_region)
)

log_message("Importing H3K9me3 domains")

H3K9me3_domains_gr <- rtracklayer::import(
  H3K9me3_domains_file
)

log_message("Importing TRIM28 peaks")

TRIM28_gr <- rtracklayer::import(
  TRIM28_peaks_file
)

log_message("Importing FOXL2 peaks")

FOXL2_gr <- rtracklayer::import(
  FOXL2_peaks_file
)

log_message("Importing CRE annotation")

CRE_gr <- rtracklayer::import(
  CRE_file
)

###########################################
# Prepare tracks
###########################################

log_message("Preparing genomic tracks")

bigwig_track <- make_bigwig_track(
  bigwig_file = H3K9me3_bigwig,
  plot_region = plot_region
)

H3K9me3_domains_track <- make_annotation_track(
  regions = H3K9me3_domains_gr,
  plot_region = plot_region,
  track_name = "H3K9me3 domains",
  colour = "#98d344"
)

gene_track <- make_gene_track(
  genome_file = genome_file,
  plot_region = plot_region
)

TRIM28_track <- make_annotation_track(
  regions = TRIM28_gr,
  plot_region = plot_region,
  track_name = "TRIM28",
  colour = "#fbc62f"
)

FOXL2_track <- make_annotation_track(
  regions = FOXL2_gr,
  plot_region = plot_region,
  track_name = "FOXL2",
  colour = "#b741b7"
)

repeat_track <- make_repeatmasker_track(
  repeat_masker_file = repeat_masker_file,
  plot_region = plot_region
)

CRE_track <- make_annotation_track(
  regions = CRE_gr,
  plot_region = plot_region,
  track_name = "CRE",
  colour = "#3a42d2"
)

###########################################
# Generate figure
###########################################

log_message("Generating genomic track figure")

title <- paste0(
  GenomicRanges::seqnames(plot_region),
  ":",
  GenomicRanges::start(plot_region),
  "-",
  GenomicRanges::end(plot_region)
)

figure <- make_gtrack_plot(
  tracks = list(
    bigwig_track,
    H3K9me3_domains_track,
    gene_track,
    TRIM28_track,
    FOXL2_track,
    CRE_track,
    repeat_track
  ),
  plot_region = plot_region,
  title = title
)

###########################################
# Save figures
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 20,
  base_height = 9,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 20,
  base_height = 9,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")