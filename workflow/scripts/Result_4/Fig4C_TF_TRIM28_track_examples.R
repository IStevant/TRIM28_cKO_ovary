source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting ATAC / TRIM28 / TF genomic track plotting script")

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
# Inputs
###########################################

ATAC_ctrl_bigwig <- snakemake@input[["ATAC_ctrl_bigwig"]]
ATAC_cKO_bigwig <- snakemake@input[["ATAC_cKO_bigwig"]]

ATAC_peaks_file <- snakemake@input[["ATAC_peaks"]]
ATAC_clustering_file <- snakemake@input[["ATAC_clustering"]]

TRIM28_peaks_file <- snakemake@input[["TRIM28_peaks"]]
FOXL2_peaks_file <- snakemake@input[["FOXL2_peaks"]]
NR5A2_peaks_file <- snakemake@input[["NR5A2_peaks"]]
ESR2_peaks_file <- snakemake@input[["ESR2_peaks"]]
RUNX_peaks_file <- snakemake@input[["RUNX_peaks"]]

genome_file <- snakemake@input[["genome"]]

region <- snakemake@params[["region"]]

###########################################
# Functions
###########################################

parse_region <- function(region) {
  region <- as.character(region)

  if (!grepl("^.+:[0-9]+-[0-9]+$", region)) {
    log_error("Region must be formatted as 'chr:start-end'. Current value: ", region)
  }

  chromosome <- sub(":.*", "", region)
  coordinates <- sub(".*:", "", region)

  GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(
      start = as.numeric(sub("-.*", "", coordinates)),
      end = as.numeric(sub(".*-", "", coordinates))
    )
  )
}

import_bigwig_signal <- function(bigwig_file, plot_region) {
  signal <- rtracklayer::import(
    bigwig_file,
    which = plot_region
  )

  if (length(signal) == 0) {
    signal <- plot_region
    signal$score <- 0
  }

  signal
}

make_shared_bigwig_ylim <- function(bigwig_files, plot_region) {
  max_values <- vapply(
    bigwig_files,
    function(bigwig_file) {
      signal <- import_bigwig_signal(bigwig_file, plot_region)
      max(signal$score, na.rm = TRUE)
    },
    numeric(1)
  )

  ylim_max <- round(max(max_values, na.rm = TRUE), 2)

  if (!is.finite(ylim_max) || ylim_max == 0) {
    ylim_max <- 1
  }

  ylim_max
}

make_bigwig_track <- function(
  bigwig_file,
  plot_region,
  track_name,
  ylim_max
) {
  signal <- import_bigwig_signal(
    bigwig_file = bigwig_file,
    plot_region = plot_region
  )

  Gviz::DataTrack(
    range = signal,
    ucscChromosomeNames = FALSE,
    type = "hist",
    baseline = 0,
    lwd.baseline = 1,
    chromosome = as.character(GenomicRanges::seqnames(plot_region)),
    name = track_name,
    col.baseline = "#333333",
    col.histogram = NA,
    fill.histogram = "#a23aa6",
    ylim = c(0, ylim_max),
    yTicksAt = c(0, ylim_max),
    rotation.title = 0,
    lwd = 0,
    sizes = 0.5
  )
}

make_empty_annotation_track <- function(
  plot_region,
  track_name
) {
  empty_region <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(plot_region),
    ranges = IRanges::IRanges(
      start = GenomicRanges::start(plot_region),
      width = 1
    )
  )

  Gviz::AnnotationTrack(
    range = empty_region,
    chromosome = as.character(GenomicRanges::seqnames(plot_region)),
    name = track_name,
    stacking = "dense",
    col = "transparent",
    col.line = "transparent",
    fill = "transparent",
    alpha = 0,
    rotation.title = 0
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
    return(
      make_empty_annotation_track(
        plot_region = plot_region,
        track_name = track_name
      )
    )
  }

  Gviz::AnnotationTrack(
    range = regions,
    chromosome = as.character(
      GenomicRanges::seqnames(plot_region)
    ),
    start = GenomicRanges::start(regions),
    end = GenomicRanges::end(regions),
    strand = "*",
    name = track_name,
    stacking = "dense",
    col.line = NA,
    col = NA,
    fill = colour,
    sizes = 0.4,
    rotation.title = 0
  )
}

make_ATAC_peak_track <- function(
  regions,
  plot_region,
  track_name
) {
  regions <- IRanges::subsetByOverlaps(
    regions,
    plot_region,
    ignore.strand = TRUE
  )

  if (length(regions) == 0) {
    return(NULL)
  }

  Gviz::AnnotationTrack(
    range = regions,
    chromosome = as.character(GenomicRanges::seqnames(plot_region)),
    start = GenomicRanges::start(regions),
    end = GenomicRanges::end(regions),
    strand = "*",
    name = track_name,
    stacking = "dense",
    col.line = NA,
    col = NA,
    fill = "#a23aa6",
    sizes = 0.4,
    rotation.title = 0
  )
}

make_DAR_highlight_track <- function(
  track_list,
  DAR_a,
  DAR_b,
  plot_region
) {
  track_list <- track_list[
    !vapply(track_list, is.null, logical(1))
  ]

  DAR_a <- IRanges::subsetByOverlaps(
    DAR_a,
    plot_region,
    ignore.strand = TRUE
  )

  DAR_b <- IRanges::subsetByOverlaps(
    DAR_b,
    plot_region,
    ignore.strand = TRUE
  )

  highlighted_tracks <- track_list

  if (length(DAR_a) > 0) {
    highlighted_tracks <- Gviz::HighlightTrack(
      trackList = highlighted_tracks,
      range = DAR_a + 100,
      chromosome = as.character(GenomicRanges::seqnames(plot_region)),
      start = GenomicRanges::start(DAR_a),
      end = GenomicRanges::end(DAR_a),
      col = "transparent",
      fill = "#7FB3FF33",
      inBackground = TRUE
    )
  }

  if (length(DAR_b) > 0) {
    highlighted_tracks <- Gviz::HighlightTrack(
      trackList = highlighted_tracks,
      range = DAR_b + 100,
      chromosome = as.character(GenomicRanges::seqnames(plot_region)),
      start = GenomicRanges::start(DAR_b),
      end = GenomicRanges::end(DAR_b),
      col = "transparent",
      fill = "#D94A6A33",
      inBackground = TRUE
    )
  }

  highlighted_tracks
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
  title,
  track_sizes
) {
  if (is.list(tracks) && !inherits(tracks, "HighlightTrack")) {
    tracks <- tracks[
      !vapply(
        tracks,
        is.null,
        logical(1)
      )
    ]
  }

  grid::grid.grabExpr({
    Gviz::plotTracks(
      tracks,
      chromosome = as.character(GenomicRanges::seqnames(plot_region)),
      from = GenomicRanges::start(plot_region),
      to = GenomicRanges::end(plot_region),
      transcriptAnnotation = "symbol",
      background.title = "transparent",
      col.border.title = "transparent",
      col.title = "#333333",
      col.axis = "#333333",
      alpha.title = 1,
      alpha.axis = 1,
      cex.axis = 0.7,
      cex.title = 0.8,
      cex.id = 0.5,
      title.width = 1.7,
      main = title,
      cex.main = 1,
      col.main = "#333333",
      sizes = track_sizes
    )
  })
}

###########################################
# Parse region
###########################################

plot_region <- parse_region(region)

title <- paste0(
  GenomicRanges::seqnames(plot_region),
  ":",
  GenomicRanges::start(plot_region),
  "-",
  GenomicRanges::end(plot_region)
)

###########################################
# Import tracks
###########################################

log_message("Importing ATAC peak files")

ATAC_peaks <- rtracklayer::import(
  ATAC_peaks_file
)

log_message("Importing DAR clustering table")

ATAC_clustering <- read.csv(
  ATAC_clustering_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("region", "cluster") %in% colnames(ATAC_clustering))) {
  log_error("ATAC_clustering must contain 'region' and 'cluster' columns")
}

DAR_regions <- GenomicRanges::GRanges(
  ATAC_clustering$region
)

DAR_regions$cluster <- ATAC_clustering$cluster

DAR_a <- DAR_regions[
  DAR_regions$cluster == "a"
]

DAR_b <- DAR_regions[
  DAR_regions$cluster == "b"
]

log_message("DAR down in cKO / cluster a: ", length(DAR_a))
log_message("DAR up in cKO / cluster b: ", length(DAR_b))

log_message("Importing TRIM28 and TF peak files")

TRIM28_peaks <- rtracklayer::import(TRIM28_peaks_file)
FOXL2_peaks <- rtracklayer::import(FOXL2_peaks_file)
NR5A2_peaks <- rtracklayer::import(NR5A2_peaks_file)
ESR2_peaks <- rtracklayer::import(ESR2_peaks_file)
RUNX_peaks <- rtracklayer::import(RUNX_peaks_file)

log_message("Import done")


###########################################
# Prepare tracks
###########################################

ATAC_ylim <- make_shared_bigwig_ylim(
  bigwig_files = c(
    ATAC_ctrl_bigwig,
    ATAC_cKO_bigwig
  ),
  plot_region = plot_region
)

ATAC_ctrl_track <- make_bigwig_track(
  bigwig_file = ATAC_ctrl_bigwig,
  plot_region = plot_region,
  track_name = "ATAC Ctrl",
  ylim_max = ATAC_ylim
)

ATAC_cKO_track <- make_bigwig_track(
  bigwig_file = ATAC_cKO_bigwig,
  plot_region = plot_region,
  track_name = "ATAC cKO",
  ylim_max = ATAC_ylim
)

ATAC_peak_track <- make_ATAC_peak_track(
  regions = ATAC_peaks,
  plot_region = plot_region,
  track_name = "ATAC peaks"
)

TRIM28_track <- make_annotation_track(
  regions = TRIM28_peaks,
  plot_region = plot_region,
  track_name = "TRIM28",
  colour = "#fbc62f"
)

FOXL2_track <- make_annotation_track(
  regions = FOXL2_peaks,
  plot_region = plot_region,
  track_name = "FOXL2",
  colour = "#e874e8"
)

NR5A2_track <- make_annotation_track(
  regions = NR5A2_peaks,
  plot_region = plot_region,
  track_name = "NR5A2",
  colour = "#329ca3"
)

ESR2_track <- make_annotation_track(
  regions = ESR2_peaks,
  plot_region = plot_region,
  track_name = "ESR2",
  colour = "#d8581c"
)

RUNX_track <- make_annotation_track(
  regions = RUNX_peaks,
  plot_region = plot_region,
  track_name = "RUNX",
  colour = "#86b524"
)

gene_track <- make_gene_track(
  genome_file = genome_file,
  plot_region = plot_region
)

track_list <- list(
  ATAC_ctrl_track,
  ATAC_cKO_track,
  ATAC_peak_track,
  TRIM28_track,
  FOXL2_track,
  NR5A2_track,
  ESR2_track,
  RUNX_track,
  gene_track
)

track_sizes <- c(
  0.3,  # ATAC ctrl
  0.3,  # ATAC cKO
  0.15, # ATAC peaks
  0.15, # TRIM28
  0.15, # FOXL2
  0.15, # NR5A2
  0.15, # ESR2
  0.15, # RUNX
  0.5   # genes
)

highlighted_tracks <- make_DAR_highlight_track(
  track_list = track_list,
  DAR_a = DAR_a,
  DAR_b = DAR_b,
  plot_region = plot_region
)


###########################################
# Generate figure
###########################################

figure <- make_gtrack_plot(
  tracks = highlighted_tracks,
  plot_region = plot_region,
  title = title,
  track_sizes = track_sizes
)

###########################################
# Save outputs
###########################################

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 10,
  base_height = 8,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 10,
  base_height = 8,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")