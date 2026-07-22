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

log_message("Starting H3K9me3 annotation overlap heatmap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("S4Vectors")
  library("rtracklayer")
  library("txdbmaker")
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
  library("grid")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

distance_to_H3K9me3 <- as.numeric(
  snakemake@params[["distance_to_H3K9me3"]]
)

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
  row.names = 1,
  check.names = FALSE
)

if (!"x" %in% colnames(H3K9me3_peaks)) {
  log_error("Column 'x' not found in H3K9me3 DER table")
}

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Extending TRIM28 peaks by ", distance_to_H3K9me3, " bp")

TRIM28_peaks <- TRIM28_peaks + distance_to_H3K9me3

###########################################
# Load RepeatMasker annotation
###########################################

log_message("Loading RepeatMasker annotation")

repeat_masker <- read.csv(
  snakemake@input[["repeatMasker"]],
  sep = "\t",
  header = TRUE
)

repeats <- GenomicRanges::makeGRangesFromDataFrame(
  repeat_masker,
  keep.extra.columns = TRUE
)

ERVs <- repeats[
  grep("ERV", repeats$repFamily, ignore.case = TRUE)
]

L1 <- repeats[
  grep("L1", repeats$repFamily, ignore.case = TRUE)
]

log_message("Number of repeats: ", length(repeats))
log_message("Number of ERV repeats: ", length(ERVs))
log_message("Number of L1 repeats: ", length(L1))

###########################################
# Load enhancer annotation
###########################################

log_message("Importing enhancer annotation")

enhancers <- rtracklayer::import(
  snakemake@input[["enhancers"]]
)

log_message("Number of enhancers: ", length(enhancers))

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

log_message("Building TxDb object")

txdb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

log_message("Extracting promoters and gene bodies")

promoters <- GenomicFeatures::promoters(
  GenomicFeatures::genes(txdb),
  upstream = 2000,
  downstream = 1
)

gene_bodies <- GenomicFeatures::genes(txdb)

log_message("Number of promoters: ", length(promoters))
log_message("Number of gene bodies: ", length(gene_bodies))

###########################################
# Functions
###########################################

get_overlap_binary <- function(
  regions,
  annotation
) {
  GenomicRanges::countOverlaps(
    regions,
    annotation,
    ignore.strand = TRUE
  ) > 0
}

overlap_percent <- function(
  regions,
  promoters,
  enhancers,
  gene_bodies,
  repeats,
  ERVs,
  L1
) {
  feature_order <- c(
    "Enhancer",
    "Promoter",
    "Gene",
    "Intergenic",
    "Repeats",
    "ERV",
    "L1"
  )

  n_regions <- length(regions)

  if (n_regions == 0) {
    return(setNames(
      rep(NA_real_, length(feature_order)),
      feature_order
    ))
  }

  overlap_enhancer <- get_overlap_binary(
    regions,
    enhancers
  )

  overlap_promoter <- get_overlap_binary(
    regions,
    promoters
  )

  overlap_gene <- get_overlap_binary(
    regions,
    gene_bodies
  )

  overlap_repeat <- get_overlap_binary(
    regions,
    repeats
  )

  overlap_ERV <- get_overlap_binary(
    regions,
    ERVs
  )

  overlap_L1 <- get_overlap_binary(
    regions,
    L1
  )

  overlap_intergenic <- !overlap_enhancer &
    !overlap_promoter &
    !overlap_gene

  c(
    Enhancer = 100 * sum(overlap_enhancer) / n_regions,
    Promoter = 100 * sum(overlap_promoter) / n_regions,
    Gene = 100 * sum(overlap_gene) / n_regions,
    Intergenic = 100 * sum(overlap_intergenic) / n_regions,
    Repeats = 100 * sum(overlap_repeat) / n_regions,
    ERV = 100 * sum(overlap_ERV) / n_regions,
    L1 = 100 * sum(overlap_L1) / n_regions
  )
}

draw_heatmap_separators <- function(
  heatmap_name,
  heatmap_matrix,
  thick_column_after
) {
  ComplexHeatmap::decorate_heatmap_body(
    heatmap_name,
    {
      n_rows <- nrow(heatmap_matrix)
      n_cols <- ncol(heatmap_matrix)

      if (n_rows > 1) {
        row_boundaries <- 1 - seq_len(n_rows - 1) / n_rows

        for (y_position in row_boundaries) {
          grid::grid.lines(
            x = grid::unit(c(0, 1), "npc"),
            y = grid::unit(c(y_position, y_position), "npc"),
            gp = grid::gpar(
              col = "white",
              lwd = 2.2
            )
          )
        }
      }

      if (!is.na(thick_column_after)) {
        x_position <- thick_column_after / n_cols

        grid::grid.lines(
          x = grid::unit(c(x_position, x_position), "npc"),
          y = grid::unit(c(0, 1), "npc"),
          gp = grid::gpar(
            col = "white",
            lwd = 3
          )
        )
      }
    }
  )
}

make_heatmap <- function(
  annotation_matrix,
  title
) {
  log_message("Generating heatmap: ", title)

  feature_order <- c(
    "Enhancer",
    "Promoter",
    "Gene",
    "Intergenic",
    "Repeats",
    "ERV",
    "L1"
  )

  annotation_matrix <- annotation_matrix[
    ,
    feature_order,
    drop = FALSE
  ]

  colour_function <- circlize::colorRamp2(
    c(0, 25, 50, 75, 100),
    rev(c(
      "#77BFA3",
      "#98C9A3",
      "#BFD8BD",
      "#DDE7C7",
      "#EDEEC9"
    ))
  )

  heatmap_name <- paste0(
    "overlap_",
    gsub("[^A-Za-z0-9]", "_", title)
  )

  heatmap <- ComplexHeatmap::Heatmap(
    annotation_matrix,
    name = heatmap_name,
    col = colour_function,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    column_names_rot = 45,
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.5
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        paste0(annotation_matrix[i, j], "%"),
        x,
        y,
        gp = grid::gpar(fontsize = 10)
      )
    },
    heatmap_legend_param = list(
      title = "% peaks",
      at = c(0, 25, 50, 75, 100)
    )
  )

  grid::grid.grabExpr(
    {
      ComplexHeatmap::draw(heatmap)

      draw_heatmap_separators(
        heatmap_name = heatmap_name,
        heatmap_matrix = annotation_matrix,
        thick_column_after = which(colnames(annotation_matrix) == "Intergenic")
      )
    }
  )
}

get_cluster_heatmap <- function(
  cluster_id
) {
  log_message("Processing cluster: ", cluster_id)

  peaks <- rownames(
    H3K9me3_peaks[
      H3K9me3_peaks$x == cluster_id,
      ,
      drop = FALSE
    ]
  )

  H3K9me3_regions <- GenomicRanges::GRanges(
    peaks
  )

  hits <- GenomicRanges::findOverlaps(
    H3K9me3_regions,
    TRIM28_peaks,
    ignore.strand = TRUE
  )

  query_hits <- unique(
    S4Vectors::queryHits(hits)
  )

  H3K9me3_only <- H3K9me3_regions[
    -query_hits
  ]

  H3K9me3_ov_TRIM28 <- H3K9me3_regions[
    query_hits
  ]

  log_message("H3K9me3-only peaks: ", length(H3K9me3_only))
  log_message("H3K9me3 peaks overlapping TRIM28: ", length(H3K9me3_ov_TRIM28))

  annotation_matrix <- rbind(
    "H3K9me3 only" = overlap_percent(
      H3K9me3_only,
      promoters,
      enhancers,
      gene_bodies,
      repeats,
      ERVs,
      L1
    ),
    "H3K9me3 + TRIM28" = overlap_percent(
      H3K9me3_ov_TRIM28,
      promoters,
      enhancers,
      gene_bodies,
      repeats,
      ERVs,
      L1
    )
  )

  annotation_matrix <- round(
    annotation_matrix,
    1
  )

  make_heatmap(
    annotation_matrix = annotation_matrix,
    title = paste0("Cluster ", cluster_id)
  )
}

###########################################
# Run analysis
###########################################

log_message("Preparing cluster-level annotation heatmaps")

clusters <- sort(
  unique(H3K9me3_peaks$x)
)

plots <- lapply(
  clusters,
  get_cluster_heatmap
)

names(plots) <- clusters

figure <- cowplot::plot_grid(
  plotlist = plots,
  labels = clusters,
  ncol = 1
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 14,
  base_height = 5.5 * length(plots),
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 14,
  base_height = 5.5 * length(plots),
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")