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

log_message("Starting genomic annotation overlap analysis script")

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

H3K9me3_peaks <- rtracklayer::import(
  snakemake@input[["H3K9me3_peaks"]]
)

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Number of H3K9me3 peaks: ", length(H3K9me3_peaks))
log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

log_message("Extending TRIM28 peaks by ", distance_to_H3K9me3, " bp")

TRIM28_peaks_extended <- TRIM28_peaks + distance_to_H3K9me3

###########################################
# Load repeat annotation
###########################################

log_message("Loading repeat annotation")

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

log_message("Loading enhancer annotation")

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

TxDb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

log_message("Extracting promoters and gene bodies")

promoters <- GenomicFeatures::promoters(
  GenomicFeatures::genes(TxDb),
  upstream = promoter,
  downstream = 1
)

gene_bodies <- GenomicFeatures::genes(TxDb)

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
  n_regions <- length(regions)

  feature_order <- c(
    "Enhancer",
    "Promoter",
    "Gene",
    "Intergenic",
    "Repeats",
    "ERV",
    "L1"
  )

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

###########################################
# Identify peak groups
###########################################

log_message("Finding overlaps between H3K9me3 and TRIM28 peaks")

hits <- GenomicRanges::findOverlaps(
  H3K9me3_peaks,
  TRIM28_peaks_extended,
  ignore.strand = TRUE
)

query_hits <- S4Vectors::queryHits(hits)
subject_hits <- S4Vectors::subjectHits(hits)

H3K9me3_only <- H3K9me3_peaks[
  -unique(query_hits)
]

TRIM28_only <- TRIM28_peaks_extended[
  -unique(subject_hits)
]

TRIM28_only <- TRIM28_only - distance_to_H3K9me3

both <- GenomicRanges::pintersect(
  H3K9me3_peaks[query_hits],
  TRIM28_peaks_extended[subject_hits]
)

H3K9me3_ov_TRIM28 <- H3K9me3_peaks[
  unique(query_hits)
]

TRIM28_ov_H3K9me3 <- TRIM28_peaks_extended[
  unique(subject_hits)
]

TRIM28_ov_H3K9me3 <- TRIM28_ov_H3K9me3 - distance_to_H3K9me3

log_message("Number of H3K9me3-only peaks: ", length(H3K9me3_only))
log_message("Number of TRIM28-only peaks: ", length(TRIM28_only))
log_message("Number of intersected peak regions: ", length(both))
log_message("Number of H3K9me3 peaks overlapping TRIM28: ", length(H3K9me3_ov_TRIM28))
log_message("Number of TRIM28 peaks overlapping H3K9me3: ", length(TRIM28_ov_H3K9me3))

###########################################
# Generate overlap heatmap
###########################################

log_message("Computing annotation overlap percentages")

hm <- rbind(
  "H3K9me3 only" = overlap_percent(
    regions = H3K9me3_only,
    promoters = promoters,
    enhancers = enhancers,
    gene_bodies = gene_bodies,
    repeats = repeats,
    ERVs = ERVs,
    L1 = L1
  ),
  "TRIM28 + H3K9me3" = overlap_percent(
    regions = TRIM28_ov_H3K9me3,
    promoters = promoters,
    enhancers = enhancers,
    gene_bodies = gene_bodies,
    repeats = repeats,
    ERVs = ERVs,
    L1 = L1
  ),
  "TRIM28 only" = overlap_percent(
    regions = TRIM28_only,
    promoters = promoters,
    enhancers = enhancers,
    gene_bodies = gene_bodies,
    repeats = repeats,
    ERVs = ERVs,
    L1 = L1
  )
)

hm <- round(hm, 1)

feature_order <- c(
  "Enhancer",
  "Promoter",
  "Gene",
  "Intergenic",
  "Repeats",
  "ERV",
  "L1"
)

hm <- hm[
  ,
  feature_order,
  drop = FALSE
]

log_message("Generating annotation overlap heatmap")

col_fun <- circlize::colorRamp2(
  c(0, 25, 50, 75, 100),
  rev(c("#77BFA3", "#98C9A3", "#BFD8BD", "#DDE7C7", "#EDEEC9"))
)

heatmap_name <- "overlap"

ht <- ComplexHeatmap::Heatmap(
  hm,
  name = heatmap_name,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  rect_gp = grid::gpar(
    col = "white",
    lwd = 0.5
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(
      paste0(hm[i, j], "%"),
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

figure <- grid::grid.grabExpr(
  {
    ComplexHeatmap::draw(ht)

    draw_heatmap_separators(
      heatmap_name = heatmap_name,
      heatmap_matrix = hm,
      thick_column_after = which(colnames(hm) == "Intergenic")
    )
  }
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 15,
  base_height = 6,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 15,
  base_height = 6,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")