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

log_message("Starting ATAC subcluster annotation overlap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("rtracklayer")
  library("txdbmaker")
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
  library("grid")
})

###########################################
# Load input data
###########################################

log_message("Loading ATAC clustering table")

atac_clustering <- read.csv(
  snakemake@input[["atac_clustering"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

if (!all(c("region", "cluster") %in% colnames(atac_clustering))) {
  log_error("Input table must contain 'region' and 'cluster' columns")
}

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
  annotation_matrix
) {
  ComplexHeatmap::decorate_heatmap_body(
    heatmap_name,
    {
      n_rows <- nrow(annotation_matrix)
      n_cols <- ncol(annotation_matrix)

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

      split_after <- which(
        colnames(annotation_matrix) == "Intergenic"
      )

      if (length(split_after) == 1) {
        x_position <- split_after / n_cols

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

make_overlap_heatmap <- function(
  annotation_matrix
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

  heatmap_name <- "annotation_overlap"

  heatmap <- ComplexHeatmap::Heatmap(
    annotation_matrix,
    name = heatmap_name,
    col = colour_function,
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

  grid::grid.grabExpr({
    ComplexHeatmap::draw(heatmap)

    draw_heatmap_separators(
      heatmap_name = heatmap_name,
      annotation_matrix = annotation_matrix
    )
  })
}

###########################################
# Run analysis
###########################################

log_message("Computing overlap percentages per ATAC subcluster")

clusters <- sort(
  unique(atac_clustering$cluster)
)

annotation_list <- lapply(clusters, function(cluster_id) {
  log_message("Processing cluster: ", cluster_id)

  regions <- GenomicRanges::GRanges(
    atac_clustering[
      atac_clustering$cluster == cluster_id,
      "region"
    ]
  )

  overlap_percent(
    regions = regions,
    promoters = promoters,
    enhancers = enhancers,
    gene_bodies = gene_bodies,
    repeats = repeats,
    ERVs = ERVs,
    L1 = L1
  )
})

annotation_matrix <- do.call(
  rbind,
  annotation_list
)

annotation_matrix <- round(
  annotation_matrix,
  1
)

rownames(annotation_matrix) <- clusters

###########################################
# Draw heatmap
###########################################

log_message("Drawing overlap heatmap")

figure <- make_overlap_heatmap(
  annotation_matrix = annotation_matrix
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 12,
  base_height = 6.5,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 12,
  base_height = 6.5,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")