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

log_message("Starting combined ATAC / chromatin enriched heatmap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("EnrichedHeatmap")
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
  library("grid")
  library("gUtils")
  library("randomcoloR")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

ATAC_bw_folder <- snakemake@params[["ATAC_bw"]]
CR_bw_folder <- snakemake@params[["CR_bw"]]
CR_AB <- snakemake@params[["AB"]]

distance_to_feature <- 1000
top_annotation_quantile <- 0.99
heatmap_quantile <- 0.99

###########################################
# Load input data
###########################################

log_message("Loading samplesheet")

samplesheet <- read.csv(
  file = snakemake@input[["samplesheet"]],
  stringsAsFactors = FALSE
)

samples <- samplesheet$sample

conditions <- unique(
  gsub("_R.*", "", samples)
)

log_message("Conditions detected: ", paste(conditions, collapse = ", "))

log_message("Loading DAR clusters")

DAR_clusters <- read.csv(
  snakemake@input[["sig_DARs"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

log_message("Loading DAR annotation")

DAR_anno <- read.csv(
  snakemake@input[["DAR_anno"]],
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

###########################################
# Find bigWig files
###########################################

log_message("Finding ATAC bigWig files")

ATAC_bw_files <- sapply(conditions, function(condition) {
  matched_files <- list.files(
    path = ATAC_bw_folder,
    pattern = condition,
    full.names = TRUE
  )

  if (length(matched_files) == 0) {
    log_error("No ATAC bigWig found for condition: ", condition)
  }

  if (length(matched_files) > 1) {
    log_error(
      "Multiple ATAC bigWigs found for condition ",
      condition,
      ": ",
      paste(basename(matched_files), collapse = ", ")
    )
  }

  matched_files
})

names(ATAC_bw_files) <- conditions

log_message("Finding chromatin bigWig files")

CR_bw_files <- lapply(CR_AB, function(antibody) {
  files <- sapply(conditions, function(condition) {
    matched_files <- list.files(
      path = CR_bw_folder,
      pattern = paste0(condition, "_", antibody),
      full.names = TRUE
    )

    if (length(matched_files) == 0) {
      log_error(
        "No bigWig found for condition ",
        condition,
        " and antibody ",
        antibody
      )
    }

    if (length(matched_files) > 1) {
      log_error(
        "Multiple bigWigs found for condition ",
        condition,
        " and antibody ",
        antibody,
        ": ",
        paste(basename(matched_files), collapse = ", ")
      )
    }

    matched_files
  })

  names(files) <- conditions
  files
})

names(CR_bw_files) <- CR_AB

bw_files <- c(
  list(ATAC = ATAC_bw_files),
  CR_bw_files
)

###########################################
# Colours
###########################################

log_message("Preparing colours")

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(conditions)
)

names(conditions_color) <- conditions

gradients <- list(
  H3K9me3 = c("#ffffff", "#98d344", "#629b0f"),
  H3K27ac = c("#ffffff", "#79b9e2", "#316695"),
  ATAC = c("#ffffff", "#d267d9", "#8b228d")
)

###########################################
# Functions
###########################################

get_coverage <- function(
  conditions,
  bigwig_files,
  peak_GR,
  distance_to_feature
) {
  coverage_matrices <- lapply(conditions, function(condition) {
    log_message("Importing bigWig for ", condition)

    reads <- rtracklayer::import(
      bigwig_files[[condition]]
    )

    EnrichedHeatmap::normalizeToMatrix(
      reads,
      peak_GR,
      background = 0,
      extend = distance_to_feature,
      w = 50,
      mean_mode = "coverage",
      value_column = "score",
      smooth = TRUE
    )
  })

  names(coverage_matrices) <- conditions

  coverage_matrices
}

get_heatmap_scale <- function(
  matrices,
  quantile_cutoff = 0.98
) {
  scale_value <- quantile(
    unlist(lapply(matrices, as.vector)),
    probs = quantile_cutoff,
    na.rm = TRUE
  )

  scale_value <- ceiling(scale_value)

  if (!is.finite(scale_value) || scale_value == 0) {
    scale_value <- 1
  }

  scale_value
}

get_profile_scale <- function(
  matrices,
  clusters,
  quantile_cutoff = 0.95,
  margin_factor = 1.25
) {
  cluster_levels <- levels(clusters)

  cluster_scales <- lapply(cluster_levels, function(cluster_id) {
    cluster_idx <- clusters == cluster_id

    if (sum(cluster_idx, na.rm = TRUE) == 0) {
      return(NA_real_)
    }

    cluster_profile_values <- unlist(
      lapply(matrices, function(mat) {
        mat <- as.matrix(mat)

        colMeans(
          mat[cluster_idx, , drop = FALSE],
          na.rm = TRUE
        )
      })
    )

    cluster_profile_values <- cluster_profile_values[
      is.finite(cluster_profile_values) &
        cluster_profile_values > 0
    ]

    if (length(cluster_profile_values) == 0) {
      return(NA_real_)
    }

    stats::quantile(
      cluster_profile_values,
      probs = quantile_cutoff,
      na.rm = TRUE
    )
  })

  scale_value <- max(
    unlist(cluster_scales),
    na.rm = TRUE
  )

  scale_value <- scale_value * margin_factor

  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }

  ceiling(scale_value)
}

make_axis_param <- function(
  show_profile_scale,
  profile_scale
) {
  if (show_profile_scale) {
    return(list(
      side = "left",
      facing = "outside",
      at = c(0, round(profile_scale / 2, 2), profile_scale),
      labels = c(0, round(profile_scale / 2, 2), profile_scale),
      labels_rot = 0,
      gp = grid::gpar(
        fontsize = 8,
        col = "#333333"
      )
    ))
  }

  list(
    side = "left",
    facing = "outside",
    at = c(0, round(profile_scale / 2, 2), profile_scale),
    labels = c("", "", ""),
    labels_rot = 0,
    gp = grid::gpar(
      fontsize = 8,
      col = "#333333"
    )
  )
}

build_heatmap <- function(
  matrices,
  heatmap_group,
  conditions,
  heatmap_scale,
  profile_scale,
  distance_to_feature,
  clusters,
  cluster_colours
) {
  log_message("Building heatmap for: ", heatmap_group)

  gradient <- gradients[[heatmap_group]]

  if (is.null(gradient)) {
    gradient <- c("#ffffff", "#999999", "#333333")
  }

  colour_function <- circlize::colorRamp2(
    c(
      0,
      round(heatmap_scale / 2),
      heatmap_scale
    ),
    gradient
  )

  axis_name <- c(
    paste0("-", distance_to_feature / 1000, " kb"),
    "start",
    "end",
    paste0("+", distance_to_feature / 1000, " kb")
  )

  ht_opt$TITLE_PADDING <- grid::unit(2, "mm")

  ht_list <- NULL

  for (condition in conditions) {
    show_profile_scale <- condition == conditions[1]
    show_heatmap_legend <- is.null(ht_list)

    ht <- EnrichedHeatmap::EnrichedHeatmap(
      matrices[[condition]],
      name = if (show_heatmap_legend) {
        heatmap_group
      } else {
        paste0(heatmap_group, "_", condition)
      },
      column_title = paste(heatmap_group, condition),
      column_title_gp = grid::gpar(
        fontface = "bold",
        fill = conditions_color[condition],
        col = "white",
        lwd = 0
      ),
      use_raster = TRUE,
      axis_name = axis_name,
      axis_name_rot = 90,
      raster_quality = 5,
      col = colour_function,
      row_split = clusters,
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        lines = EnrichedHeatmap::anno_enriched(
          ylim = c(0, profile_scale),
          gp = grid::gpar(
            col = cluster_colours,
            lwd = 1.5
          ),
          axis_param = make_axis_param(
            show_profile_scale = show_profile_scale,
            profile_scale = profile_scale
          )
        )
      ),
      show_heatmap_legend = show_heatmap_legend,
      heatmap_legend_param = list(
        title = heatmap_group
      )
    )

    if (is.null(ht_list)) {
      ht_list <- ht
    } else {
      ht_list <- ht_list + ht
    }
  }

  ht_list
}

draw_combined_heatmap <- function(
  ht_list,
  clusters,
  cluster_colours
) {
  cluster_legend <- ComplexHeatmap::Legend(
    at = paste("cluster", levels(clusters)),
    title = "Clusters",
    type = "lines",
    legend_gp = grid::gpar(
      col = cluster_colours
    )
  )

  ComplexHeatmap::draw(
    ht_list,
    ht_gap = grid::unit(3, "mm"),
    annotation_legend_list = list(cluster_legend),
    merge_legend = TRUE,
    heatmap_legend_side = "right"
  )
}

###########################################
# Prepare peak regions
###########################################

log_message("Preparing peak regions")

peaks <- GenomicRanges::GRanges(
  rownames(DAR_clusters)
)

names(peaks) <- rownames(DAR_clusters)
peaks$cluster <- DAR_clusters[, 1]

###########################################
# Compute coverage matrices
###########################################

log_message("Computing coverage matrices")

matrix_list <- lapply(names(bw_files), function(group_name) {
  get_coverage(
    conditions = conditions,
    bigwig_files = bw_files[[group_name]],
    peak_GR = peaks,
    distance_to_feature = distance_to_feature
  )
})

names(matrix_list) <- names(bw_files)

###########################################
# Sub-cluster DAR clusters
###########################################

log_message("Sub-clustering DAR clusters")

big_matrix <- do.call(
  cbind,
  lapply(matrix_list, function(sublist) {
    do.call(
      cbind,
      lapply(sublist, as.matrix)
    )
  })
)

cluster_a_idx <- peaks$cluster == "a"
cluster_b_idx <- peaks$cluster == "b"

clusters_a <- character(0)
clusters_b <- character(0)

if (sum(cluster_a_idx) > 1) {
  matrix_a <- big_matrix[cluster_a_idx, , drop = FALSE]
  matrix_a <- t(scale(t(matrix_a)))
  matrix_a[is.na(matrix_a)] <- 0

  hc_a <- stats::hclust(
    stats::dist(matrix_a),
    method = "ward.D2"
  )

  clusters_a <- paste0(
    "a",
    stats::cutree(hc_a, k = 2)
  )

  clusters_a[clusters_a == "a1"] <- "tmp"
  clusters_a[clusters_a == "a2"] <- "a1"
  clusters_a[clusters_a == "tmp"] <- "a2"
} else if (sum(cluster_a_idx) == 1) {
  clusters_a <- "a1"
}

if (sum(cluster_b_idx) > 1) {
  matrix_b <- big_matrix[cluster_b_idx, , drop = FALSE]
  matrix_b <- t(scale(t(matrix_b)))
  matrix_b[is.na(matrix_b)] <- 0

  hc_b <- stats::hclust(
    stats::dist(matrix_b),
    method = "ward.D2"
  )

  clusters_b <- paste0(
    "b",
    stats::cutree(hc_b, k = 2)
  )
} else if (sum(cluster_b_idx) == 1) {
  clusters_b <- "b1"
}

clusters <- character(length(peaks))
clusters[cluster_a_idx] <- clusters_a
clusters[cluster_b_idx] <- clusters_b

clusters <- factor(
  clusters,
  levels = c("a1", "a2", "b1", "b2")
)

cluster_colours <- randomcoloR::distinctColorPalette(
  length(levels(clusters))
)

names(cluster_colours) <- levels(clusters)

###########################################
# Save clustering table
###########################################

log_message("Writing ATAC/DAR clustering table")

ATAC_DAR_clustering <- data.frame(
  region = names(peaks),
  cluster = as.character(clusters),
  DAR_anno[names(peaks), , drop = FALSE],
  stringsAsFactors = FALSE
)

write.table(
  ATAC_DAR_clustering,
  file = snakemake@output[["atac_clustering"]],
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

###########################################
# Build heatmaps
###########################################

log_message("Building heatmaps")

heatmap_list <- lapply(names(matrix_list), function(group_name) {
  heatmap_scale <- get_heatmap_scale(
    matrices = matrix_list[[group_name]]
  )

  profile_scale <- get_profile_scale(
    matrices = matrix_list[[group_name]],
    clusters = clusters,
    quantile_cutoff = top_annotation_quantile,
    margin_factor = 1.1
  )

  log_message(group_name, " heatmap scale: ", heatmap_scale)
  log_message(group_name, " top annotation profile scale: ", profile_scale)

  build_heatmap(
    matrices = matrix_list[[group_name]],
    heatmap_group = group_name,
    conditions = conditions,
    heatmap_scale = heatmap_scale,
    profile_scale = profile_scale,
    distance_to_feature = distance_to_feature,
    clusters = clusters,
    cluster_colours = cluster_colours
  )
})

ht_list <- NULL

for (ht in heatmap_list) {
  if (is.null(ht_list)) {
    ht_list <- ht
  } else {
    ht_list <- ht_list + ht
  }
}

###########################################
# Add TRIM28 binary annotation
###########################################

log_message("Adding TRIM28 binary annotation")

overlaps <- GenomicRanges::findOverlaps(
  peaks,
  TRIM28_peaks,
  ignore.strand = TRUE
)

binary_vector <- as.integer(
  seq_along(peaks) %in% S4Vectors::queryHits(overlaps)
)

binary_matrix <- matrix(
  binary_vector,
  ncol = 1
)

rownames(binary_matrix) <- names(peaks)
colnames(binary_matrix) <- "TRIM28"

TRIM28_heatmap <- ComplexHeatmap::Heatmap(
  binary_matrix,
  name = "TRIM28",
  col = c(
    "0" = "white",
    "1" = "black"
  ),
  show_row_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "TRIM28",
  width = grid::unit(8, "mm"),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("No", "Yes")
  )
)

ht_list <- ht_list + TRIM28_heatmap

###########################################
# Draw heatmap
###########################################

log_message("Drawing final heatmap")

ATAC_heatmap <- grid::grid.grabExpr(
  draw_combined_heatmap(
    ht_list = ht_list,
    clusters = clusters,
    cluster_colours = cluster_colours
  )
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  snakemake@output[["pdf"]],
  ATAC_heatmap,
  base_width = 5 * length(heatmap_list),
  base_height = 15,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  snakemake@output[["png"]],
  ATAC_heatmap,
  base_width = 5 * length(heatmap_list),
  base_height = 15,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")