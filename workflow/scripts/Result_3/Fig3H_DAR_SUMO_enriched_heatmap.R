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

log_message("Starting ATAC / SUMO enriched heatmap split by TRIM28")

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
  library("randomcoloR")
})

###########################################
# Parameters
###########################################

ATAC_bw_folder <- snakemake@params[["ATAC_bw"]]
ChIP_bw_folder <- snakemake@params[["ChIP_bw"]]
SUMO <- snakemake@params[["SUMO"]]

distance_to_feature <- 1000
top_annotation_quantile <- 0.99
heatmap_quantile <- 0.99

###########################################
# Load data
###########################################

samplesheet <- read.csv(
  snakemake@input[["samplesheet"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

conditions <- unique(
  gsub("_R.*", "", samplesheet$sample)
)

atac_clustering <- read.csv(
  snakemake@input[["atac_clustering"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("region", "cluster") %in% colnames(atac_clustering))) {
  log_error("Input 'atac_clustering' must contain columns: region, cluster")
}

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

###########################################
# Find bigWig files
###########################################

ATAC_bw_files <- sapply(conditions, function(condition) {
  matched_files <- list.files(
    ATAC_bw_folder,
    pattern = condition,
    full.names = TRUE
  )

  if (length(matched_files) != 1) {
    log_error(
      "Expected exactly one ATAC bigWig for ",
      condition,
      ", found: ",
      paste(basename(matched_files), collapse = ", ")
    )
  }

  matched_files
})

names(ATAC_bw_files) <- conditions

SUMO_bw_files <- lapply(SUMO, function(sumo_name) {
  files <- sapply(conditions, function(condition) {
    matched_files <- list.files(
      ChIP_bw_folder,
      pattern = paste0(condition, "_", sumo_name),
      full.names = TRUE
    )

    if (length(matched_files) != 1) {
      log_error(
        "Expected exactly one bigWig for ",
        condition,
        " / ",
        sumo_name,
        ", found: ",
        paste(basename(matched_files), collapse = ", ")
      )
    }

    matched_files
  })

  names(files) <- conditions
  files
})

names(SUMO_bw_files) <- SUMO

bw_files <- c(
  list(ATAC = ATAC_bw_files),
  SUMO_bw_files
)

###########################################
# Colours
###########################################

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(conditions)
)

names(conditions_color) <- conditions

gradients <- list(
  ATAC = c("#ffffff", "#d267d9", "#8b228d"),
  SUMO1 = c("#ffffff", "#f2b6c6", "#b81f5f"),
  SUMO2 = c("#ffffff", "#9fd3e6", "#246a9b")
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
  quantile_cutoff = 0.99
) {
  scale_value <- stats::quantile(
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
  quantile_cutoff = 0.99,
  margin_factor = 1.1
) {
  cluster_scales <- lapply(levels(clusters), function(cluster_id) {
    cluster_idx <- clusters == cluster_id

    if (sum(cluster_idx, na.rm = TRUE) == 0) {
      return(NA_real_)
    }

    profile_values <- unlist(
      lapply(matrices, function(mat) {
        mat <- as.matrix(mat)

        colMeans(
          mat[cluster_idx, , drop = FALSE],
          na.rm = TRUE
        )
      })
    )

    profile_values <- profile_values[
      is.finite(profile_values) &
        profile_values > 0
    ]

    if (length(profile_values) == 0) {
      return(NA_real_)
    }

    stats::quantile(
      profile_values,
      probs = quantile_cutoff,
      na.rm = TRUE
    )
  })

  scale_value <- max(
    unlist(cluster_scales),
    na.rm = TRUE
  ) * margin_factor

  if (!is.finite(scale_value) || scale_value <= 0) {
    scale_value <- 1
  }

  ceiling(scale_value)
}

make_axis_param <- function(
  show_profile_scale,
  profile_scale
) {
  labels <- if (show_profile_scale) {
    c(0, round(profile_scale / 2, 2), profile_scale)
  } else {
    c("", "", "")
  }

  list(
    side = "left",
    facing = "outside",
    at = c(0, round(profile_scale / 2, 2), profile_scale),
    labels = labels,
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

draw_heatmap_page <- function(
  ht_list,
  clusters,
  cluster_colours,
  title
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
    heatmap_legend_side = "right",
    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 14,
      fontface = "bold"
    )
  )
}

make_heatmap_for_subset <- function(
  peaks_subset,
  clusters_subset,
  page_title
) {
  if (length(peaks_subset) == 0) {
    return(
      grid::textGrob(
        paste0(page_title, "\nNo peaks found"),
        gp = grid::gpar(
          fontsize = 14,
          fontface = "bold"
        )
      )
    )
  }

  matrix_list <- lapply(names(bw_files), function(group_name) {
    get_coverage(
      conditions = conditions,
      bigwig_files = bw_files[[group_name]],
      peak_GR = peaks_subset,
      distance_to_feature = distance_to_feature
    )
  })

  names(matrix_list) <- names(bw_files)

  heatmap_list <- lapply(names(matrix_list), function(group_name) {
    heatmap_scale <- get_heatmap_scale(
      matrices = matrix_list[[group_name]],
      quantile_cutoff = heatmap_quantile
    )

    profile_scale <- get_profile_scale(
      matrices = matrix_list[[group_name]],
      clusters = clusters_subset,
      quantile_cutoff = top_annotation_quantile,
      margin_factor = 1.1
    )

    build_heatmap(
      matrices = matrix_list[[group_name]],
      heatmap_group = group_name,
      conditions = conditions,
      heatmap_scale = heatmap_scale,
      profile_scale = profile_scale,
      distance_to_feature = distance_to_feature,
      clusters = clusters_subset,
      cluster_colours = cluster_colours[
        levels(clusters_subset)
      ]
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

  grid::grid.grabExpr(
    draw_heatmap_page(
      ht_list = ht_list,
      clusters = clusters_subset,
      cluster_colours = cluster_colours[
        levels(clusters_subset)
      ],
      title = page_title
    )
  )
}

###########################################
# Prepare peaks and split TRIM28+/-
###########################################

peaks <- GenomicRanges::GRanges(
  atac_clustering$region
)

names(peaks) <- atac_clustering$region

cluster_levels <- c(
  "a1",
  "a2",
  "b1",
  "b2"
)

extra_levels <- setdiff(
  unique(as.character(atac_clustering$cluster)),
  cluster_levels
)

cluster_levels <- c(
  cluster_levels,
  sort(extra_levels)
)

clusters <- factor(
  as.character(atac_clustering$cluster),
  levels = cluster_levels
)

cluster_levels <- levels(clusters)[
  levels(clusters) %in% unique(as.character(clusters))
]

clusters <- factor(
  as.character(clusters),
  levels = cluster_levels
)

peaks$cluster <- clusters

cluster_colours <- randomcoloR::distinctColorPalette(
  length(levels(clusters))
)

names(cluster_colours) <- levels(clusters)

trim28_hits <- GenomicRanges::findOverlaps(
  peaks,
  TRIM28_peaks,
  ignore.strand = TRUE
)

is_trim28_positive <- seq_along(peaks) %in% S4Vectors::queryHits(trim28_hits)

peaks_trim28_positive <- peaks[is_trim28_positive]
clusters_trim28_positive <- droplevels(clusters[is_trim28_positive])

peaks_trim28_negative <- peaks[!is_trim28_positive]
clusters_trim28_negative <- droplevels(clusters[!is_trim28_positive])

log_message("TRIM28+ peaks: ", length(peaks_trim28_positive))
log_message("TRIM28- peaks: ", length(peaks_trim28_negative))

###########################################
# Compute global coverage matrices
###########################################

log_message("Computing global coverage matrices")

global_matrix_list <- lapply(names(bw_files), function(group_name) {
  get_coverage(
    conditions = conditions,
    bigwig_files = bw_files[[group_name]],
    peak_GR = peaks,
    distance_to_feature = distance_to_feature
  )
})

names(global_matrix_list) <- names(bw_files)

global_heatmap_scales <- lapply(names(global_matrix_list), function(group_name) {
  get_heatmap_scale(
    matrices = global_matrix_list[[group_name]],
    quantile_cutoff = heatmap_quantile
  )
})

names(global_heatmap_scales) <- names(global_matrix_list)

global_profile_scales <- lapply(names(global_matrix_list), function(group_name) {
  get_profile_scale(
    matrices = global_matrix_list[[group_name]],
    clusters = clusters,
    quantile_cutoff = top_annotation_quantile,
    margin_factor = 1.1
  )
})

names(global_profile_scales) <- names(global_matrix_list)

###########################################
# Draw heatmaps
###########################################

make_heatmap_for_subset <- function(
  matrix_list,
  row_index,
  clusters_subset,
  page_title
) {
  if (sum(row_index) == 0) {
    return(
      grid::textGrob(
        paste0(page_title, "\nNo peaks found"),
        gp = grid::gpar(
          fontsize = 14,
          fontface = "bold"
        )
      )
    )
  }

  subset_matrix_list <- lapply(matrix_list, function(group_matrices) {
    lapply(group_matrices, function(mat) {
      mat[row_index, , drop = FALSE]
    })
  })

  heatmap_list <- lapply(names(subset_matrix_list), function(group_name) {
    build_heatmap(
      matrices = subset_matrix_list[[group_name]],
      heatmap_group = group_name,
      conditions = conditions,
      heatmap_scale = global_heatmap_scales[[group_name]],
      profile_scale = global_profile_scales[[group_name]],
      distance_to_feature = distance_to_feature,
      clusters = clusters_subset,
      cluster_colours = cluster_colours[
        levels(clusters_subset)
      ]
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

  grid::grid.grabExpr(
    draw_heatmap_page(
      ht_list = ht_list,
      clusters = clusters_subset,
      cluster_colours = cluster_colours[
        levels(clusters_subset)
      ],
      title = page_title
    )
  )
}

TRIM28_positive_heatmap <- make_heatmap_for_subset(
  matrix_list = global_matrix_list,
  row_index = is_trim28_positive,
  clusters_subset = clusters_trim28_positive,
  page_title = "ATAC DARs overlapping TRIM28 peaks"
)

TRIM28_negative_heatmap <- make_heatmap_for_subset(
  matrix_list = global_matrix_list,
  row_index = !is_trim28_positive,
  clusters_subset = clusters_trim28_negative,
  page_title = "ATAC DARs not overlapping TRIM28 peaks"
)

###########################################
# Save PDF with two pages
###########################################

pdf(
  file = snakemake@output[["pdf"]],
  width = 5 * length(bw_files) / 2.54,
  height = 15 / 2.54
)

grid::grid.newpage()
grid::grid.draw(TRIM28_positive_heatmap)

grid::grid.newpage()
grid::grid.draw(TRIM28_negative_heatmap)

dev.off()

###########################################
# Save PNG summary
###########################################

png(
  filename = snakemake@output[["png"]],
  width = 5 * length(bw_files),
  height = 30,
  units = "cm",
  res = 300,
  bg = "white"
)

grid::grid.newpage()
grid::grid.draw(
  cowplot::plot_grid(
    TRIM28_positive_heatmap,
    TRIM28_negative_heatmap,
    ncol = 1
  )
)

dev.off()

log_message("Analysis completed successfully")