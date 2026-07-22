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

log_message("Starting enriched heatmap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
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

bigwig_folder <- snakemake@params[["merged_bigwig_folder"]]
distance_to_feature <- 10000

###########################################
# Load input data
###########################################

log_message("Loading samplesheet")

samplesheet <- read.csv(
  file = snakemake@input[["samplesheet"]],
  stringsAsFactors = FALSE
)

samples <- samplesheet$sample

if (length(samples) == 0) {
  log_error("No samples found in samplesheet")
}

conditions <- unique(
  gsub("_R.*", "", samples)
)

log_message("Conditions detected: ", paste(conditions, collapse = ", "))

log_message("Loading DER clusters")

DER_clusters <- read.csv(
  file = snakemake@input[["clusters"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (nrow(DER_clusters) == 0) {
  log_error("Cluster file has zero rows")
}

###########################################
# Find bigWig files
###########################################

log_message("Finding bigWig files")

bigwig_files <- sapply(conditions, function(condition) {
  matched_files <- list.files(
    path = bigwig_folder,
    pattern = condition,
    full.names = TRUE
  )

  if (length(matched_files) == 0) {
    log_error("No bigWig file found for condition: ", condition)
  }

  if (length(matched_files) > 1) {
    log_error(
      "Multiple bigWig files found for condition ",
      condition,
      ": ",
      paste(basename(matched_files), collapse = ", ")
    )
  }

  matched_files
})

names(bigwig_files) <- conditions

log_message("BigWig files selected: ", paste(basename(bigwig_files), collapse = ", "))

###########################################
# Prepare colours
###########################################

log_message("Preparing colours")

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(conditions)
)

names(conditions_color) <- conditions

heatmap_gradient <- c(
  "#ffffff",
  "#98d344",
  "#629b0f"
)

profile_colour <- "#629b0f"

###########################################
# Functions
###########################################

get_cluster_labels <- function(clusters) {
  cluster_ids <- unique(as.character(clusters))
  cluster_labels <- letters[seq_along(cluster_ids)]
  names(cluster_labels) <- cluster_ids

  cluster_labels
}

get_coverage <- function(
  conditions,
  bigwig_files,
  peak_GR,
  distance_to_feature
) {
  log_message("Computing coverage matrices")

  coverage_matrices <- lapply(conditions, function(condition) {
    log_message("Importing bigWig for condition: ", condition)

    reads <- rtracklayer::import(
      bigwig_files[[condition]]
    )

    log_message("Normalising coverage around peaks for condition: ", condition)

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
  matrices
) {
  scale_value <- max(
    unlist(
      lapply(
        matrices,
        function(mat) {
          colMeans(
            as.matrix(mat),
            na.rm = TRUE
          )
        }
      )
    ),
    na.rm = TRUE
  )

  scale_value <- ceiling(scale_value * 1.1)

  if (!is.finite(scale_value) || scale_value == 0) {
    scale_value <- 1
  }

  scale_value
}

draw_cluster_labels <- function(
  cluster_labels,
  cluster_counts
) {
  for (i in seq_along(cluster_labels)) {
    cluster_id <- names(cluster_labels)[i]

    ComplexHeatmap::decorate_annotation(
      "Clusters",
      slice = i,
      {
        grid::grid.circle(
          x = grid::unit(4, "mm"),
          r = grid::unit(3.5, "mm"),
          gp = grid::gpar(
            fill = "black",
            col = NA
          )
        )

        grid::grid.text(
          cluster_labels[i],
          x = grid::unit(4, "mm"),
          just = "center",
          gp = grid::gpar(
            fontsize = 15,
            col = "white",
            fontface = "bold"
          )
        )

        grid::grid.text(
          paste0(cluster_counts[cluster_id], " regions"),
          x = grid::unit(11, "mm"),
          just = "center",
          rot = 90,
          gp = grid::gpar(
            fontsize = 12,
            col = "black"
          )
        )
      }
    )
  }
}

draw_enriched_heatmap <- function(
  matrices,
  peak_GR,
  conditions,
  conditions_color,
  heatmap_gradient,
  profile_colour,
  max_heatmap_coverage,
  max_profile_coverage,
  distance_to_feature
) {
  log_message("Drawing enriched heatmap")

  colour_function <- circlize::colorRamp2(
    c(
      0,
      round(max_heatmap_coverage / 2),
      max_heatmap_coverage
    ),
    heatmap_gradient
  )

  axis_name <- c(
    paste0("-", distance_to_feature / 1000, " kb"),
    "center",
    paste0("+", distance_to_feature / 1000, " kb")
  )

  cluster_labels <- get_cluster_labels(
    peak_GR$cluster
  )

  peak_GR$cluster <- factor(
    as.character(peak_GR$cluster),
    levels = names(cluster_labels)
  )

  cluster_counts <- table(
    peak_GR$cluster
  )

  cluster_lty <- rep(
    "solid",
    length(cluster_labels)
  )

  if (length(cluster_lty) > 1) {
    cluster_lty[2:length(cluster_lty)] <- "dashed"
  }

  cluster_annotation <- ComplexHeatmap::rowAnnotation(
    Clusters = ComplexHeatmap::anno_empty(
      border = FALSE,
      width = grid::unit(24, "mm")
    ),
    show_annotation_name = FALSE
  )

  ht_opt$TITLE_PADDING <- grid::unit(2, "mm")

  ht_list <- NULL

  for (condition in conditions) {
    show_profile_scale <- condition == conditions[1]
    show_heatmap_legend <- condition == conditions[1]

    axis_param <- if (show_profile_scale) {
      list(
        side = "left",
        facing = "outside",
        gp = grid::gpar(
          fontsize = 8,
          col = "#333333"
        ),
        labels_rot = 0,
        direction = "normal"
      )
    } else {
      list(
        side = "left",
        facing = "outside",
        at = c(
          0,
          round(max_profile_coverage / 2),
          max_profile_coverage
        ),
        labels = c("", "", ""),
        gp = grid::gpar(
          fontsize = 8,
          col = "#333333"
        ),
        labels_rot = 0,
        direction = "normal"
      )
    }

    condition_heatmap <- EnrichedHeatmap::EnrichedHeatmap(
      matrices[[condition]],
      name = paste0("coverage_", condition),
      column_title = condition,
      column_title_gp = grid::gpar(
        col = "black",
        lwd = 0
      ),
      use_raster = TRUE,
      axis_name = axis_name,
      raster_quality = 5,
      col = colour_function,
      axis_name_rot = 90,
      row_split = peak_GR$cluster,
      left_annotation = if (is.null(ht_list)) cluster_annotation else NULL,
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        lines = EnrichedHeatmap::anno_enriched(
          ylim = c(0, max_profile_coverage),
          gp = grid::gpar(
            col = rep(profile_colour, length(cluster_labels)),
            lwd = 1.5,
            lty = cluster_lty
          ),
          axis_param = axis_param
        )
      ),
      show_heatmap_legend = show_heatmap_legend,
      heatmap_legend_param = list(
        title = NULL
      )
    )

    if (is.null(ht_list)) {
      ht_list <- condition_heatmap
    } else {
      ht_list <- ht_list + condition_heatmap
    }
  }

  heatmap <- ComplexHeatmap::draw(
    ht_list,
    ht_gap = grid::unit(3, "mm"),
    merge_legend = TRUE,
    heatmap_legend_side = "right"
  )

  draw_cluster_labels(
    cluster_labels = cluster_labels,
    cluster_counts = cluster_counts
  )

  heatmap
}

###########################################
# Prepare peak regions
###########################################

log_message("Preparing peak regions")

peaks <- GenomicRanges::GRanges(
  rownames(DER_clusters)
)

peaks$cluster <- DER_clusters[, 1]

peaks <- gUtils::gr.mid(peaks)

log_message("Number of peaks: ", length(peaks))

###########################################
# Compute coverage
###########################################

CR_mat <- get_coverage(
  conditions = conditions,
  bigwig_files = bigwig_files,
  peak_GR = peaks,
  distance_to_feature = distance_to_feature
)

###########################################
# Compute shared scales
###########################################

log_message("Computing shared heatmap coverage scale")

max_heatmap_coverage <- get_heatmap_scale(
  matrices = CR_mat,
  quantile_cutoff = 0.98
)

log_message("Computing shared profile coverage scale")

max_profile_coverage <- get_profile_scale(
  matrices = CR_mat
)

log_message("Maximum heatmap coverage scale: ", max_heatmap_coverage)
log_message("Maximum top profile coverage scale: ", max_profile_coverage)

###########################################
# Draw heatmap
###########################################

CR_heatmap <- grid::grid.grabExpr(
  draw_enriched_heatmap(
    matrices = CR_mat,
    peak_GR = peaks,
    conditions = conditions,
    conditions_color = conditions_color,
    heatmap_gradient = heatmap_gradient,
    profile_colour = profile_colour,
    max_heatmap_coverage = max_heatmap_coverage,
    max_profile_coverage = max_profile_coverage,
    distance_to_feature = distance_to_feature
  )
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = CR_heatmap,
  base_width = 9,
  base_height = 15,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = CR_heatmap,
  base_width = 9,
  base_height = 15,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")