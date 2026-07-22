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

log_message("Starting AB / SUMO enriched heatmap split by TRIM28")

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

log_message("Loading parameters")

AB <- snakemake@wildcards[["AB"]]

CR_bigwig_folder <- snakemake@params[["CR_bw"]]
ChIP_bigwig_folder <- snakemake@params[["ChIP_bw"]]
SUMO <- snakemake@params[["SUMO"]]

distance_to_feature <- 1000
top_annotation_quantile <- 0.99
heatmap_quantile <- 0.99


###########################################
# Load data
###########################################

log_message("Loading samplesheet")

samplesheet <- read.csv(
  snakemake@input[["samplesheet"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

samples <- samplesheet$sample

sample_conditions <- unique(
  gsub("_R.*", "", samples)
)

base_conditions <- unique(
  sub(paste0("_", AB, "$"), "", sample_conditions)
)

conditions <- sample_conditions

log_message("AB conditions detected: ", paste(sample_conditions, collapse = ", "))
log_message("Base conditions detected for SUMO: ", paste(base_conditions, collapse = ", "))

log_message("Loading DER cluster table")

DER_clusters <- read.csv(
  snakemake@input[["clusters"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (nrow(DER_clusters) == 0) {
  log_error("Cluster file has zero rows")
}

log_message("Number of DER peaks: ", nrow(DER_clusters))

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

###########################################
# Find bigWig files
###########################################

find_bigwig_files <- function(
  conditions,
  bigwig_folder,
  mark
) {
  log_message("Finding bigWig files for: ", mark)

  bigwig_files <- sapply(conditions, function(condition) {
    matched_files <- list.files(
      path = bigwig_folder,
      pattern = condition,
      full.names = TRUE
    )

    matched_files <- matched_files[
      grepl(mark, basename(matched_files))
    ]

    if (length(matched_files) == 0) {
      log_error(
        "No bigWig file found for condition ",
        condition,
        " and mark ",
        mark
      )
    }

    if (length(matched_files) > 1) {
      log_error(
        "Multiple bigWig files found for condition ",
        condition,
        " and mark ",
        mark,
        ": ",
        paste(basename(matched_files), collapse = ", ")
      )
    }

    matched_files
  })

  names(bigwig_files) <- conditions
  bigwig_files
}

AB_bw_files <- find_bigwig_files(
  conditions = sample_conditions,
  bigwig_folder = CR_bigwig_folder,
  mark = AB
)

SUMO_bw_files <- lapply(SUMO, function(sumo_name) {
  sumo_files <- find_bigwig_files(
    conditions = base_conditions,
    bigwig_folder = ChIP_bigwig_folder,
    mark = sumo_name
  )

  names(sumo_files) <- sample_conditions

  sumo_files
})

names(SUMO_bw_files) <- SUMO

bw_files <- c(
  setNames(list(AB_bw_files), AB),
  SUMO_bw_files
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
  SUMO1 = c("#ffffff", "#f2b6c6", "#b81f5f"),
  SUMO2 = c("#ffffff", "#9fd3e6", "#246a9b")
)

gradients[[AB]] <- c(
  "#ffffff",
  "#98d344",
  "#629b0f"
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

  log_message("Computing coverage matrices")

  coverage_matrices <- lapply(conditions, function(condition) {

    log_message("Importing bigWig for condition: ", condition)

    reads <- rtracklayer::import(
      bigwig_files[[condition]]
    )

    log_message(
      "Normalising signal around regions for condition: ",
      condition
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

  return(coverage_matrices)
}


get_heatmap_scale <- function(
  matrices,
  quantile_cutoff = 0.99
) {
  log_message("Calculating heatmap colour scale")

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
  log_message("Calculating top annotation profile scale")

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
  log_message("Building heatmap for signal group: ", heatmap_group)

  gradient <- gradients[[heatmap_group]]

  if (is.null(gradient)) {
    log_message("No predefined gradient found for ", heatmap_group, "; using grey scale")
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
  log_message("Drawing heatmap page: ", title)

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
  matrix_list,
  row_index,
  clusters_subset,
  page_title,
  global_heatmap_scales,
  global_profile_scales
) {
  log_message("Preparing subset heatmap: ", page_title)

  if (sum(row_index) == 0) {
    log_message("No regions found for subset: ", page_title)

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

  log_message("Number of regions in subset: ", sum(row_index))

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

###########################################
# Prepare DER peaks and clustering
###########################################

log_message("Preparing DER regions and cluster labels")

peaks <- GenomicRanges::GRanges(
  rownames(DER_clusters)
)

names(peaks) <- rownames(DER_clusters)

clusters <- factor(
  DER_clusters[, 1]
)

peaks$cluster <- clusters

log_message("Clusters detected: ", paste(levels(clusters), collapse = ", "))

cluster_colours <- randomcoloR::distinctColorPalette(
  length(levels(clusters))
)

names(cluster_colours) <- levels(clusters)

###########################################
# Split TRIM28+ / TRIM28-
###########################################

log_message("Splitting DER peaks by TRIM28 overlap")

trim28_hits <- GenomicRanges::findOverlaps(
  peaks,
  TRIM28_peaks,
  ignore.strand = TRUE
)

is_trim28_positive <- seq_along(peaks) %in% S4Vectors::queryHits(trim28_hits)

clusters_trim28_positive <- droplevels(
  clusters[is_trim28_positive]
)

clusters_trim28_negative <- droplevels(
  clusters[!is_trim28_positive]
)

log_message("TRIM28+ DER peaks: ", sum(is_trim28_positive))
log_message("TRIM28- DER peaks: ", sum(!is_trim28_positive))

###########################################
# Compute global coverage matrices
###########################################

log_message("Computing global coverage matrices for all DER peaks")

global_matrix_list <- lapply(names(bw_files), function(group_name) {
  log_message("Processing signal group: ", group_name)

  get_coverage(
    conditions = conditions,
    bigwig_files = bw_files[[group_name]],
    peak_GR = peaks,
    distance_to_feature = distance_to_feature
  )
})

names(global_matrix_list) <- names(bw_files)

log_message("Computing global heatmap scales")

global_heatmap_scales <- lapply(names(global_matrix_list), function(group_name) {
  scale <- get_heatmap_scale(
    matrices = global_matrix_list[[group_name]],
    quantile_cutoff = heatmap_quantile
  )

  log_message(group_name, " heatmap scale: ", scale)

  scale
})

names(global_heatmap_scales) <- names(global_matrix_list)

log_message("Computing global top profile scales")

global_profile_scales <- lapply(names(global_matrix_list), function(group_name) {
  scale <- get_profile_scale(
    matrices = global_matrix_list[[group_name]],
    clusters = clusters,
    quantile_cutoff = top_annotation_quantile,
    margin_factor = 1.1
  )

  log_message(group_name, " profile scale: ", scale)

  scale
})

names(global_profile_scales) <- names(global_matrix_list)

###########################################
# Draw heatmaps
###########################################

log_message("Drawing TRIM28+ heatmap")

TRIM28_positive_heatmap <- make_heatmap_for_subset(
  matrix_list = global_matrix_list,
  row_index = is_trim28_positive,
  clusters_subset = clusters_trim28_positive,
  page_title = paste0(AB, " DER peaks overlapping TRIM28 peaks"),
  global_heatmap_scales = global_heatmap_scales,
  global_profile_scales = global_profile_scales
)

log_message("Drawing TRIM28- heatmap")

TRIM28_negative_heatmap <- make_heatmap_for_subset(
  matrix_list = global_matrix_list,
  row_index = !is_trim28_positive,
  clusters_subset = clusters_trim28_negative,
  page_title = paste0(AB, " DER peaks not overlapping TRIM28 peaks"),
  global_heatmap_scales = global_heatmap_scales,
  global_profile_scales = global_profile_scales
)

###########################################
# Save PDF with two pages
###########################################

log_message("Saving PDF with two pages")

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

log_message("Saving PNG summary")

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