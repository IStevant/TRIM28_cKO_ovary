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

log_message("Starting differential region heatmap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("ComplexHeatmap")
  library("cowplot")
  library("grid")
  library("scales")
  library("randomcoloR")
  library("JLutils")
})

###########################################
# Load input data
###########################################

log_message("Loading normalised count matrix")

norm_counts <- read.csv(
  file = snakemake@input[["norm_counts"]],
  row.names = 1,
  check.names = FALSE
)

norm_counts <- as.matrix(norm_counts)
storage.mode(norm_counts) <- "numeric"

if (nrow(norm_counts) == 0) {
  log_error("Normalised count matrix has zero rows")
}

if (ncol(norm_counts) == 0) {
  log_error("Normalised count matrix has zero columns")
}

if (anyNA(norm_counts)) {
  log_error("Normalised count matrix contains NA values")
}

log_message("Number of regions: ", nrow(norm_counts))
log_message("Number of samples: ", ncol(norm_counts))

log_message("Loading significant differential regions")

load(snakemake@input[["sig_DERs"]])

if (!exists("filtered_DER_ids")) {
  log_error("Object 'filtered_DER_ids' was not found in sig_DERs file")
}

if (length(filtered_DER_ids) == 0) {
  log_error("No significant differential regions found")
}

log_message("Number of significant differential regions: ", length(filtered_DER_ids))

###########################################
# Prepare sample annotations
###########################################

log_message("Preparing sample annotations")

conditions <- gsub(
  "_R.*",
  "",
  colnames(norm_counts)
)

conditions <- factor(
  conditions,
  levels = unique(conditions)
)

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(levels(conditions))
)

names(conditions_color) <- levels(conditions)

###########################################
# Functions
###########################################

clip_matrix <- function(
  matrix,
  min_value,
  max_value
) {
  matrix[matrix < min_value] <- min_value
  matrix[matrix > max_value] <- max_value

  matrix
}

prepare_heatmap_matrix <- function(
  counts,
  selected_regions
) {
  log_message("Preparing heatmap matrix")

  selected_regions <- intersect(
    selected_regions,
    rownames(counts)
  )

  if (length(selected_regions) == 0) {
    log_error("None of the significant differential regions were found in the count matrix")
  }

  matrix_selected <- counts[selected_regions, , drop = FALSE]

  log_message("Number of regions retained for heatmap: ", nrow(matrix_selected))

  zscore_matrix <- t(
    scale(
      t(matrix_selected)
    )
  )

  zscore_matrix[is.na(zscore_matrix)] <- 0

  clip_matrix(
    zscore_matrix,
    min_value = -2,
    max_value = 2
  )
}

cluster_heatmap_rows <- function(
  matrix
) {
  log_message("Clustering heatmap rows")

  row_dend <- stats::hclust(
    stats::dist(matrix),
    method = "ward.D"
  )

  n_clusters <- JLutils::best.cutree(
    row_dend,
    min = 2,
    max = 15
  )

  log_message("Number of clusters selected: ", n_clusters)

  clusters <- stats::cutree(
    row_dend,
    k = n_clusters
  )

  clusters <- factor(
    clusters,
    levels = sort(unique(clusters))
  )

  names(clusters) <- rownames(matrix)

  list(
    row_dend = row_dend,
    clusters = clusters
  )
}

make_condition_annotation <- function(
  conditions,
  colours
) {
  ComplexHeatmap::HeatmapAnnotation(
    Conditions = ComplexHeatmap::anno_block(
      gp = grid::gpar(
        fill = colours,
        col = 0
      ),
      labels = names(colours),
      labels_gp = grid::gpar(
        col = "white",
        fontsize = 17,
        fontface = "bold"
      ),
      height = grid::unit(6.5, "mm")
    )
  )
}

make_cluster_annotation <- function() {
  ComplexHeatmap::rowAnnotation(
    Clusters = ComplexHeatmap::anno_empty(
      border = FALSE,
      width = grid::unit(10, "mm")
    )
  )
}

make_heatmap_colours <- function() {
  cold_palette <- grDevices::colorRampPalette(
    c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8")
  )

  warm_palette <- grDevices::colorRampPalette(
    c("#fffee8", "#f5da9f", "#ebb655", "#d18f43", "#b56832")
  )

  c(
    cold_palette(12),
    warm_palette(12)
  )
}

save_clusters <- function(
  clusters,
  heatmap_object,
  output_file
) {
  log_message("Writing cluster assignment file")

  split_order <- names(
    ComplexHeatmap::row_order(heatmap_object)
  )

  cluster_names <- letters[seq_along(split_order)]
  names(cluster_names) <- split_order

  cluster_labels <- cluster_names[as.character(clusters)]
  names(cluster_labels) <- names(clusters)

  write.csv(
    cluster_labels,
    file = output_file,
    quote = FALSE
  )

  cluster_labels
}

draw_cluster_labels <- function(
  cluster_names
) {
  for (i in seq_along(cluster_names)) {
    ComplexHeatmap::decorate_annotation(
      "Clusters",
      slice = i,
      {
        grid::grid.rect(
          x = 0.9,
          width = grid::unit(0.7, "mm"),
          gp = grid::gpar(
            fill = "black",
            col = NA
          ),
          just = "right"
        )

        grid::grid.circle(
          x = 0.3,
          r = grid::unit(3.8, "mm"),
          gp = grid::gpar(fill = "black")
        )

        grid::grid.text(
          x = 0.3,
          cluster_names[i],
          just = "center",
          gp = grid::gpar(
            fontsize = 17,
            col = "white"
          )
        )
      }
    )
  }
}

plot_simple_heatmap <- function(
  matrix,
  selected_regions,
  conditions,
  colours,
  cluster_file
) {
  heatmap_matrix <- prepare_heatmap_matrix(
    counts = matrix,
    selected_regions = selected_regions
  )

  clustering_results <- cluster_heatmap_rows(
    matrix = heatmap_matrix
  )

  conditions <- factor(
    conditions,
    levels = unique(conditions)
  )

  condition_annotation <- make_condition_annotation(
    conditions = conditions,
    colours = colours
  )

  cluster_annotation <- make_cluster_annotation()

  heatmap_colours <- make_heatmap_colours()

  heatmap <- ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = "z-score",
    top_annotation = condition_annotation,
    left_annotation = cluster_annotation,
    row_title_rot = 0,
    row_split = clustering_results$clusters,
    column_split = conditions,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = order(clustering_results$clusters),
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    cluster_column_slices = FALSE,
    col = heatmap_colours,
    heatmap_legend_param = list(direction = "vertical"),
    row_title = NULL,
    column_title = NULL,
    column_gap = grid::unit(0.4, "mm")
  )

  grid::grid.grabExpr(
    {
      drawn_heatmap <- ComplexHeatmap::draw(
        heatmap,
        column_title_gp = grid::gpar(
          fontsize = 12,
          fontface = "bold"
        ),
        row_title = paste(
          scales::comma(nrow(heatmap_matrix)),
          "differentially enriched regions"
        ),
        row_title_gp = grid::gpar(fontsize = 19),
        merge_legend = TRUE,
        use_raster = TRUE,
        raster_quality = 5
      )

      split_order <- names(
        ComplexHeatmap::row_order(drawn_heatmap)
      )

      cluster_names <- letters[seq_along(split_order)]
      names(cluster_names) <- split_order

      save_clusters(
        clusters = clustering_results$clusters,
        heatmap_object = drawn_heatmap,
        output_file = cluster_file
      )

      draw_cluster_labels(cluster_names)
    }
  )
}

###########################################
# Draw heatmap
###########################################

log_message("Drawing heatmap")

dynamic_peaks <- plot_simple_heatmap(
  matrix = norm_counts,
  selected_regions = filtered_DER_ids,
  conditions = conditions,
  colours = conditions_color,
  cluster_file = snakemake@output[["clusters"]]
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")