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

log_message("Starting DAR heatmap script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
  library("grid")
  library("ggplot2")
  library("randomcoloR")
  library("scales")
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

log_message("Loading significant DARs")

load(
  snakemake@input[["sig_DARs"]]
)

if (!exists("filtered_DARs")) {
  log_error("Object 'filtered_DARs' was not found in sig_DARs file")
}

if (length(filtered_DARs) == 0) {
  log_error("No significant DARs found")
}

conditions <- gsub(
  "_R.*",
  "",
  colnames(norm_counts)
)

###########################################
# Prepare colours
###########################################

log_message("Preparing condition colours")

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(unique(conditions))
)

names(conditions_color) <- unique(conditions)

###########################################
# Functions
###########################################

make_condition_annotation <- function(
  conditions,
  conditions_color
) {
  ComplexHeatmap::HeatmapAnnotation(
    Conditions = ComplexHeatmap::anno_block(
      gp = grid::gpar(
        fill = conditions_color,
        col = NA
      ),
      labels = names(conditions_color),
      labels_gp = grid::gpar(
        col = "white",
        fontsize = 15,
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
    ),
    show_annotation_name = FALSE
  )
}

write_cluster_table <- function(
  clustering,
  row_order_list,
  output_file,
  heatmap_matrix
) {
  split_order <- names(row_order_list)

  wt_cols <- grep(
    "WT",
    colnames(heatmap_matrix),
    ignore.case = TRUE
  )

  ko_cols <- grep(
    "KO",
    colnames(heatmap_matrix),
    ignore.case = TRUE
  )

  if (length(wt_cols) == 0 || length(ko_cols) == 0) {
    log_message(
      "Could not identify WT/KO columns. Keeping heatmap split order for cluster naming."
    )

    ordered_cluster_ids <- split_order
  } else {
    cluster_score <- sapply(split_order, function(cluster_id) {
      rows <- names(clustering)[
        as.character(clustering) == cluster_id
      ]

      mean(
        heatmap_matrix[rows, wt_cols, drop = FALSE],
        na.rm = TRUE
      ) -
        mean(
          heatmap_matrix[rows, ko_cols, drop = FALSE],
          na.rm = TRUE
        )
    })

    ordered_cluster_ids <- names(
      sort(
        cluster_score,
        decreasing = TRUE
      )
    )
  }

  cluster_names <- letters[
    seq_along(ordered_cluster_ids)
  ]

  names(cluster_names) <- ordered_cluster_ids

  ordered_regions <- unlist(
    lapply(
      row_order_list[ordered_cluster_ids],
      function(idx) {
        names(clustering)[idx]
      }
    ),
    use.names = FALSE
  )

  output_clusters <- cluster_names[
    as.character(clustering[ordered_regions])
  ]

  output_table <- data.frame(
    region = ordered_regions,
    cluster = as.character(output_clusters),
    stringsAsFactors = FALSE
  )

  write.csv(
    output_table,
    file = output_file,
    quote = FALSE,
    row.names = FALSE
  )

  cluster_names
}

decorate_cluster_annotation <- function(
  cluster_names
) {
  for (i in seq_along(cluster_names)) {
    ComplexHeatmap::decorate_annotation(
      "Clusters",
      slice = i,
      {
        grid::grid.rect(
          x = grid::unit(9.5, "mm"),
          width = grid::unit(0.7, "mm"),
          just = "right",
          gp = grid::gpar(
            fill = "black",
            col = NA
          )
        )

        grid::grid.circle(
          x = grid::unit(3, "mm"),
          r = grid::unit(3.8, "mm"),
          gp = grid::gpar(
            fill = "black",
            col = NA
          )
        )

        grid::grid.text(
          cluster_names[i],
          x = grid::unit(3, "mm"),
          just = "center",
          gp = grid::gpar(
            fontsize = 17,
            col = "white",
            fontface = "bold"
          )
        )
      }
    )
  }
}

plot_simple_heatmap <- function(
  data,
  de_features,
  colours,
  cluster_file
) {
  log_message("Preparing DAR matrix")

  matrix_DAR <- data[
    rownames(data) %in% de_features,
    ,
    drop = FALSE
  ]

  if (nrow(matrix_DAR) == 0) {
    log_error("No DARs from sig_DARs were found in norm_counts")
  }

  log_message("Number of DARs in heatmap: ", nrow(matrix_DAR))

  log_message("Computing z-scores")

  heatmap_matrix <- t(
    scale(
      t(matrix_DAR)
    )
  )

  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  set.seed(654)

  log_message("Clustering rows")

  row_dend <- stats::hclust(
    stats::dist(heatmap_matrix),
    method = "ward.D"
  )

  n_clusters <- JLutils::best.cutree(
    row_dend,
    min = 2,
    max = 15
  )

  log_message("Number of clusters selected: ", n_clusters)

  clustering <- stats::cutree(
    row_dend,
    k = n_clusters
  )

  names(clustering) <- rownames(heatmap_matrix)

  clustering <- factor(
    clustering,
    levels = sort(unique(clustering))
  )

  heatmap_matrix[heatmap_matrix > 2] <- 2
  heatmap_matrix[heatmap_matrix < -2] <- -2

  sample_conditions <- gsub(
    "_R.*",
    "",
    colnames(heatmap_matrix)
  )

  condition_annotation <- make_condition_annotation(
    conditions = sample_conditions,
    conditions_color = colours
  )

  cluster_annotation <- make_cluster_annotation()

  cold <- grDevices::colorRampPalette(
    c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8")
  )

  warm <- grDevices::colorRampPalette(
    c("#fffee8", "#f5da9f", "#ebb655", "#d18f43", "#b56832")
  )

  colour_function <- circlize::colorRamp2(
    breaks = seq(-2, 2, length.out = 24),
    colors = c(cold(12), warm(12))
  )

  heatmap <- ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = "z-score",
    top_annotation = condition_annotation,
    left_annotation = cluster_annotation,
    row_title_rot = 0,
    row_split = clustering,
    column_split = factor(
      sample_conditions,
      levels = unique(sample_conditions)
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    col = colour_function,
    heatmap_legend_param = list(
      direction = "vertical"
    ),
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
          "differentially accessible regions"
        ),
        row_title_gp = grid::gpar(
          fontsize = 19
        ),
        merge_legend = TRUE,
        use_raster = TRUE,
        raster_quality = 5
      )

      cluster_names <- write_cluster_table(
        clustering = clustering,
        row_order_list = ComplexHeatmap::row_order(drawn_heatmap),
        output_file = cluster_file,
        heatmap_matrix = heatmap_matrix
      )

      decorate_cluster_annotation(
        cluster_names = cluster_names
      )
    },
    height = 9,
    width = 5
  )
}

###########################################
# Draw heatmap
###########################################

dynamic_peaks <- plot_simple_heatmap(
  data = norm_counts,
  de_features = filtered_DARs,
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