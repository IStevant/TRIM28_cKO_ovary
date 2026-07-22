source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("ComplexHeatmap")
  library("cowplot")
  library("gtable")
  library("grid")
  library("ggplot2")
  library("clusterProfiler")
  library("JLutils")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

norm_counts <- read.csv(file = snakemake@input[["norm_counts"]], row.names = 1)
load(snakemake@input[["sig_DERs"]])
conditions <- gsub("_R.*", "\\1", colnames(norm_counts))

# Pick colors
set.seed(1234)
conditions_color <- randomcoloR::distinctColorPalette(length(unique(conditions)))
names(conditions_color) <- unique(conditions)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the heatmap of the z-scores from the regions differentially accessible along developmental conds
#' @param data Full expression matrix (TPM or normalized counts).
#' @param de_feature Vector containing the names of the regions found as differentially accessible.
#' @param colors Vector containing the hexadecimal colours corresponding to each conditions.
#' @param res_file Name of the output file containing the list of regions per cluster.
#' @return Pheatmap object.
plot_simple_heatmap <- function(data, de_feature, colors, res_file) {
  matrix_DER <- data[rownames(data) %in% de_feature, ]


  # Calculate z-scores
  matrix <- t(scale(t(matrix_DER)))

  set.seed(654)

  # Cluster matrix
  row_dend <- hclust(dist(matrix), method = "ward.D")

  # Calculate the optimal number of clusters
  clusters <- best.cutree(row_dend, min = 2, max = 15)

  clustering <- cutree(row_dend, k = clusters)

  clustering <- as.factor(cutree(row_dend, k = clusters))
  
  # Limit zscore to |2|
  matrix[matrix > 2] <- 2
  matrix[matrix < (-2)] <- (-2)

  # Prepare top annotation
  conditions <- gsub("_R.*", "\\1", colnames(matrix))

  # annotation_col <- data.frame(
  #   Conditions = conditions
  # )

  # colnames(annotation_col) <- conditions

  # Conditions palette
  col_cond <- colors


  # Color palette for the heatmap
  # blue-yellow-red
  cold <- colorRampPalette(c("#4677b7", "#709eca", "#9ac4dd", "#cce1e3", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#fdd4ab", "#fbab70", "#e96e33", "#d83329"))
  BYR <- c(cold(12), warm(12))
  # green-yellow-brown
  cold <- colorRampPalette(c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#f5da9f", "#ebb655", "#d18f43", "#b56832"))
  GYB <- c(cold(12), warm(12))
  # green-yellow-purple
  cold <- colorRampPalette(c("#138586", "#44aa9b", "#74cfb1", "#bbe7ce", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#f3c2d3", "#d981cf", "#b452c1", "#9030b4"))
  GYP <- c(cold(12), warm(12))

  mypalette <- GYB

  # Top annotation (conds)
  cond_anno <- HeatmapAnnotation(
    Conditions = anno_block(
      gp = gpar(fill = col_cond, col = 0),
      labels = names(col_cond),
      # at = names(col_cond),
      labels_gp = gpar(col = "white", fontsize = 17, fontface = "bold"),
      height = unit(6.5, "mm")
    )
  )

  # Row annotation (peak clusters)
  cluster_anno <- rowAnnotation(
    Clusters = anno_empty(
      border = FALSE,
      width = unit(10, "mm")
    )
  )


  # Make the heatmap
  ht_list <- Heatmap(
    matrix,
    name = "z-score",
    top_annotation = cond_anno,
    left_annotation = cluster_anno,
    row_title_rot = 0,
    row_split = clustering,
    column_split = factor(conditions, levels=unique(conditions)),
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    cluster_column_slices = FALSE,
    col = mypalette,
    heatmap_legend_param = list(direction = "vertical"),
    row_title = NULL,
    column_title = NULL,
    column_gap = unit(0.4, "mm")
  )

  height <- 9
  width <- 5

  gTree <- grid.grabExpr(
    {
      # Draw the heatmap
      ht <- draw(
        ht_list,
        column_title_gp=grid::gpar(fontsize=12, fontface="bold"),
        row_title = paste(scales::comma(nrow(matrix)), "differentially enriched regions"),
        row_title_gp = gpar(fontsize = 19),
        # annotation_legend_list=cluster_legend,
        merge_legend = TRUE,
        use_raster = TRUE,
        raster_quality = 5
      )

      # Add peak cluster annotation as extra-rectangles to avoid bad rasterization (faded colors)
      split_order <- names(row_order(ht))
      cluster_names <- letters[seq_along(split_order)]
      # levels(split_order) <- cluster_names
      names(cluster_names) <- split_order
      peaks <- names(clustering)
      clustering <- cluster_names[as.character(clustering)]
      names(clustering) <- peaks
      # Save clustering
      write.csv(clustering, file = res_file, quote = FALSE)

      for (i in seq_along(cluster_names)) {
        decorate_annotation(
          "Clusters",
          slice = i,
          {
            grid.rect(x = 0.9, width = unit(0.7, "mm"), gp = gpar(fill = "black", col = NA), just = "right")
            grid.circle(x = 0.3, r = unit(3.8, "mm"), gp = gpar(fill = "black"))
            grid.text(x = 0.3, cluster_names[i], just = "center", gp = gpar(fontsize = 17, col = "white"))
          }
        )
      }
    },
    height = height,
    width = width
  )

  return(gTree)
}

###########################################
#                                         #
#              Draw heatmap               #
#                                         #
###########################################

matrix <- norm_counts

dynamic_peaks <- plot_simple_heatmap(
  matrix,
  filtered_DERs,
  conditions_color,
  snakemake@output[["clusters"]]
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  dynamic_peaks,
  base_width = 14,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
