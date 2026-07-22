source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("EnrichedHeatmap")
  library("cowplot")
  library("gtable")
  library("grid")
  library("ggplot2")
  library("gUtils")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################


samples <- read.csv(file = snakemake@input[["samplesheet"]])$sample
DER_clusters <- read.csv(snakemake@input[["sig_DERs"]], header = TRUE, row.names = 1)
bigwig_folder <- snakemake@params[["merged_bigwig_folder"]]

conditions <- unique(gsub("_R.*", "\\1", samples))

# Get the bigwig file names in the same order as the samples in the samplesheet
bigwig <- sapply(conditions,function(cond) {list.files(path = bigwig_folder, pattern = cond)})
bigwig_files <- paste0(bigwig_folder, "/", bigwig)

# Pick colors
set.seed(1234)
conditions_color <- randomcoloR::distinctColorPalette(length(unique(conditions)))
names(conditions_color) <- unique(conditions)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Get coverage information around the peaks from the bigwig files
get_coverage <- function(conditions, path_bw, peak_GR, dist){
  coverage_matrices <- lapply(path_bw, function(bw){
    reads <- rtracklayer::import(bw)
    coverage_matrix <- normalizeToMatrix(
      reads,
      peak_GR,
      background = 0,
      extend = dist, 
      w = 50,
      mean_mode = "coverage",
      value_column = "score",
      smooth = TRUE
    )
    return(coverage_matrix)
  })

  names(coverage_matrices) <- conditions
  return(coverage_matrices)
}

draw_heatmap <- function(matrices, peak_GR, conditions, max_coverage, dist){
  gradient <- c("#ffffff", "#98d344", "#629b0f")

  grad <- circlize::colorRamp2(c(0, round(max_coverage/2), max_coverage), rev(c("#138586", "#74cfb1", "#fffee8")))
  grad <- circlize::colorRamp2(c(0, round(max_coverage/2), max_coverage), gradient)

  axis_name <- c(
    paste0("-", dist/1000 ,"kb"), 
    "center", 
    paste0("+", dist/1000, "kb")
  )

  ht_list <- NULL
  ht_opt$TITLE_PADDING = unit(2, "mm")

  for (cond in conditions) {
    ht_list <- ht_list + EnrichedHeatmap(
        matrices[[cond]], 
        name = cond,
        column_title = cond,
        column_title_gp = gpar(
          fontface = "bold",
          fill = conditions_color[cond],
          col = "white",
          lwd = 0
        ),
        use_raster = TRUE, 
        axis_name = axis_name,
        raster_quality = 5,
        col = grad,
        axis_name_rot = 90,
        row_split = peak_GR$cluster,
        top_annotation = HeatmapAnnotation(
          lines = anno_enriched(
            ylim=c(0, max_coverage),
            gp = gpar(col = conditions_color, lex = 1.5),
            axis_param = list(side = "left", facing = "inside")
          )
        )
    )
  }
  heatmap <- draw(
    ht_list,
    ht_gap = unit(3, "mm")
  )
  return(heatmap)
}

###########################################
#                                         #
#              Draw heatmap               #
#                                         #
###########################################

distance_to_feature <- 10000

peaks <- GenomicRanges::GRanges(rownames(DER_clusters))
peaks$cluster <- DER_clusters

peaks <- gr.mid(peaks)

ATAC_mat <- get_coverage(conditions, bigwig_files, peaks, distance_to_feature)

# Get the max mean value of the covverage to set the legend scale comparable between condition
max_coverage <- max(unlist(lapply(ATAC_mat, colMeans)))
# Add 50% of the value to compensate for the aggregation of the signals in the heatmap
max_coverage <- ceiling(max_coverage *2)

ATAC_heatmap <- grid::grid.grabExpr(
  draw_heatmap(ATAC_mat, peaks, conditions, max_coverage, distance_to_feature)
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  ATAC_heatmap,
  base_width = 9,
  base_height = 15,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  ATAC_heatmap,
  base_width = 9,
  base_height = 15,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
