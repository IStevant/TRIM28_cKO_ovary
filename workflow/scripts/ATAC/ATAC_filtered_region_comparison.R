source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("cowplot")
  library("grid")
  library("viridis")
  library("ggplot2")
  library("ComplexHeatmap")
  library("ChIPpeakAnno")
  library("futile.logger")
  library("eulerr")
  library("dplyr")
  library("shades")
  library("ggrepel")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

counts <- read.csv(file = snakemake@input[["counts"]], header=TRUE, row.names = 1)
norm_data <- read.csv(file = snakemake@input[["norm_data"]], row.names = 1)
raw_peak_folder <- snakemake@params[["raw_peak_folder"]]
min_peak_rep <- snakemake@params[["min_peak_rep"]]
minReads <- snakemake@params[["minReads"]]
promoter <- snakemake@params[["promoter"]]
conditions <- gsub("_REP.*", "\\1", colnames(norm_data))

# Pick colors
set.seed(1234)
conditions_color <- randomcoloR::distinctColorPalette(length(unique(conditions)))

names(conditions_color) <- unique(conditions)
# conditions_color <- conditions_color[order(names(conditions_color))]

corr_method <- snakemake@params[["corr_method"]]

###########################################
#                                         #
#          ChIPseeker options             #
#                                         #
###########################################

# Ignore unnecessary annotation
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

#################################################################################################################################

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Draw the correlation matrix between the samples
#' @param matrix Expression matrix.
#' @param conditions Vector containing the names of the conditions of the samples. The length should be the same as the number of samples.
#' @param method Method for the correlation (Example: "Spearman" or "Pearson"). Default is "Spearman".
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @return Pheatmap object.
correlation <- function(matrix, method = "Spearman", colours) {
  # matrix <- matrix[, order(names(matrix))]
  cor_data <- cor(matrix, method = method)
  cor_data[cor_data == 1.000] <- NA
  conditions <- gsub("_REP.*", "\\1", colnames(matrix))

  anno <- data.frame(
    Samples = conditions
  )

  colors <- list(Samples=conditions_color)

  heatmap <- Heatmap(
    cor_data,
    name = "Correlation",
    cluster_rows = F, 
    cluster_columns = F, 
    col = viridis::viridis(n = 100,option = 'C'),
    left_annotation = rowAnnotation(df = anno, col=colors, annotation_legend_param = list(at = names(conditions_color))),
    top_annotation = columnAnnotation(df = anno, col=colors, show_legend=c(F,F,F)),
    width = unit(0.65, "snpc"),
    height = unit(0.65, "snpc")
  )

  heatmap <- draw(
      heatmap,
      column_title=paste0("Pairwise sample correlation (",method,")"),
      column_title_gp=grid::gpar(fontsize=12, fontface="bold")
  )

  plot <- grid.grabExpr(
    draw(heatmap)
  )

  return(plot)
}

#' Proceed to the PCA using prcomp.
#' @param matrix Expression matrix.
#' @return prcomp object.
run_pca <- function(matrix) {
  print("Calculating the PCA...")
  t.matrix <- t(matrix)
  t.matrix.no0 <- t.matrix[, colSums(t.matrix) != 0]
  pca <- prcomp(
    t.matrix.no0,
    center = TRUE,
    scale. = TRUE
  )
  return(pca)
}

#' Plot the PCA.
#' @param matrix Expression matrix.
#' @param conditions Vector containing the names of the conditions of the samples. The length should be the same as the number of samples.
#' @param colours Vector containing the hexadecimal colours corresponding to each conditions.
#' @param PCs Vector PCs to plot. Default is "PC1" vs "PC2". If more than two PCs provided, all the possible combinations are drawn.
#' @return Pheatmap object.
plot.pca <- function(matrix, colours, PCs = c("PC1", "PC2")) {
  pca <- run_pca(matrix)
  print("Generating the plots...")
  percent_var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  cond <- factor(conditions)
  col <- factor(conditions)
  levels(col) <- colours
  col <- as.vector(col)
  conditions <- gsub("_REP.*", "\\1", colnames(matrix))
  replicates <- stringr::str_to_sentence(sub(".*_(REP[0-9]+)", "\\1", colnames(matrix)))
  scores <- as.data.frame(pca$x)
  PCs.combinations <- combn(PCs, 2)
  plots <- apply(
    PCs.combinations,
    2,
    function(combination) {
      data <- scores[, c(combination[1], combination[2])]
      data$cond <- conditions
      data$rep <- replicates
      colnames(data) <- c("PC_x", "PC_y", "cond")
      plot <- ggplot(data, aes(x = PC_x, y = PC_y, fill = cond, color = cond)) +
        geom_label_repel(
          aes(
            label = replicates
            # color = cond
          ),
          box.padding = 0.7,
          color = "#444444",
          size = 4,
          fill = alpha("white", 0.5),
          force = 1,
          max.overlaps = 25
        ) +
        geom_point(shape = 21, size = 6, stroke = 0.5, color = "#333333") +
        scale_fill_manual(values = conditions_color, breaks=names(conditions_color)) +
        ggtitle("PCA on open chromatin region accessibility")+
        theme_bw() +
        xlab(paste(combination[1], " ", "(", round(percent_var_explained[as.numeric(gsub("PC", "", combination[1]))], digit = 2), "%)", sep = "")) +
        ylab(paste(combination[2], " ", "(", round(percent_var_explained[as.numeric(gsub("PC", "", combination[2]))], digit = 2), "%)", sep = "")) +
        theme(
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          aspect.ratio = 1,
          legend.title = element_blank(),
          legend.text = element_text(size = 12)
        )
      return(plot)
    }
  )
  print("Done.")
  return(plots)
}


#' Plot the number of peaks per condition
#' @param counts list of peaks from each condition.
#' @return Return a ggplot object.
plot_venn <- function(counts, colours) {

  conditions <- unique(gsub("_REP.*", "\\1", colnames(counts)))

  split_by_conditions <- function(df, conditions) {
    lapply(conditions, function(condition) {
      df[, grep(condition, colnames(df)), drop = FALSE]
    }) |> setNames(conditions)
  }

  split_conds <- split_by_conditions(counts, conditions)

  print("filter replicated peaks")
  regions <- lapply(split_conds, function(cond){
      is_open <- ifelse(cond >= minReads,1,0)
      is_open <- is_open[rowSums(is_open) >= 1, , drop = FALSE]
      open_regions <- rownames(is_open)
      return(open_regions)
    }) |> setNames(conditions)

  print("Overlap OCR per condition")

  venn_data <- euler(regions)

  # Get the sample order to apply the correct colors
  sample_order <- grep("&", names(venn_data$original.values), invert = TRUE, value= TRUE)
  colours_order <- colours[sample_order]

  # Get the total number of open chromatin region per condition
  # set_names <- names(venn)
  totals <- unlist(lapply(sample_order, function(set) {
    total <- sum(venn_data$original.values[grep(set, names(venn_data$original.values))])
    return(total)
  }))
  names(totals) <- sample_order


  venn_plot <- plot(
    venn_data,
    labels = list(
      labels = paste0(sample_order, "\n(", totals, ")"),
      font = 2,
      cex = 1.2
    ),
    quantities = list(
      fontsize = 12
    ),
    edges = list(
      col = colours_order, 
      lex = 2
    ),
    fills = list(
      fill= colours_order,
      alpha=0.45
    )
  )

  venn_plot$vp$width <- unit(0.8, "npc")
  venn_plot$vp$height <- unit(0.8, "npc")

  title <- ggdraw() + draw_label("Overlap of the open chromatin regions\nbetween conditions", fontface='bold', size=12)

  venn_plot <- plot_grid(
    title,
    venn_plot, 
    ncol=1,
    nrow=2,
    rel_heights = c(0.1, 1)
  )

  return(venn_plot)
}

#' Generate the read count matrix
#' @param anno annotation list containing outputs from ChIPseeker::annotatePeak.
#' @return Return a dataframe.
plot_anno <- function(counts, colours) {

  # Load mouse genome

  gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
  alt_gtf_url <- "http://ftp.cbi.pku.edu.cn/pub/mirror/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

  Genes <- tryCatch({
    rtracklayer::import(gtf_url)
   }, error = function(e) {
    rtracklayer::import(alt_gtf_url)
   })

  txdb <- txdbmaker::makeTxDbFromGRanges(
    Genes,
    drop.stop.codons = FALSE
  )

  conditions <- unique(gsub("_REP.*", "\\1", colnames(counts)))

  split_by_conditions <- function(df, conditions) {
    lapply(conditions, function(condition) {
      df[, grep(condition, colnames(df)), drop = FALSE]
    }) |> setNames(conditions)
  }

  split_conds <- split_by_conditions(counts, conditions)

  regions <- lapply(split_conds, function(cond){
      is_open <- ifelse(cond > minReads,1,0)
      is_open <- is_open[rowSums(is_open) >= min_peak_rep, , drop = FALSE]
      open_gr <- GRanges(rownames(is_open))
      open_gr$names <- rownames(is_open)
      return(open_gr)
    }) |> setNames(conditions)

  anno <- lapply(
    regions,
    function(gr) {
      ChIPseeker::annotatePeak(
        gr,
        tssRegion = c(-promoter, 0),
        TxDb = txdb,
        overlap = "all"
      )
    }
  ) |> setNames(conditions)


  anno <- lapply(1:length(anno), function(x) {
    stat <- anno[[x]]@annoStat
    stat$Conditions <- rep(conditions[x], nrow(stat))
    return(stat)
  })

  data <- data.table::rbindlist(anno)

  labels <- round(data$Frequency, digit = 2)
  labels[labels < 4] <- 0
  labels <- paste0(labels, "%")
  labels[labels == "0%"] <- " "

  colors <- c(
    "#f8777c",
    "#0e1b47",
    "#4461a8",
    "#21a3ea",
    "#3bc9d6",
    "#b886da"
  )

  plot <- ggplot(data, aes(fill = Feature, x = factor(Conditions, level = unique(conditions)), y = Frequency)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = labels, color = Feature), size = 5, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = alpha(colors, 0.8)) +
    scale_color_manual(values = c("black", "white", "white", "black", "black", "black"), guide = "none") +
    xlab("Conditions") +
    ggtitle("Open chromatin region\ngenomic location") +
    # facet_wrap(~Sex, nrow = 1) +
    theme_light() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, margin = margin(r = 10, unit = "pt")),
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.box.spacing = unit(0, "mm"),
    )
  return(plot)
}

#################################################################################################################################

###########################################
#                                         #
#        Plot correlation and PCA         #
#                                         #
###########################################

corr_plot <- correlation(norm_data, corr_method, conditions_color)
pca_plot <- plot.pca(norm_data, conditions_color, c("PC1", "PC2"))

venn_plot <- plot_venn(counts, conditions_color)
anno_plot <- plot_anno(counts, conditions_color)

# Combine the two plots
figure <- plot_grid(
  plotlist = list(corr_plot, pca_plot[[1]], venn_plot, anno_plot),
  # plotlist = list(corr_plot, pca_plot[[1]]),
  labels = "AUTO",
  ncol = 2, 
  align = "h"
)

##########################################
#                                        #
#               Save plots               #
#                                        #
##########################################

# As PDF
save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 30,
  base_height = 20,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 30,
  base_height = 20,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
