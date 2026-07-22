source(".Rprofile")

suppressPackageStartupMessages({
  library("ComplexHeatmap")
  library("cowplot")
  library("eulerr")
  library("ggplot2")
  library("ggrepel")
  library("grid")
  library("viridis")
})


###########################################
# Logging functions
###########################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}


###########################################
# Snakemake inputs, outputs and parameters
###########################################

log_message("Reading Snakemake inputs, outputs and parameters.")

counts_file <- snakemake@input[["counts"]]
norm_data_file <- snakemake@input[["norm_data"]]

raw_peak_folder <- snakemake@params[["raw_peak_folder"]]
min_peak_rep <- snakemake@params[["min_peak_rep"]]
minReads <- snakemake@params[["minReads"]]
promoter <- snakemake@params[["promoter"]]
corr_method <- snakemake@params[["corr_method"]]

output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]


###########################################
# Load data
###########################################

log_message("Loading raw count matrix: ", counts_file)

if (!file.exists(counts_file)) {
  log_error("Raw count matrix does not exist: ", counts_file)
}

counts <- read.csv(
  file = counts_file,
  header = TRUE,
  row.names = 1
)

log_message("Loading normalised count matrix: ", norm_data_file)

if (!file.exists(norm_data_file)) {
  log_error("Normalised count matrix does not exist: ", norm_data_file)
}

norm_data <- read.csv(
  file = norm_data_file,
  header = TRUE,
  row.names = 1
)

if (!identical(colnames(counts), colnames(norm_data))) {
  log_error(
    "Raw and normalised count matrices do not contain samples in the same order."
  )
}

conditions <- gsub(
  "_REP.*",
  "\\1",
  colnames(norm_data)
)

log_message(
  "Detected conditions: ",
  paste(unique(conditions), collapse = ", ")
)


###########################################
# Plot parameters
###########################################

set.seed(1234)

conditions_color <- randomcoloR::distinctColorPalette(
  length(unique(conditions))
)

names(conditions_color) <- unique(conditions)


###########################################
# ChIPseeker options
###########################################

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)


###########################################
# Functions
###########################################

#' Draw the correlation matrix between samples
#'
#' @param matrix Expression matrix.
#' @param method Correlation method.
#' @param colours Named vector containing one colour per condition.
#'
#' @return A graphical object containing the correlation heatmap.
correlation <- function(matrix, method = "spearman", colours) {

  log_message(
    "Calculating pairwise sample correlations using ",
    method,
    "."
  )

  method <- tolower(method)

  cor_data <- cor(
    matrix,
    method = method
  )

  cor_data[cor_data == 1] <- NA

  sample_conditions <- gsub(
    "_REP.*",
    "\\1",
    colnames(matrix)
  )

  annotation <- data.frame(
    Samples = sample_conditions,
    row.names = colnames(matrix)
  )

  annotation_colours <- list(
    Samples = colours
  )

  heatmap <- ComplexHeatmap::Heatmap(
    cor_data,
    name = "Correlation",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = viridis::viridis(
      n = 100,
      option = "C"
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      df = annotation,
      col = annotation_colours,
      annotation_legend_param = list(
        at = names(colours)
      )
    ),
    top_annotation = ComplexHeatmap::columnAnnotation(
      df = annotation,
      col = annotation_colours,
      show_legend = FALSE
    ),
    width = grid::unit(
      0.65,
      "snpc"
    ),
    height = grid::unit(
      0.65,
      "snpc"
    )
  )

  correlation_plot <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      column_title = paste0(
        "Pairwise sample correlation (",
        method,
        ")"
      ),
      column_title_gp = grid::gpar(
        fontsize = 12,
        fontface = "bold"
      )
    )
  )

  return(correlation_plot)
}


#' Perform principal component analysis
#'
#' @param matrix Expression matrix.
#'
#' @return A prcomp object.
run_pca <- function(matrix) {

  log_message("Calculating principal component analysis.")

  transposed_matrix <- t(matrix)

  transposed_matrix <- transposed_matrix[
    ,
    colSums(transposed_matrix) != 0,
    drop = FALSE
  ]

  if (ncol(transposed_matrix) == 0) {
    log_error(
      "PCA cannot be performed because all regions have zero counts."
    )
  }

  pca <- prcomp(
    transposed_matrix,
    center = TRUE,
    scale. = TRUE
  )

  return(pca)
}


#' Plot principal component analysis
#'
#' @param matrix Expression matrix.
#' @param colours Named vector containing one colour per condition.
#' @param PCs Principal components to plot.
#'
#' @return A list of ggplot objects.
plot_pca <- function(
    matrix,
    colours,
    PCs = c("PC1", "PC2")) {

  pca <- run_pca(matrix)

  log_message("Generating PCA plots.")

  percent_var_explained <- (
    pca$sdev^2 / sum(pca$sdev^2)
  ) * 100

  sample_conditions <- gsub(
    "_REP.*",
    "\\1",
    colnames(matrix)
  )

  replicates <- stringr::str_to_sentence(
    sub(
      ".*_(REP[0-9]+)",
      "\\1",
      colnames(matrix)
    )
  )

  scores <- as.data.frame(pca$x)

  missing_PCs <- setdiff(
    PCs,
    colnames(scores)
  )

  if (length(missing_PCs) > 0) {
    log_error(
      "The following principal components are unavailable: ",
      paste(missing_PCs, collapse = ", ")
    )
  }

  PC_combinations <- combn(
    PCs,
    2
  )

  plots <- apply(
    PC_combinations,
    2,
    function(combination) {

      data <- data.frame(
        PC_x = scores[[combination[1]]],
        PC_y = scores[[combination[2]]],
        condition = sample_conditions,
        replicate = replicates
      )

      ggplot(
        data,
        aes(
          x = PC_x,
          y = PC_y,
          fill = condition
        )
      ) +
        geom_label_repel(
          aes(label = replicate),
          box.padding = 0.7,
          color = "#444444",
          size = 4,
          fill = alpha("white", 0.5),
          force = 1,
          max.overlaps = 25
        ) +
        geom_point(
          shape = 21,
          size = 6,
          stroke = 0.5,
          color = "#333333"
        ) +
        scale_fill_manual(
          values = colours,
          breaks = names(colours)
        ) +
        xlab(
          paste0(
            combination[1],
            " (",
            round(
              percent_var_explained[
                as.numeric(
                  gsub(
                    "PC",
                    "",
                    combination[1]
                  )
                )
              ],
              digits = 2
            ),
            "%)"
          )
        ) +
        ylab(
          paste0(
            combination[2],
            " (",
            round(
              percent_var_explained[
                as.numeric(
                  gsub(
                    "PC",
                    "",
                    combination[2]
                  )
                )
              ],
              digits = 2
            ),
            "%)"
          )
        ) +
        ggtitle(
          "PCA on open chromatin region accessibility"
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(
            size = 12,
            face = "bold",
            hjust = 0.5
          ),
          axis.text = element_text(
            size = 12
          ),
          axis.title = element_text(
            size = 12
          ),
          aspect.ratio = 1,
          legend.title = element_blank(),
          legend.text = element_text(
            size = 12
          )
        )
    }
  )

  log_message("PCA plot generation completed.")

  return(plots)
}


#' Plot the overlap of open chromatin regions between conditions
#'
#' @param counts Raw count matrix.
#' @param colours Named vector containing one colour per condition.
#'
#' @return A graphical object containing an Euler diagram.
plot_venn <- function(counts, colours) {

  log_message("Identifying open chromatin regions per condition.")

  conditions <- unique(
    gsub(
      "_REP.*",
      "\\1",
      colnames(counts)
    )
  )

  split_by_conditions <- function(data, conditions) {

    condition_data <- lapply(
      conditions,
      function(condition) {
        data[
          ,
          grep(
            condition,
            colnames(data)
          ),
          drop = FALSE
        ]
      }
    )

    names(condition_data) <- conditions

    return(condition_data)
  }

  split_counts <- split_by_conditions(
    counts,
    conditions
  )

  regions <- lapply(
    split_counts,
    function(condition_counts) {

      is_open <- ifelse(
        condition_counts >= minReads,
        1,
        0
      )

      is_open <- is_open[
        rowSums(is_open) >= 1,
        ,
        drop = FALSE
      ]

      rownames(is_open)
    }
  )

  names(regions) <- conditions

  log_message("Calculating overlap between conditions.")

  venn_data <- eulerr::euler(regions)

  sample_order <- grep(
    "&",
    names(venn_data$original.values),
    invert = TRUE,
    value = TRUE
  )

  colours_order <- colours[sample_order]

  totals <- unlist(
    lapply(
      sample_order,
      function(condition) {
        sum(
          venn_data$original.values[
            grep(
              condition,
              names(venn_data$original.values)
            )
          ]
        )
      }
    )
  )

  names(totals) <- sample_order

  venn_plot <- plot(
    venn_data,
    labels = list(
      labels = paste0(
        sample_order,
        "\n(",
        totals,
        ")"
      ),
      font = 2,
      cex = 1.2
    ),
    quantities = list(
      fontsize = 12
    ),
    edges = list(
      col = colours_order,
      lwd = 2
    ),
    fills = list(
      fill = colours_order,
      alpha = 0.45
    )
  )

  venn_plot$vp$width <- grid::unit(
    0.8,
    "npc"
  )

  venn_plot$vp$height <- grid::unit(
    0.8,
    "npc"
  )

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Overlap of the open chromatin regions\nbetween conditions",
      fontface = "bold",
      size = 12
    )

  venn_plot <- cowplot::plot_grid(
    title,
    venn_plot,
    ncol = 1,
    nrow = 2,
    rel_heights = c(0.1, 1)
  )

  return(venn_plot)
}


#' Plot the genomic annotation of open chromatin regions
#'
#' @param counts Raw count matrix.
#' @param colours Named vector containing one colour per condition.
#'
#' @return A ggplot object.
plot_anno <- function(counts, colours) {

  log_message("Loading GENCODE mouse genome annotation.")

  gtf_url <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/",
    "Gencode_mouse/release_M25/",
    "gencode.vM25.annotation.gtf.gz"
  )

  alt_gtf_url <- paste0(
    "http://ftp.cbi.pku.edu.cn/pub/mirror/GENCODE/",
    "Gencode_mouse/release_M25/",
    "gencode.vM25.annotation.gtf.gz"
  )

  genes <- tryCatch(
    {
      rtracklayer::import(gtf_url)
    },
    error = function(e) {
      log_message(
        "Primary GENCODE mirror unavailable. Trying alternative mirror."
      )

      rtracklayer::import(alt_gtf_url)
    }
  )

  txdb <- txdbmaker::makeTxDbFromGRanges(
    genes,
    drop.stop.codons = FALSE
  )

  conditions <- unique(
    gsub(
      "_REP.*",
      "\\1",
      colnames(counts)
    )
  )

  split_by_conditions <- function(data, conditions) {

    condition_data <- lapply(
      conditions,
      function(condition) {
        data[
          ,
          grep(
            condition,
            colnames(data)
          ),
          drop = FALSE
        ]
      }
    )

    names(condition_data) <- conditions

    return(condition_data)
  }

  split_counts <- split_by_conditions(
    counts,
    conditions
  )

  log_message("Selecting replicated open chromatin regions.")

  regions <- lapply(
    split_counts,
    function(condition_counts) {

      is_open <- ifelse(
        condition_counts > minReads,
        1,
        0
      )

      is_open <- is_open[
        rowSums(is_open) >= min_peak_rep,
        ,
        drop = FALSE
      ]

      open_regions <- GenomicRanges::GRanges(
        rownames(is_open)
      )

      open_regions$names <- rownames(is_open)

      return(open_regions)
    }
  )

  names(regions) <- conditions

  log_message("Annotating open chromatin regions.")

  annotations <- lapply(
    regions,
    function(regions_gr) {
      ChIPseeker::annotatePeak(
        regions_gr,
        tssRegion = c(-promoter, 0),
        TxDb = txdb,
        overlap = "all"
      )
    }
  )

  names(annotations) <- conditions

  annotation_stats <- lapply(
    seq_along(annotations),
    function(index) {

      stats <- annotations[[index]]@annoStat

      stats$Conditions <- rep(
        conditions[index],
        nrow(stats)
      )

      return(stats)
    }
  )

  data <- data.table::rbindlist(
    annotation_stats
  )

  labels <- round(
    data$Frequency,
    digits = 2
  )

  labels[labels < 4] <- 0
  labels <- paste0(labels, "%")
  labels[labels == "0%"] <- " "

  annotation_colours <- c(
    "#f8777c",
    "#0e1b47",
    "#4461a8",
    "#21a3ea",
    "#3bc9d6",
    "#b886da"
  )

  annotation_plot <- ggplot(
    data,
    aes(
      fill = Feature,
      x = factor(
        Conditions,
        levels = unique(conditions)
      ),
      y = Frequency
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
    geom_text(
      aes(
        label = labels,
        color = Feature
      ),
      size = 5,
      position = position_stack(
        vjust = 0.5
      )
    ) +
    scale_y_continuous(
      labels = scales::comma
    ) +
    scale_fill_manual(
      values = alpha(
        annotation_colours,
        0.8
      )
    ) +
    scale_color_manual(
      values = c(
        "black",
        "white",
        "white",
        "black",
        "black",
        "black"
      ),
      guide = "none"
    ) +
    xlab("Conditions") +
    ylab("Frequency") +
    ggtitle(
      "Open chromatin region\ngenomic location"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = 12,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = 12
      ),
      axis.title = element_text(
        size = 12
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 12,
        margin = margin(
          r = 10,
          unit = "pt"
        )
      ),
      strip.text.x = element_text(
        size = 12,
        face = "bold"
      ),
      legend.box.spacing = grid::unit(
        0,
        "mm"
      )
    )

  return(annotation_plot)
}


###########################################
# Generate plots
###########################################

log_message("Generating correlation heatmap.")

corr_plot <- correlation(
  norm_data,
  corr_method,
  conditions_color
)

log_message("Generating PCA plot.")

pca_plot <- plot_pca(
  norm_data,
  conditions_color,
  c("PC1", "PC2")
)

log_message("Generating open chromatin region overlap plot.")

venn_plot <- plot_venn(
  counts,
  conditions_color
)

log_message("Generating genomic annotation plot.")

anno_plot <- plot_anno(
  counts,
  conditions_color
)


###########################################
# Assemble figure
###########################################

log_message("Assembling QC figure.")

figure <- cowplot::plot_grid(
  plotlist = list(
    corr_plot,
    pca_plot[[1]],
    venn_plot,
    anno_plot
  ),
  labels = "AUTO",
  ncol = 2,
  align = "h"
)


###########################################
# Export plots
###########################################

dir.create(
  dirname(output_pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(output_png),
  recursive = TRUE,
  showWarnings = FALSE
)

log_message("Saving PDF figure: ", output_pdf)

cowplot::save_plot(
  output_pdf,
  figure,
  base_width = 30,
  base_height = 20,
  units = "cm",
  dpi = 300
)

log_message("Saving PNG figure: ", output_png)

cowplot::save_plot(
  output_png,
  figure,
  base_width = 30,
  base_height = 20,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("QC figure generation completed.")