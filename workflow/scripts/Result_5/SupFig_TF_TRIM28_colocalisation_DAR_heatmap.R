source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting DAR multi-TF composition analysis")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("ggplot2")
  library("cowplot")
  library("scales")
})

###########################################
# Load inputs
###########################################

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28"]]
)

TF_list <- list(
  FOXL2 = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2 = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2 = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX = rtracklayer::import(snakemake@input[["RUNX"]])
)

ATAC_clustering <- read.csv(
  snakemake@input[["ATAC_clustering"]],
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("region", "cluster") %in% colnames(ATAC_clustering))) {
  log_error("ATAC_clustering must contain 'region' and 'cluster' columns")
}

if (!all(ATAC_clustering$cluster %in% c("a", "b"))) {
  log_error("ATAC_clustering cluster column must contain only 'a' and/or 'b'")
}

DAR_regions <- GenomicRanges::GRanges(
  ATAC_clustering$region
)

DAR_regions$cluster <- ATAC_clustering$cluster

top_n <- if ("top_n" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["top_n"]])
} else {
  15
}

if (is.na(top_n) || top_n < 1) {
  log_error("Parameter 'top_n' must be a positive numeric value")
}

###########################################
# Functions
###########################################

make_region_id <- function(regions) {
  paste0(
    as.character(GenomicRanges::seqnames(regions)),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

count_overlap_binary <- function(query, subject) {
  GenomicRanges::countOverlaps(
    query,
    subject,
    ignore.strand = TRUE
  ) > 0
}

annotate_DAR_with_TRIM28_TFs <- function(
  DAR_regions,
  TRIM28_peaks,
  TF_list
) {
  log_message("Annotating DAR regions with TRIM28 and TF overlaps")

  overlap_TRIM28 <- count_overlap_binary(
    query = DAR_regions,
    subject = TRIM28_peaks
  )

  TF_overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      count_overlap_binary(
        query = DAR_regions,
        subject = TF_list[[TF_name]]
      )
    }
  )

  n_TFs <- rowSums(
    TF_overlap_matrix
  )

  TF_combination <- apply(
    TF_overlap_matrix,
    1,
    function(x) {
      TFs <- names(x)[x]

      if (length(TFs) == 0) {
        return("No TF")
      }

      paste(TFs, collapse = " + ")
    }
  )

  DAR_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(DAR_regions)),
    start = GenomicRanges::start(DAR_regions),
    end = GenomicRanges::end(DAR_regions),
    region = make_region_id(DAR_regions),
    DAR_cluster = DAR_regions$cluster,
    overlap_TRIM28 = overlap_TRIM28,
    TF_overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = TF_combination,
    TRIM28_only = overlap_TRIM28 & n_TFs == 0,
    TRIM28_plus_1TF = overlap_TRIM28 & n_TFs == 1,
    TRIM28_plus_2TFs = overlap_TRIM28 & n_TFs >= 2,
    TRIM28_plus_3TFs = overlap_TRIM28 & n_TFs >= 3,
    TRIM28_plus_4TFs = overlap_TRIM28 & n_TFs >= 4,
    stringsAsFactors = FALSE
  )

  DAR_table$TRIM28_TF_class <- ifelse(
    !DAR_table$overlap_TRIM28,
    "No TRIM28",
    ifelse(
      DAR_table$n_TFs == 0,
      "TRIM28 only",
      paste0("TRIM28 + ", DAR_table$n_TFs, " TF")
    )
  )

  DAR_table$TRIM28_TF_class[
    DAR_table$overlap_TRIM28 & DAR_table$n_TFs > 1
  ] <- paste0(
    "TRIM28 + ",
    DAR_table$n_TFs[DAR_table$overlap_TRIM28 & DAR_table$n_TFs > 1],
    " TFs"
  )

  DAR_table
}

make_composition_heatmap_table <- function(
  DAR_table,
  TF_names
) {
  log_message("Preparing DAR TF composition heatmap table")

  class_order <- c(
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )

  DAR_table <- DAR_table[
    DAR_table$overlap_TRIM28,
    ,
    drop = FALSE
  ]

  heatmap_table <- do.call(
    rbind,
    lapply(
      sort(unique(DAR_table$DAR_cluster)),
      function(cluster_id) {
        do.call(
          rbind,
          lapply(
            class_order,
            function(class_name) {
              subset_table <- DAR_table[
                DAR_table$DAR_cluster == cluster_id &
                  DAR_table$TRIM28_TF_class == class_name,
                ,
                drop = FALSE
              ]

              n_regions <- nrow(subset_table)

              do.call(
                rbind,
                lapply(
                  TF_names,
                  function(TF_name) {
                    n_with_TF <- if (n_regions == 0) {
                      0
                    } else {
                      sum(subset_table[[TF_name]])
                    }

                    data.frame(
                      DAR_cluster = cluster_id,
                      TRIM28_TF_class = class_name,
                      TF = TF_name,
                      n_regions = n_regions,
                      n_with_TF = n_with_TF,
                      pct_with_TF = ifelse(
                        n_regions > 0,
                        100 * n_with_TF / n_regions,
                        NA_real_
                      ),
                      stringsAsFactors = FALSE
                    )
                  }
                )
              )
            }
          )
        )
      }
    )
  )

  heatmap_table$TRIM28_TF_class <- factor(
    heatmap_table$TRIM28_TF_class,
    levels = class_order
  )

  heatmap_table$TF <- factor(
    heatmap_table$TF,
    levels = TF_names
  )

  heatmap_table
}

make_combination_summary <- function(
  DAR_table,
  top_n
) {
  log_message("Summarising most frequent TF combinations in DAR regions")

  DAR_table <- DAR_table[
    DAR_table$overlap_TRIM28 &
      DAR_table$TF_combination != "No TF",
    ,
    drop = FALSE
  ]

  combination_table <- as.data.frame(
    table(
      DAR_table$DAR_cluster,
      DAR_table$TF_combination
    ),
    stringsAsFactors = FALSE
  )

  colnames(combination_table) <- c(
    "DAR_cluster",
    "TF_combination",
    "n_regions"
  )

  combination_table <- combination_table[
    combination_table$n_regions > 0,
    ,
    drop = FALSE
  ]

  total_by_cluster <- aggregate(
    n_regions ~ DAR_cluster,
    data = combination_table,
    FUN = sum
  )

  colnames(total_by_cluster) <- c(
    "DAR_cluster",
    "n_TRIM28_TF_regions"
  )

  combination_table <- merge(
    combination_table,
    total_by_cluster,
    by = "DAR_cluster",
    all.x = TRUE,
    sort = FALSE
  )

  combination_table$pct_TRIM28_TF_regions <- 100 *
    combination_table$n_regions /
    combination_table$n_TRIM28_TF_regions

  combination_table <- combination_table[
    order(
      combination_table$DAR_cluster,
      -combination_table$n_regions
    ),
    ,
    drop = FALSE
  ]

  top_combination_table <- do.call(
    rbind,
    lapply(
      split(combination_table, combination_table$DAR_cluster),
      function(x) {
        x <- x[
          order(x$n_regions, decreasing = TRUE),
          ,
          drop = FALSE
        ]

        head(x, top_n)
      }
    )
  )

  list(
    full = combination_table,
    top = top_combination_table
  )
}

plot_composition_heatmap <- function(
  heatmap_table,
  cluster_id
) {
  plot_table <- heatmap_table[
    heatmap_table$DAR_cluster == cluster_id,
    ,
    drop = FALSE
  ]

  ggplot(
    plot_table,
    aes(
      x = TRIM28_TF_class,
      y = TF,
      fill = pct_with_TF
    )
  ) +
    geom_tile(
      colour = "white",
      linewidth = 0.6
    ) +
    geom_text(
      aes(
        label = ifelse(
          is.na(pct_with_TF),
          "",
          paste0(round(pct_with_TF, 0), "%")
        )
      ),
      size = 4
    ) +
    scale_fill_gradient(
      low = "#FFF9EC",
      high = "#7D1128",
      na.value = "grey90",
      limits = c(0, 100),
      labels = function(x) paste0(x, "%")
    ) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle(
      paste0("TF composition of TRIM28-positive DAR ", cluster_id, " regions")
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = 12,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = 12,
        angle = 35,
        hjust = 1
      ),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "right"
    ) +
    labs(
      fill = "% regions\nwith TF"
    )
}

plot_top_combinations <- function(
  top_combination_table,
  cluster_id
) {
  plot_table <- top_combination_table[
    top_combination_table$DAR_cluster == cluster_id,
    ,
    drop = FALSE
  ]

  plot_table <- plot_table[
    order(plot_table$n_regions, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  plot_table$TF_combination <- factor(
    plot_table$TF_combination,
    levels = rev(plot_table$TF_combination)
  )

  ggplot(
    plot_table,
    aes(
      x = TF_combination,
      y = n_regions
    )
  ) +
    geom_col(
      fill = "#7D1128",
      colour = "#7D1128",
      linewidth = 0.3,
      width = 0.75
    ) +
    geom_text(
      aes(
        label = paste0(
          scales::comma(n_regions),
          " (",
          round(pct_TRIM28_TF_regions, 1),
          "%)"
        )
      ),
      hjust = -0.05,
      size = 4
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.25))
    ) +
    xlab(NULL) +
    ylab("Number of DAR regions") +
    ggtitle(
      paste0("Most frequent TF combinations in DAR ", cluster_id)
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = 12,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(size = 12)
    )
}

###########################################
# Run analysis
###########################################

DAR_table <- annotate_DAR_with_TRIM28_TFs(
  DAR_regions = DAR_regions,
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

heatmap_table <- make_composition_heatmap_table(
  DAR_table = DAR_table,
  TF_names = names(TF_list)
)

combination_summary <- make_combination_summary(
  DAR_table = DAR_table,
  top_n = top_n
)

###########################################
# Save tables
###########################################

write.table(
  DAR_table,
  file = snakemake@output[["regions"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  heatmap_table,
  file = snakemake@output[["heatmap_table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  combination_summary$full,
  file = snakemake@output[["combination_summary"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Save figure
###########################################

figure_heatmap_a <- plot_composition_heatmap(
  heatmap_table = heatmap_table,
  cluster_id = "a"
)

figure_heatmap_b <- plot_composition_heatmap(
  heatmap_table = heatmap_table,
  cluster_id = "b"
)

figure_combinations_a <- plot_top_combinations(
  top_combination_table = combination_summary$top,
  cluster_id = "a"
)

figure_combinations_b <- plot_top_combinations(
  top_combination_table = combination_summary$top,
  cluster_id = "b"
)

figure <- cowplot::plot_grid(
  figure_heatmap_a,
  figure_heatmap_b,
  figure_combinations_a,
  figure_combinations_b,
  ncol = 1,
  rel_heights = c(1, 1, 1.2, 1.2),
  labels = c("a", "b", "c", "d")
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 36,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 36,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")