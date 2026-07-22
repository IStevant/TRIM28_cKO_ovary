source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting ATAC multi-TF TRIM28 dependency analysis")

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

ATAC_all <- rtracklayer::import(
  snakemake@input[["ATAC_all"]]
)

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

DAR_a <- DAR_regions[
  DAR_regions$cluster == "a"
]

DAR_b <- DAR_regions[
  DAR_regions$cluster == "b"
]

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

annotate_ATAC_regions <- function(
  ATAC_all,
  TRIM28_peaks,
  TF_list,
  DAR_a,
  DAR_b
) {
  log_message("Annotating all ATAC peaks with TRIM28, TFs and DAR clusters")

  overlap_TRIM28 <- count_overlap_binary(
    query = ATAC_all,
    subject = TRIM28_peaks
  )

  TF_overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      count_overlap_binary(
        query = ATAC_all,
        subject = TF_list[[TF_name]]
      )
    }
  )

  n_TFs <- rowSums(TF_overlap_matrix)

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

  DAR_a_overlap <- count_overlap_binary(
    query = ATAC_all,
    subject = DAR_a
  )

  DAR_b_overlap <- count_overlap_binary(
    query = ATAC_all,
    subject = DAR_b
  )

  ATAC_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(ATAC_all)),
    start = GenomicRanges::start(ATAC_all),
    end = GenomicRanges::end(ATAC_all),
    region = make_region_id(ATAC_all),
    overlap_TRIM28 = overlap_TRIM28,
    TF_overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = TF_combination,
    multiTF_2 = n_TFs >= 2,
    multiTF_3 = n_TFs >= 3,
    DAR_a = DAR_a_overlap,
    DAR_b = DAR_b_overlap,
    DAR_cluster = ifelse(
      DAR_a_overlap,
      "a",
      ifelse(DAR_b_overlap, "b", "-")
    ),
    stringsAsFactors = FALSE
  )

  ATAC_table$multiTF2_TRIM28_status <- ifelse(
    ATAC_table$n_TFs >= 2 & ATAC_table$overlap_TRIM28,
    ">=2 TFs + TRIM28",
    ifelse(
      ATAC_table$n_TFs >= 2 & !ATAC_table$overlap_TRIM28,
      ">=2 TFs without TRIM28",
      "Other ATAC"
    )
  )

  ATAC_table$multiTF3_TRIM28_status <- ifelse(
    ATAC_table$n_TFs >= 3 & ATAC_table$overlap_TRIM28,
    ">=3 TFs + TRIM28",
    ifelse(
      ATAC_table$n_TFs >= 3 & !ATAC_table$overlap_TRIM28,
      ">=3 TFs without TRIM28",
      "Other ATAC"
    )
  )

  ATAC_table
}

summarise_closing_by_status <- function(
  ATAC_table,
  status_column,
  analysis_label
) {
  log_message("Summarising closing frequency for ", analysis_label)

  keep <- ATAC_table[[status_column]] != "Other ATAC"

  summary_table <- do.call(
    rbind,
    lapply(
      unique(ATAC_table[[status_column]][keep]),
      function(status) {
        subset_table <- ATAC_table[
          ATAC_table[[status_column]] == status,
          ,
          drop = FALSE
        ]

        n_total <- nrow(subset_table)
        n_closing <- sum(subset_table$DAR_a)
        n_opening <- sum(subset_table$DAR_b)

        data.frame(
          analysis = analysis_label,
          group = status,
          n_total = n_total,
          n_closing_DAR_a = n_closing,
          pct_closing_DAR_a = 100 * n_closing / n_total,
          n_opening_DAR_b = n_opening,
          pct_opening_DAR_b = 100 * n_opening / n_total,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  summary_table
}

run_TRIM28_dependency_fisher <- function(
  ATAC_table,
  status_column,
  with_TRIM28_label,
  without_TRIM28_label,
  response_column,
  response_label,
  analysis_label
) {
  log_message("Running Fisher test for ", analysis_label, " / ", response_label)

  test_table <- ATAC_table[
    ATAC_table[[status_column]] %in% c(
      with_TRIM28_label,
      without_TRIM28_label
    ),
    ,
    drop = FALSE
  ]

  with_TRIM28 <- test_table[[status_column]] == with_TRIM28_label
  without_TRIM28 <- test_table[[status_column]] == without_TRIM28_label
  response <- test_table[[response_column]]

  contingency_table <- matrix(
    c(
      sum(with_TRIM28 & response),
      sum(with_TRIM28 & !response),
      sum(without_TRIM28 & response),
      sum(without_TRIM28 & !response)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(contingency_table) <- c(
    with_TRIM28_label,
    without_TRIM28_label
  )

  colnames(contingency_table) <- c(
    response_column,
    paste0("not_", response_column)
  )

  fisher_result <- stats::fisher.test(
    contingency_table,
    alternative = "greater"
  )

  data.frame(
    analysis = analysis_label,
    response = response_label,
    with_TRIM28_group = with_TRIM28_label,
    without_TRIM28_group = without_TRIM28_label,
    n_with_TRIM28_response = contingency_table[1, 1],
    n_with_TRIM28_total = sum(contingency_table[1, ]),
    pct_with_TRIM28_response = 100 * contingency_table[1, 1] / sum(contingency_table[1, ]),
    n_without_TRIM28_response = contingency_table[2, 1],
    n_without_TRIM28_total = sum(contingency_table[2, ]),
    pct_without_TRIM28_response = 100 * contingency_table[2, 1] / sum(contingency_table[2, ]),
    odds_ratio = unname(fisher_result$estimate),
    p_value = fisher_result$p.value,
    stringsAsFactors = FALSE
  )
}

format_pvalue <- function(p_value) {
  if (is.na(p_value)) {
    return("P = NA")
  }

  if (p_value < 2.2e-16) {
    return("P < 2.2e-16")
  }

  paste0(
    "P = ",
    formatC(
      p_value,
      format = "e",
      digits = 2
    )
  )
}

plot_dependency <- function(
  summary_table,
  fisher_table,
  response = "closing"
) {
  if (response == "closing") {
    y_column <- "pct_closing_DAR_a"
    n_column <- "n_closing_DAR_a"
    y_label <- "ATAC peaks closing in cKO (%)"
    title <- "TRIM28 increases closing of multi-TF ATAC peaks"
    fisher_response <- "DAR a closing"
  } else {
    y_column <- "pct_opening_DAR_b"
    n_column <- "n_opening_DAR_b"
    y_label <- "ATAC peaks opening in cKO (%)"
    title <- "Opening of multi-TF ATAC peaks by TRIM28 status"
    fisher_response <- "DAR b opening"
  }

  fisher_labels <- fisher_table[
    fisher_table$response == fisher_response,
    ,
    drop = FALSE
  ]

  fisher_label <- paste(
    paste0(
      fisher_labels$analysis,
      ": OR = ",
      round(fisher_labels$odds_ratio, 2),
      ", ",
      vapply(fisher_labels$p_value, format_pvalue, character(1))
    ),
    collapse = "\n"
  )

  plot_table <- summary_table

  ggplot2::ggplot(
    plot_table,
    ggplot2::aes(
      x = group,
      y = .data[[y_column]]
    )
  ) +
    ggplot2::geom_col(
      fill = "#320A28",
      colour = "#320A28",
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          round(.data[[y_column]], 1),
          "%\n(",
          .data[[n_column]],
          "/",
          n_total,
          ")"
        )
      ),
      vjust = -0.25,
      size = 3.5
    ) +
    ggplot2::facet_wrap(
      ~analysis,
      scales = "free_x"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = max(plot_table[[y_column]], na.rm = TRUE) * 1.25,
      label = fisher_label,
      size = 3.2,
      hjust = 0.5
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        0,
        max(plot_table[[y_column]], na.rm = TRUE) * 1.45
      ),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(y_label) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_light() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 12,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        size = 11,
        angle = 30,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(size = 11),
      strip.text = ggplot2::element_text(
        size = 11,
        face = "bold"
      )
    )
}

###########################################
# Run analysis
###########################################

ATAC_table <- annotate_ATAC_regions(
  ATAC_all = ATAC_all,
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list,
  DAR_a = DAR_a,
  DAR_b = DAR_b
)

summary_multiTF2 <- summarise_closing_by_status(
  ATAC_table = ATAC_table,
  status_column = "multiTF2_TRIM28_status",
  analysis_label = "ATAC peaks with >=2 TFs"
)

summary_multiTF3 <- summarise_closing_by_status(
  ATAC_table = ATAC_table,
  status_column = "multiTF3_TRIM28_status",
  analysis_label = "ATAC peaks with >=3 TFs"
)

summary_table <- rbind(
  summary_multiTF2,
  summary_multiTF3
)

fisher_table <- rbind(
  run_TRIM28_dependency_fisher(
    ATAC_table = ATAC_table,
    status_column = "multiTF2_TRIM28_status",
    with_TRIM28_label = ">=2 TFs + TRIM28",
    without_TRIM28_label = ">=2 TFs without TRIM28",
    response_column = "DAR_a",
    response_label = "DAR a closing",
    analysis_label = "ATAC peaks with >=2 TFs"
  ),
  run_TRIM28_dependency_fisher(
    ATAC_table = ATAC_table,
    status_column = "multiTF3_TRIM28_status",
    with_TRIM28_label = ">=3 TFs + TRIM28",
    without_TRIM28_label = ">=3 TFs without TRIM28",
    response_column = "DAR_a",
    response_label = "DAR a closing",
    analysis_label = "ATAC peaks with >=3 TFs"
  ),
  run_TRIM28_dependency_fisher(
    ATAC_table = ATAC_table,
    status_column = "multiTF2_TRIM28_status",
    with_TRIM28_label = ">=2 TFs + TRIM28",
    without_TRIM28_label = ">=2 TFs without TRIM28",
    response_column = "DAR_b",
    response_label = "DAR b opening",
    analysis_label = "ATAC peaks with >=2 TFs"
  ),
  run_TRIM28_dependency_fisher(
    ATAC_table = ATAC_table,
    status_column = "multiTF3_TRIM28_status",
    with_TRIM28_label = ">=3 TFs + TRIM28",
    without_TRIM28_label = ">=3 TFs without TRIM28",
    response_column = "DAR_b",
    response_label = "DAR b opening",
    analysis_label = "ATAC peaks with >=3 TFs"
  )
)

###########################################
# Save tables
###########################################

write.table(
  ATAC_table,
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# write.table(
#   summary_table,
#   file = snakemake@output[["summary"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

# write.table(
#   fisher_table,
#   file = snakemake@output[["fisher"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

###########################################
# Save figures
###########################################

figure_closing <- plot_dependency(
  summary_table = summary_table,
  fisher_table = fisher_table,
  response = "closing"
)

figure_opening <- plot_dependency(
  summary_table = summary_table,
  fisher_table = fisher_table,
  response = "opening"
)

figure <- cowplot::plot_grid(
  figure_closing,
  figure_opening,
  ncol = 1,
  labels = c("a", "b")
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 18,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 18,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")