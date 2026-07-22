source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting TRIM28 multi-TF enrichment in DAR clusters using all ATAC peaks as reference")

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

ATAC_all <- rtracklayer::import(snakemake@input[["ATAC_all"]])

ATAC_clustering <- read.table(
  snakemake@input[["ATAC_clustering"]],
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(ATAC_clustering) < 2) {
  log_error("ATAC_clustering must contain at least two columns: region and cluster")
}

colnames(ATAC_clustering)[1:2] <- c(
  "region",
  "cluster"
)

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

log_message("All ATAC peaks: ", length(ATAC_all))
log_message("DAR cluster a regions: ", length(DAR_a))
log_message("DAR cluster b regions: ", length(DAR_b))
log_message("TRIM28 peaks: ", length(TRIM28_peaks))

for (TF_name in names(TF_list)) {
  log_message(TF_name, " peaks: ", length(TF_list[[TF_name]]))
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

annotate_ATAC_with_TRIM28_TFs <- function(
  ATAC_regions,
  TRIM28_peaks,
  TF_list
) {
  log_message("Annotating all ATAC peaks with TRIM28 and TF overlaps")

  overlap_TRIM28 <- count_overlap_binary(
    query = ATAC_regions,
    subject = TRIM28_peaks
  )

  TF_overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      count_overlap_binary(
        query = ATAC_regions,
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
        return("-")
      }

      paste(TFs, collapse = "+")
    }
  )

  data.frame(
    chromosome = as.character(GenomicRanges::seqnames(ATAC_regions)),
    start = GenomicRanges::start(ATAC_regions),
    end = GenomicRanges::end(ATAC_regions),
    region = make_region_id(ATAC_regions),
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
}

add_DAR_clusters <- function(
  ATAC_table,
  ATAC_regions,
  DAR_a,
  DAR_b
) {
  log_message("Adding DAR cluster overlap to all ATAC peaks")

  ATAC_table$DAR_a <- count_overlap_binary(
    query = ATAC_regions,
    subject = DAR_a
  )

  ATAC_table$DAR_b <- count_overlap_binary(
    query = ATAC_regions,
    subject = DAR_b
  )

  ATAC_table$DAR_cluster <- "-"

  ATAC_table$DAR_cluster[
    ATAC_table$DAR_a
  ] <- "a"

  ATAC_table$DAR_cluster[
    ATAC_table$DAR_b
  ] <- "b"

  ATAC_table
}

summarise_TRIM28_TF_classes <- function(
  ATAC_table,
  response_column,
  response_label
) {
  log_message("Summarising ", response_label, " by TRIM28 / TF class")

  trim28_regions <- ATAC_table[
    ATAC_table$overlap_TRIM28,
    ,
    drop = FALSE
  ]

  trim28_regions$TRIM28_TF_class <- ifelse(
    trim28_regions$n_TFs == 0,
    "TRIM28 only",
    paste0("TRIM28 + ", trim28_regions$n_TFs, " TF")
  )

  trim28_regions$TRIM28_TF_class[
    trim28_regions$n_TFs > 1
  ] <- paste0(
    "TRIM28 + ",
    trim28_regions$n_TFs[trim28_regions$n_TFs > 1],
    " TFs"
  )

  class_order <- c(
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )

  summary_list <- lapply(
    class_order,
    function(class_name) {
      subset_table <- trim28_regions[
        trim28_regions$TRIM28_TF_class == class_name,
        ,
        drop = FALSE
      ]

      n_total <- nrow(subset_table)
      n_response <- sum(subset_table[[response_column]])

      data.frame(
        analysis = response_label,
        TRIM28_TF_class = class_name,
        n_total = n_total,
        n_response = n_response,
        pct_response = ifelse(
          n_total > 0,
          100 * n_response / n_total,
          NA_real_
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  summary_table <- do.call(
    rbind,
    summary_list
  )

  summary_table$TRIM28_TF_class <- factor(
    summary_table$TRIM28_TF_class,
    levels = class_order
  )

  summary_table
}

run_fisher_vs_all_ATAC <- function(
  ATAC_table,
  category_column,
  response_column,
  label
) {
  log_message("Running Fisher test for ", label)

  in_category <- ATAC_table[[category_column]]
  in_response <- ATAC_table[[response_column]]

  contingency_table <- matrix(
    c(
      sum(in_category & in_response),
      sum(in_category & !in_response),
      sum(!in_category & in_response),
      sum(!in_category & !in_response)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(contingency_table) <- c(
    category_column,
    paste0("not_", category_column)
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
    comparison = label,
    category = category_column,
    response = response_column,
    n_category_response = contingency_table[1, 1],
    n_category_total = sum(contingency_table[1, ]),
    pct_category_response = 100 * contingency_table[1, 1] / sum(contingency_table[1, ]),
    n_background_response = contingency_table[2, 1],
    n_background_total = sum(contingency_table[2, ]),
    pct_background_response = 100 * contingency_table[2, 1] / sum(contingency_table[2, ]),
    odds_ratio = unname(fisher_result$estimate),
    p_value = fisher_result$p.value,
    stringsAsFactors = FALSE
  )
}

run_fisher_a_vs_b <- function(
  ATAC_table,
  category_column,
  label
) {
  log_message("Running DAR a vs DAR b Fisher test for ", label)

  DAR_only <- ATAC_table[
    ATAC_table$DAR_a | ATAC_table$DAR_b,
    ,
    drop = FALSE
  ]

  in_category <- DAR_only[[category_column]]
  is_a <- DAR_only$DAR_a
  is_b <- DAR_only$DAR_b

  contingency_table <- matrix(
    c(
      sum(is_a & in_category),
      sum(is_a & !in_category),
      sum(is_b & in_category),
      sum(is_b & !in_category)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(contingency_table) <- c(
    "DAR_a",
    "DAR_b"
  )

  colnames(contingency_table) <- c(
    category_column,
    paste0("not_", category_column)
  )

  fisher_result <- stats::fisher.test(
    contingency_table,
    alternative = "greater"
  )

  data.frame(
    comparison = label,
    category = category_column,
    response = "DAR_a_vs_DAR_b",
    n_category_response = contingency_table[1, 1],
    n_category_total = sum(contingency_table[1, ]),
    pct_category_response = 100 * contingency_table[1, 1] / sum(contingency_table[1, ]),
    n_background_response = contingency_table[2, 1],
    n_background_total = sum(contingency_table[2, ]),
    pct_background_response = 100 * contingency_table[2, 1] / sum(contingency_table[2, ]),
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

plot_response_by_TF_number <- function(
  class_summary,
  fisher_summary
) {
  fisher_label <- paste(
    paste0(
      fisher_summary$comparison,
      ": OR = ",
      round(fisher_summary$odds_ratio, 2),
      ", ",
      vapply(
        fisher_summary$p_value,
        format_pvalue,
        character(1)
      )
    ),
    collapse = "\n"
  )

  ggplot(
    class_summary,
    aes(
      x = TRIM28_TF_class,
      y = pct_response
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
          round(pct_response, 1),
          "%\n(",
          n_response,
          "/",
          n_total,
          ")"
        )
      ),
      vjust = -0.25,
      size = 4
    ) +
    annotate(
      "text",
      x = 3.5,
      y = max(class_summary$pct_response, na.rm = TRUE) * 1.18,
      label = fisher_label,
      hjust = 0.5,
      size = 4
    ) +
    scale_y_continuous(
      limits = c(
        0,
        max(class_summary$pct_response, na.rm = TRUE) * 1.35
      ),
      labels = function(x) paste0(x, "%")
    ) +
    xlab(NULL) +
    ylab(unique(class_summary$analysis)) +
    ggtitle(
      paste0(unique(class_summary$analysis), " among TRIM28 regions")
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
      axis.text = element_text(size = 12)
    )
}

###########################################
# Annotate all ATAC peaks
###########################################

ATAC_table <- annotate_ATAC_with_TRIM28_TFs(
  ATAC_regions = ATAC_all,
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

ATAC_table <- add_DAR_clusters(
  ATAC_table = ATAC_table,
  ATAC_regions = ATAC_all,
  DAR_a = DAR_a,
  DAR_b = DAR_b
)

###########################################
# Summary and Fisher tests
###########################################

class_summary_a <- summarise_TRIM28_TF_classes(
  ATAC_table = ATAC_table,
  response_column = "DAR_a",
  response_label = "DAR cluster a overlap"
)

class_summary_b <- summarise_TRIM28_TF_classes(
  ATAC_table = ATAC_table,
  response_column = "DAR_b",
  response_label = "DAR cluster b overlap"
)

fisher_summary_ATAC_reference <- do.call(
  rbind,
  list(
    run_fisher_vs_all_ATAC(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_2TFs",
      response_column = "DAR_a",
      label = "TRIM28 + >=2 TFs / DAR a vs all ATAC"
    ),
    run_fisher_vs_all_ATAC(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_3TFs",
      response_column = "DAR_a",
      label = "TRIM28 + >=3 TFs / DAR a vs all ATAC"
    ),
    run_fisher_vs_all_ATAC(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_2TFs",
      response_column = "DAR_b",
      label = "TRIM28 + >=2 TFs / DAR b vs all ATAC"
    ),
    run_fisher_vs_all_ATAC(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_3TFs",
      response_column = "DAR_b",
      label = "TRIM28 + >=3 TFs / DAR b vs all ATAC"
    )
  )
)

fisher_summary_a_vs_b <- do.call(
  rbind,
  list(
    run_fisher_a_vs_b(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_2TFs",
      label = "DAR a enriched over DAR b for TRIM28 + >=2 TFs"
    ),
    run_fisher_a_vs_b(
      ATAC_table = ATAC_table,
      category_column = "TRIM28_plus_3TFs",
      label = "DAR a enriched over DAR b for TRIM28 + >=3 TFs"
    )
  )
)

fisher_summary <- rbind(
  fisher_summary_ATAC_reference,
  fisher_summary_a_vs_b
)

class_summary <- rbind(
  class_summary_a,
  class_summary_b
)

###########################################
# Save tables
###########################################

write.table(
  ATAC_table,
  file = snakemake@output[["regions"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  class_summary,
  file = snakemake@output[["class_summary"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  fisher_summary,
  file = snakemake@output[["fisher_summary"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Save figure
###########################################

figure_a <- plot_response_by_TF_number(
  class_summary = class_summary_a,
  fisher_summary = fisher_summary_ATAC_reference[
    fisher_summary_ATAC_reference$response == "DAR_a",
    ,
    drop = FALSE
  ]
)

figure_b <- plot_response_by_TF_number(
  class_summary = class_summary_b,
  fisher_summary = fisher_summary_ATAC_reference[
    fisher_summary_ATAC_reference$response == "DAR_b",
    ,
    drop = FALSE
  ]
)

figure <- cowplot::plot_grid(
  figure_a,
  figure_b,
  ncol = 1,
  labels = c("a", "b")
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 16,
  base_height = 22,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 16,
  base_height = 22,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")