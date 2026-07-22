source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting genome-wide TRIM28 multi-TF co-occupation analysis")

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
  ESR2  = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX = rtracklayer::import(snakemake@input[["RUNX"]])
)

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

make_TRIM28_multiTF_table <- function(
  TRIM28_peaks,
  TF_list
) {
  log_message("Annotating TRIM28 peaks with ovarian TF overlaps")

  TF_overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      count_overlap_binary(
        query = TRIM28_peaks,
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

  TRIM28_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(TRIM28_peaks)),
    start = GenomicRanges::start(TRIM28_peaks),
    end = GenomicRanges::end(TRIM28_peaks),
    region = make_region_id(TRIM28_peaks),
    TF_overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = TF_combination,
    TRIM28_only = n_TFs == 0,
    TRIM28_plus_1TF = n_TFs == 1,
    TRIM28_plus_2TFs = n_TFs >= 2,
    TRIM28_plus_3TFs = n_TFs >= 3,
    TRIM28_plus_4TFs = n_TFs >= 4,
    stringsAsFactors = FALSE
  )

  TRIM28_table$TRIM28_TF_class <- ifelse(
    TRIM28_table$n_TFs == 0,
    "TRIM28 only",
    paste0("TRIM28 + ", TRIM28_table$n_TFs, " TF")
  )

  TRIM28_table$TRIM28_TF_class[
    TRIM28_table$n_TFs > 1
  ] <- paste0(
    "TRIM28 + ",
    TRIM28_table$n_TFs[TRIM28_table$n_TFs > 1],
    " TFs"
  )

  TRIM28_table
}

summarise_TRIM28_multiTF_classes <- function(
  TRIM28_table
) {
  log_message("Summarising TRIM28 co-occupation classes")

  class_order <- c(
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )

  summary_table <- do.call(
    rbind,
    lapply(
      class_order,
      function(class_name) {
        n_regions <- sum(
          TRIM28_table$TRIM28_TF_class == class_name
        )

        data.frame(
          TRIM28_TF_class = class_name,
          n_regions = n_regions,
          pct_TRIM28_peaks = 100 * n_regions / nrow(TRIM28_table),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  summary_table$TRIM28_TF_class <- factor(
    summary_table$TRIM28_TF_class,
    levels = class_order
  )

  summary_table
}

summarise_TF_combinations <- function(
  TRIM28_table
) {
  log_message("Summarising TF combinations at TRIM28 peaks")

  combination_table <- as.data.frame(
    table(TRIM28_table$TF_combination),
    stringsAsFactors = FALSE
  )

  colnames(combination_table) <- c(
    "TF_combination",
    "n_regions"
  )

  combination_table$pct_TRIM28_peaks <- 100 *
    combination_table$n_regions /
    nrow(TRIM28_table)

  combination_table$n_TFs <- vapply(
    strsplit(as.character(combination_table$TF_combination), " \\+ "),
    function(x) {
      if (identical(x, "No TF")) {
        return(0)
      }

      length(x)
    },
    numeric(1)
  )

  combination_table <- combination_table[
    order(
      combination_table$n_TFs,
      combination_table$n_regions,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  combination_table
}

plot_TRIM28_multiTF_classes <- function(
  summary_table
) {
  ggplot(
    summary_table,
    aes(
      x = TRIM28_TF_class,
      y = pct_TRIM28_peaks
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
          round(pct_TRIM28_peaks, 1),
          "%\n(",
          scales::comma(n_regions),
          ")"
        )
      ),
      vjust = -0.25,
      size = 4
    ) +
    scale_y_continuous(
      limits = c(
        0,
        max(summary_table$pct_TRIM28_peaks, na.rm = TRUE) * 1.25
      ),
      labels = function(x) paste0(x, "%")
    ) +
    xlab(NULL) +
    ylab("Proportion of TRIM28 peaks") +
    ggtitle(
      "Genome-wide TRIM28 co-occupation with ovarian TFs"
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
      axis.text.y = element_text(size = 12)
    )
}

plot_top_multiTF_combinations <- function(
  combination_table,
  top_n = 15
) {
  plot_table <- combination_table[
    combination_table$n_TFs >= 2,
    ,
    drop = FALSE
  ]

  plot_table <- plot_table[
    order(plot_table$n_regions, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  plot_table <- head(
    plot_table,
    top_n
  )

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
      width = 0.75
    ) +
    geom_text(
      aes(
        label = paste0(
          scales::comma(n_regions),
          " (",
          round(pct_TRIM28_peaks, 1),
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
    ylab("Number of TRIM28 peaks") +
    ggtitle(
      "Most frequent multi-TF combinations at TRIM28 peaks"
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

TRIM28_table <- make_TRIM28_multiTF_table(
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

summary_table <- summarise_TRIM28_multiTF_classes(
  TRIM28_table = TRIM28_table
)

combination_table <- summarise_TF_combinations(
  TRIM28_table = TRIM28_table
)

multiTF_table <- TRIM28_table[
  TRIM28_table$n_TFs >= 2,
  ,
  drop = FALSE
]

log_message("TRIM28 peaks with at least 2 TFs: ", nrow(multiTF_table))
log_message(
  "Percentage of TRIM28 peaks with at least 2 TFs: ",
  round(100 * nrow(multiTF_table) / nrow(TRIM28_table), 2),
  "%"
)

###########################################
# Save tables
###########################################

write.table(
  multiTF_table,
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Save figure
###########################################

figure_classes <- plot_TRIM28_multiTF_classes(
  summary_table = summary_table
)

figure_combinations <- plot_top_multiTF_combinations(
  combination_table = combination_table,
  top_n = 15
)

figure <- cowplot::plot_grid(
  figure_classes,
  figure_combinations,
  ncol = 1,
  rel_heights = c(1, 1.3),
  labels = c("a", "b")
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 22,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 22,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")