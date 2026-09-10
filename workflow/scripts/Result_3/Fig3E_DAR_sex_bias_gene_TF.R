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

log_message("Starting ATAC DAR summary analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
  library("scales")
})

###########################################
# Load input data
###########################################

log_message("Loading ATAC DAR summary table")

ATAC_summary <- read.delim(
  snakemake@input[["ATAC_summary"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "ATAC change",
  "FOXL2 overlap",
  "SOX9 overlap",
  "DMRT1 overlap",
  "Nearest expressed gene 50 kb",
  "Expression change 8w",
  "Gene sex bias"
)

missing_columns <- required_columns[
  !required_columns %in% colnames(ATAC_summary)
]

if (length(missing_columns) > 0) {
  log_error(
    "Missing columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

###########################################
# Functions
###########################################

run_fisher_test <- function(group, overlap) {
  valid <- !is.na(group) & !is.na(overlap)

  group <- group[valid]
  overlap <- overlap[valid]

  contingency_table <- matrix(
    c(
      sum(group == "b" & overlap),
      sum(group == "b" & !overlap),
      sum(group == "a" & overlap),
      sum(group == "a" & !overlap)
    ),
    nrow = 2,
    byrow = TRUE
  )

  rownames(contingency_table) <- c(
    "b",
    "a"
  )

  colnames(contingency_table) <- c(
    "Overlap",
    "No overlap"
  )

  fisher.test(
    contingency_table,
    alternative = "two.sided"
  )
}

make_gene_sex_bias_plot <- function(plot_data) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = ATAC_group,
      y = percent,
      fill = Gene_sex_bias
    )
  ) +
    ggplot2::geom_bar(
      width = 0.7,
      stat = "identity",
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = label
      ),
      position = ggplot2::position_stack(
        vjust = 0.5
      ),
      size = 4.2,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Granulosa" = "#b741b7",
        "Sertoli" = "#316695"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        "a" = "a\nLoss",
        "b" = "b\nGain"
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(
        0,
        100
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Sex-biased nearby genes"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 11
      ),
      axis.text = ggplot2::element_text(
        size = 10,
        colour = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        size = 10
      )
    )
}

make_expression_concordance_plot <- function(plot_data) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = ATAC_group,
      y = percent,
      fill = Expression_direction
    )
  ) +
    ggplot2::geom_bar(
      width = 0.7,
      stat = "identity",
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = label
      ),
      position = ggplot2::position_stack(
        vjust = 0.5
      ),
      size = 4.2,
      lineheight = 0.9
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Down" = "#b741b7",
        "Up" = "#316695"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = c(
        "a" = "a\nLoss",
        "b" = "b\nGain"
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(
        0,
        100
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Differentially expressed nearby genes"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 11
      ),
      axis.text = ggplot2::element_text(
        size = 10,
        colour = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        size = 10
      )
    )
}

make_TF_occupancy_plot <- function(plot_data) {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = TF,
      y = percent,
      fill = ATAC_group
    )
  ) +
    ggplot2::geom_bar(
      width = 0.7,
      stat = "identity",
      position = ggplot2::position_dodge(
        width = 0.75
      ),
      colour = "white",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = label
      ),
      position = ggplot2::position_dodge(
        width = 0.75
      ),
      vjust = -0.4,
      size = 3.8
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "a" = "#b741b7",
        "b" = "#316695"
      ),
      labels = c(
        "a" = "a - Loss",
        "b" = "b - Gain"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(
        0,
        105
      ),
      breaks = seq(
        0,
        100,
        20
      ),
      expand = c(
        0,
        0
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Regions overlapping TF ChIP-seq peaks"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = 11
      ),
      axis.text = ggplot2::element_text(
        size = 10,
        colour = "black"
      ),
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(
        size = 10
      )
    )
}

###########################################
# Define ATAC groups
###########################################

log_message("Defining ATAC groups a and b")

ATAC_summary$ATAC_group <- NA_character_

ATAC_summary$ATAC_group[
  ATAC_summary[["ATAC change"]] %in% c(
    "a1",
    "a2"
  )
] <- "a"

ATAC_summary$ATAC_group[
  ATAC_summary[["ATAC change"]] %in% c(
    "b1",
    "b2"
  )
] <- "b"

ATAC_summary$ATAC_group <- factor(
  ATAC_summary$ATAC_group,
  levels = c(
    "a",
    "b"
  )
)

log_message(
  "Group a regions: ",
  sum(
    ATAC_summary$ATAC_group == "a",
    na.rm = TRUE
  )
)

log_message(
  "Group b regions: ",
  sum(
    ATAC_summary$ATAC_group == "b",
    na.rm = TRUE
  )
)

###########################################
# Prepare nearby gene table
###########################################

log_message("Preparing nearby gene table")

valid_gene <- !is.na(
  ATAC_summary[["Nearest expressed gene 50 kb"]]
) &
  ATAC_summary[["Nearest expressed gene 50 kb"]] != "None" &
  !is.na(
    ATAC_summary$ATAC_group
  )

gene_group_table <- unique(
  ATAC_summary[
    valid_gene,
    c(
      "Nearest expressed gene 50 kb",
      "ATAC_group"
    )
  ]
)

gene_group_count <- table(
  gene_group_table[["Nearest expressed gene 50 kb"]]
)

genes_in_both_groups <- names(
  gene_group_count[
    gene_group_count > 1
  ]
)

log_message(
  "Genes associated with both groups a and b: ",
  length(
    genes_in_both_groups
  )
)

keep_gene <- valid_gene &
  !ATAC_summary[["Nearest expressed gene 50 kb"]] %in%
    genes_in_both_groups

gene_table <- ATAC_summary[
  keep_gene,
  c(
    "Nearest expressed gene 50 kb",
    "ATAC_group",
    "Expression change 8w",
    "Gene sex bias"
  )
]

gene_id <- paste(
  gene_table[["Nearest expressed gene 50 kb"]],
  gene_table$ATAC_group,
  sep = "_"
)

gene_table <- gene_table[
  !duplicated(
    gene_id
  ),
  ,
  drop = FALSE
]

###########################################
# Nearby gene sex-bias analysis
###########################################

log_message("Calculating nearby gene sex-bias statistics")

sex_bias_gene_table <- gene_table[
  gene_table[["Gene sex bias"]] %in% c(
    "Granulosa",
    "Sertoli"
  ),
  ,
  drop = FALSE
]

log_message(
  "Sex-biased genes retained for analysis: ",
  nrow(
    sex_bias_gene_table
  )
)

gene_count_table <- table(
  sex_bias_gene_table$ATAC_group,
  sex_bias_gene_table[["Gene sex bias"]]
)

gene_count_table <- gene_count_table[
  c(
    "a",
    "b"
  ),
  c(
    "Granulosa",
    "Sertoli"
  ),
  drop = FALSE
]

gene_fisher_table <- matrix(
  c(
    gene_count_table[
      "b",
      "Sertoli"
    ],
    gene_count_table[
      "b",
      "Granulosa"
    ],
    gene_count_table[
      "a",
      "Sertoli"
    ],
    gene_count_table[
      "a",
      "Granulosa"
    ]
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(gene_fisher_table) <- c(
  "b",
  "a"
)

colnames(gene_fisher_table) <- c(
  "Sertoli",
  "Granulosa"
)

gene_fisher_test <- fisher.test(
  gene_fisher_table
)

gene_statistics <- data.frame(
  comparison = "Group b vs a: Sertoli vs Granulosa nearby genes",
  group_b_Sertoli = gene_fisher_table[
    "b",
    "Sertoli"
  ],
  group_b_Granulosa = gene_fisher_table[
    "b",
    "Granulosa"
  ],
  group_a_Sertoli = gene_fisher_table[
    "a",
    "Sertoli"
  ],
  group_a_Granulosa = gene_fisher_table[
    "a",
    "Granulosa"
  ],
  odds_ratio = unname(
    gene_fisher_test$estimate
  ),
  confidence_low = gene_fisher_test$conf.int[1],
  confidence_high = gene_fisher_test$conf.int[2],
  p_value = gene_fisher_test$p.value,
  stringsAsFactors = FALSE
)

###########################################
# Prepare nearby gene sex-bias plot
###########################################

gene_proportion_table <- prop.table(
  gene_count_table,
  margin = 1
)

gene_plot_data <- data.frame(
  ATAC_group = rep(
    rownames(
      gene_count_table
    ),
    each = ncol(
      gene_count_table
    )
  ),
  Gene_sex_bias = rep(
    colnames(
      gene_count_table
    ),
    times = nrow(
      gene_count_table
    )
  ),
  count = as.vector(
    t(
      gene_count_table
    )
  ),
  percent = 100 * as.vector(
    t(
      gene_proportion_table
    )
  ),
  stringsAsFactors = FALSE
)

gene_plot_data$ATAC_group <- factor(
  gene_plot_data$ATAC_group,
  levels = c(
    "a",
    "b"
  )
)

gene_plot_data$Gene_sex_bias <- factor(
  gene_plot_data$Gene_sex_bias,
  levels = c(
    "Granulosa",
    "Sertoli"
  )
)

gene_plot_data$label <- paste0(
  round(
    gene_plot_data$percent,
    1
  ),
  "%\n(",
  scales::comma(
    gene_plot_data$count
  ),
  ")"
)

plot_gene_sex_bias <- make_gene_sex_bias_plot(
  gene_plot_data
)

###########################################
# Nearby gene expression-direction analysis
###########################################

log_message(
  "Calculating concordance between ATAC and gene expression changes"
)

expression_gene_table <- gene_table[
  gene_table[["Expression change 8w"]] %in% c(
    "Down",
    "Up"
  ),
  ,
  drop = FALSE
]

log_message(
  "Differentially expressed nearby genes retained: ",
  nrow(
    expression_gene_table
  )
)

expression_count_table <- table(
  expression_gene_table$ATAC_group,
  expression_gene_table[["Expression change 8w"]]
)

expression_count_table <- expression_count_table[
  c(
    "a",
    "b"
  ),
  c(
    "Down",
    "Up"
  ),
  drop = FALSE
]

expression_fisher_table <- matrix(
  c(
    expression_count_table[
      "b",
      "Up"
    ],
    expression_count_table[
      "b",
      "Down"
    ],
    expression_count_table[
      "a",
      "Up"
    ],
    expression_count_table[
      "a",
      "Down"
    ]
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(expression_fisher_table) <- c(
  "b",
  "a"
)

colnames(expression_fisher_table) <- c(
  "Up",
  "Down"
)

expression_fisher_test <- fisher.test(
  expression_fisher_table
)

expression_statistics <- data.frame(
  comparison = "Group b vs a: Up vs Down nearby genes",
  group_b_Up = expression_fisher_table[
    "b",
    "Up"
  ],
  group_b_Down = expression_fisher_table[
    "b",
    "Down"
  ],
  group_a_Up = expression_fisher_table[
    "a",
    "Up"
  ],
  group_a_Down = expression_fisher_table[
    "a",
    "Down"
  ],
  odds_ratio = unname(
    expression_fisher_test$estimate
  ),
  confidence_low = expression_fisher_test$conf.int[1],
  confidence_high = expression_fisher_test$conf.int[2],
  p_value = expression_fisher_test$p.value,
  stringsAsFactors = FALSE
)

log_message(
  "Expression-direction Fisher test: OR = ",
  round(
    expression_statistics$odds_ratio,
    2
  ),
  "; P = ",
  format(
    expression_statistics$p_value,
    scientific = TRUE,
    digits = 3
  )
)

###########################################
# Prepare expression-direction plot
###########################################

expression_proportion_table <- prop.table(
  expression_count_table,
  margin = 1
)

expression_plot_data <- data.frame(
  ATAC_group = rep(
    rownames(
      expression_count_table
    ),
    each = ncol(
      expression_count_table
    )
  ),
  Expression_direction = rep(
    colnames(
      expression_count_table
    ),
    times = nrow(
      expression_count_table
    )
  ),
  count = as.vector(
    t(
      expression_count_table
    )
  ),
  percent = 100 * as.vector(
    t(
      expression_proportion_table
    )
  ),
  stringsAsFactors = FALSE
)

expression_plot_data$ATAC_group <- factor(
  expression_plot_data$ATAC_group,
  levels = c(
    "a",
    "b"
  )
)

expression_plot_data$Expression_direction <- factor(
  expression_plot_data$Expression_direction,
  levels = c(
    "Down",
    "Up"
  )
)

expression_plot_data$label <- paste0(
  round(
    expression_plot_data$percent,
    1
  ),
  "%\n(",
  scales::comma(
    expression_plot_data$count
  ),
  ")"
)

plot_expression_concordance <- make_expression_concordance_plot(
  expression_plot_data
)

###########################################
# Prepare TF occupancy analysis
###########################################

log_message("Preparing TF occupancy analysis")

TF_names <- c(
  "FOXL2",
  "SOX9",
  "DMRT1"
)

TF_overlap_columns <- c(
  "FOXL2 overlap",
  "SOX9 overlap",
  "DMRT1 overlap"
)

TF_statistics <- data.frame(
  TF = character(),
  group_a_total = integer(),
  group_a_overlap = integer(),
  group_a_percent = numeric(),
  group_b_total = integer(),
  group_b_overlap = integer(),
  group_b_percent = numeric(),
  odds_ratio_b_vs_a = numeric(),
  confidence_low = numeric(),
  confidence_high = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

TF_plot_data <- data.frame(
  TF = character(),
  ATAC_group = character(),
  count = integer(),
  total = integer(),
  percent = numeric(),
  stringsAsFactors = FALSE
)

###########################################
# Calculate TF occupancy statistics
###########################################

for (i in seq_along(
  TF_names
)) {
  TF <- TF_names[i]
  overlap_column <- TF_overlap_columns[i]

  valid_region <- !is.na(
    ATAC_summary$ATAC_group
  ) &
    !is.na(
      ATAC_summary[[overlap_column]]
    )

  current_data <- ATAC_summary[
    valid_region,
    ,
    drop = FALSE
  ]

  overlap <- current_data[[overlap_column]] ==
    "Yes"

  group_a <- current_data$ATAC_group ==
    "a"

  group_b <- current_data$ATAC_group ==
    "b"

  group_a_total <- sum(
    group_a
  )

  group_b_total <- sum(
    group_b
  )

  group_a_overlap <- sum(
    group_a & overlap
  )

  group_b_overlap <- sum(
    group_b & overlap
  )

  if (group_a_total > 0) {
    group_a_percent <- 100 *
      group_a_overlap /
      group_a_total
  } else {
    group_a_percent <- NA_real_
  }

  if (group_b_total > 0) {
    group_b_percent <- 100 *
      group_b_overlap /
      group_b_total
  } else {
    group_b_percent <- NA_real_
  }

  fisher_result <- run_fisher_test(
    as.character(
      current_data$ATAC_group
    ),
    overlap
  )

  current_statistics <- data.frame(
    TF = TF,
    group_a_total = group_a_total,
    group_a_overlap = group_a_overlap,
    group_a_percent = group_a_percent,
    group_b_total = group_b_total,
    group_b_overlap = group_b_overlap,
    group_b_percent = group_b_percent,
    odds_ratio_b_vs_a = unname(
      fisher_result$estimate
    ),
    confidence_low = fisher_result$conf.int[1],
    confidence_high = fisher_result$conf.int[2],
    p_value = fisher_result$p.value,
    stringsAsFactors = FALSE
  )

  TF_statistics <- rbind(
    TF_statistics,
    current_statistics
  )

  current_plot_data <- data.frame(
    TF = c(
      TF,
      TF
    ),
    ATAC_group = c(
      "a",
      "b"
    ),
    count = c(
      group_a_overlap,
      group_b_overlap
    ),
    total = c(
      group_a_total,
      group_b_total
    ),
    percent = c(
      group_a_percent,
      group_b_percent
    ),
    stringsAsFactors = FALSE
  )

  TF_plot_data <- rbind(
    TF_plot_data,
    current_plot_data
  )
}

###########################################
# Log TF occupancy statistics
###########################################

for (i in seq_len(
  nrow(
    TF_statistics
  )
)) {
  log_message(
    TF_statistics$TF[i],
    ": group a = ",
    TF_statistics$group_a_overlap[i],
    "/",
    TF_statistics$group_a_total[i],
    " (",
    round(
      TF_statistics$group_a_percent[i],
      1
    ),
    "%); group b = ",
    TF_statistics$group_b_overlap[i],
    "/",
    TF_statistics$group_b_total[i],
    " (",
    round(
      TF_statistics$group_b_percent[i],
      1
    ),
    "%); OR b/a = ",
    round(
      TF_statistics$odds_ratio_b_vs_a[i],
      2
    ),
    " [",
    round(
      TF_statistics$confidence_low[i],
      2
    ),
    "-",
    round(
      TF_statistics$confidence_high[i],
      2
    ),
    "]; P = ",
    format(
      TF_statistics$p_value[i],
      scientific = TRUE,
      digits = 3
    )
  )
}

###########################################
# Prepare TF occupancy plot
###########################################

TF_plot_data$TF <- factor(
  TF_plot_data$TF,
  levels = c(
    "FOXL2",
    "SOX9",
    "DMRT1"
  )
)

TF_plot_data$ATAC_group <- factor(
  TF_plot_data$ATAC_group,
  levels = c(
    "a",
    "b"
  )
)

TF_plot_data$label <- paste0(
  round(
    TF_plot_data$percent,
    1
  ),
  "%\n(",
  scales::comma(
    TF_plot_data$count
  ),
  ")"
)

plot_TF_occupancy <- make_TF_occupancy_plot(
  TF_plot_data
)

###########################################
# Write statistics
###########################################

log_message("Writing statistics")

write.table(
  gene_statistics,
  file = snakemake@output[["gene_statistics"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  expression_statistics,
  file = snakemake@output[["expression_statistics"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  TF_statistics,
  file = snakemake@output[["TF_statistics"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Draw plots
###########################################

top_row <- cowplot::plot_grid(
  plot_gene_sex_bias,
  plot_expression_concordance,
  labels = c(
    "d1",
    "d2"
  ),
  ncol = 2,
  align = "h",
  axis = "tb"
)

bottom_row <- cowplot::plot_grid(
  plot_TF_occupancy,
  labels = "d3",
  ncol = 1
)

figure <- cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(
    1,
    1
  )
)

###########################################
# Save output files
###########################################

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 16,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 16,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")