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

log_message("Starting TE enrichment analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("scales")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

AB <- snakemake@params[["AB"]]

###########################################
# Load input data
###########################################

log_message("Importing WT peaks")

WT_peaks <- rtracklayer::import(
  snakemake@input[["WT_peaks"]]
)

log_message("Loading differential regions")

DER <- read.csv(
  snakemake@input[["DER"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"x" %in% colnames(DER)) {
  log_error("Column 'x' not found in DER table")
}

log_message("Loading RepeatMasker annotation")

repeat_masker <- read.csv(
  snakemake@input[["repeatMasker"]],
  sep = "\t",
  header = TRUE
)

rmsk <- GenomicRanges::makeGRangesFromDataFrame(
  repeat_masker,
  keep.extra.columns = TRUE
)

if (!all(c("repFamily", "repClass") %in% colnames(S4Vectors::mcols(rmsk)))) {
  log_error("RepeatMasker annotation must contain 'repFamily' and 'repClass' columns")
}

###########################################
# Keep only transposable elements
###########################################

log_message("Keeping only transposable elements")

te_keep <- !is.na(rmsk$repClass) &
  rmsk$repClass != "" &
  !grepl(
    "Simple_repeat|Low_complexity|Satellite|RNA|rRNA|tRNA|snRNA|scRNA|srpRNA",
    rmsk$repClass,
    ignore.case = TRUE
  )

rmsk <- rmsk[te_keep]

log_message("Number of TE annotations retained: ", length(rmsk))

###########################################
# Functions
###########################################

get_repeats_overlapping_regions <- function(
  regions,
  repeats
) {
  hits <- GenomicRanges::findOverlaps(
    repeats,
    regions,
    ignore.strand = TRUE
  )

  unique(
    S4Vectors::queryHits(hits)
  )
}

make_te_class_mapping <- function(
  repeats
) {
  repeat_data <- data.frame(
    repFamily = as.character(repeats$repFamily),
    repClass = as.character(repeats$repClass),
    stringsAsFactors = FALSE
  )

  repeat_data <- repeat_data[
    !is.na(repeat_data$repFamily) &
      repeat_data$repFamily != "",
    ,
    drop = FALSE
  ]

  repeat_data <- repeat_data[
    !duplicated(repeat_data$repFamily),
    ,
    drop = FALSE
  ]

  te_class <- repeat_data$repClass
  names(te_class) <- repeat_data$repFamily

  te_class
}

run_fisher_by_family <- function(
  target_repeat_idx,
  background_repeat_idx,
  rep_family
) {
  families <- sort(
    unique(as.character(rep_family[background_repeat_idx]))
  )

  background_without_target <- setdiff(
    background_repeat_idx,
    target_repeat_idx
  )

  results <- lapply(families, function(family) {
    family_idx <- which(
      as.character(rep_family) == family
    )

    a <- sum(family_idx %in% target_repeat_idx)
    b <- length(target_repeat_idx) - a
    c <- sum(family_idx %in% background_without_target)
    d <- length(background_without_target) - c

    fisher_matrix <- matrix(
      c(a, b, c, d),
      nrow = 2
    )

    fisher_result <- tryCatch(
      stats::fisher.test(
        fisher_matrix,
        alternative = "two.sided"
      ),
      error = function(e) NULL
    )

    if (is.null(fisher_result)) {
      odds <- NA_real_
      pval <- NA_real_
    } else {
      odds <- unname(fisher_result$estimate)
      pval <- fisher_result$p.value
    }

    data.frame(
      repFamily = family,
      target_TE = a,
      target_other_TE = b,
      background_TE = c,
      background_other_TE = d,
      odds = odds,
      pval = pval,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(
    rbind,
    results
  )

  results$padj <- stats::p.adjust(
    results$pval,
    method = "BH"
  )

  results$log2OR <- log2(
    results$odds
  )

  results$negLog10P <- -log10(
    results$pval
  )

  results$direction <- ifelse(
    results$log2OR >= 0,
    "Enriched",
    "Depleted"
  )

  results
}

prepare_te_enrichment <- function(
  regions,
  repeats,
  label
) {
  log_message("Computing TE enrichment for: ", label)

  target_repeat_idx <- get_repeats_overlapping_regions(
    regions = regions,
    repeats = repeats
  )

  background_repeat_idx <- seq_along(
    repeats
  )

  log_message("Number of overlapping TE elements: ", length(target_repeat_idx))

  enrichment <- run_fisher_by_family(
    target_repeat_idx = target_repeat_idx,
    background_repeat_idx = background_repeat_idx,
    rep_family = repeats$repFamily
  )

  enrichment$comparison <- label

  enrichment
}

plot_te_enrichment <- function(
  enrichment,
  te_class,
  title,
  pval_cutoff = 0.05
) {
  plot_data <- enrichment

  plot_data$class <- te_class[
    match(
      plot_data$repFamily,
      names(te_class)
    )
  ]

  plot_data$class[is.na(plot_data$class)] <- "Unknown"
  plot_data$significant <- plot_data$padj <= pval_cutoff

  plot_data <- plot_data[
    is.finite(plot_data$log2OR) &
      is.finite(plot_data$negLog10P),
    ,
    drop = FALSE
  ]

  label_data <- plot_data[
    plot_data$significant,
    ,
    drop = FALSE
  ]

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = log2OR,
      y = negLog10P,
      colour = direction
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 0,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = "#FCEBF0",
      alpha = 0.7
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = 0,
      ymin = -Inf,
      ymax = Inf,
      fill = "#EAF3FA",
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(pval_cutoff),
      linetype = "dashed",
      colour = "#666666",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "#666666",
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = Inf,
      label = "Enriched\nmore than expected",
      hjust = 1.05,
      vjust = 1.25,
      colour = "#d23a5b",
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = "Depleted\nless than expected",
      hjust = -0.05,
      vjust = 1.25,
      colour = "#4A90C2",
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = target_TE,
        alpha = significant
      ),
      stroke = 0
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(
        label = repFamily
      ),
      size = 4.3,
      colour = "#222222",
      max.overlaps = 20,
      box.padding = 0.4,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.colour = "#999999"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Enriched" = "#d23a5b",
        "Depleted" = "#4A90C2"
      )
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "TRUE" = 0.95,
        "FALSE" = 0.35
      ),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 7),
      labels = scales::comma
    ) +
    ggplot2::labs(
      title = title,
      x = "log2(odds ratio)",
      y = "-log10(p-value)",
      colour = "Direction",
      size = "Overlapping TE elements"
    ) +
    ggplot2::theme_classic(
      base_size = 15
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 17
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      axis.text = ggplot2::element_text(
        size = 13,
        colour = "#333333"
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      legend.text = ggplot2::element_text(
        size = 12
      ),
      legend.position = "right",
      panel.border = ggplot2::element_rect(
        colour = "#333333",
        fill = NA,
        linewidth = 0.4
      )
    )
}

write_enrichment_table <- function(
  enrichment,
  output_file
) {
  write.table(
    enrichment,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

###########################################
# Prepare regions
###########################################

log_message("Preparing H3K9me3 differential region groups")

peaks_lost <- GenomicRanges::GRanges(
  rownames(
    DER[
      DER$x == "a",
      ,
      drop = FALSE
    ]
  )
)

peaks_gained <- GenomicRanges::GRanges(
  rownames(
    DER[
      DER$x == "b",
      ,
      drop = FALSE
    ]
  )
)

log_message("Number of regions with H3K9me3 loss: ", length(peaks_lost))
log_message("Number of regions with H3K9me3 gain: ", length(peaks_gained))

###########################################
# Run TE enrichment
###########################################

te_class <- make_te_class_mapping(
  repeats = rmsk
)

res_WT <- prepare_te_enrichment(
  regions = WT_peaks,
  repeats = rmsk,
  label = paste0(AB, " WT H3K9me3 peaks")
)

res_lost <- prepare_te_enrichment(
  regions = peaks_lost,
  repeats = rmsk,
  label = paste0(AB, " H3K9me3 loss")
)

res_gained <- prepare_te_enrichment(
  regions = peaks_gained,
  repeats = rmsk,
  label = paste0(AB, " H3K9me3 gain")
)

###########################################
# Save result tables
###########################################

log_message("Writing TE enrichment result tables")

write_enrichment_table(
  res_WT,
  snakemake@output[["WT_table"]]
)

write_enrichment_table(
  res_lost,
  snakemake@output[["lost_table"]]
)

write_enrichment_table(
  res_gained,
  snakemake@output[["gained_table"]]
)

###########################################
# Plot results
###########################################

log_message("Plotting TE enrichment results")

p1 <- plot_te_enrichment(
  enrichment = res_WT,
  te_class = te_class,
  title = paste("TE families enriched/depleted in WT H3K9me3 peaks")
)

p2 <- plot_te_enrichment(
  enrichment = res_lost,
  te_class = te_class,
  title = paste("TE families enriched/depleted in regions with H3K9me3 loss")
)

p3 <- plot_te_enrichment(
  enrichment = res_gained,
  te_class = te_class,
  title = paste("TE families enriched/depleted in regions with H3K9me3 gain")
)

figure <- cowplot::plot_grid(
  plotlist = list(
    p1,
    p2,
    p3
  ),
  labels = "AUTO",
  ncol = 1,
  align = "v"
)

###########################################
# Save figures
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 30,
  base_height = 34,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 30,
  base_height = 34,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")