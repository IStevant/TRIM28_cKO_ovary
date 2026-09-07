source(".Rprofile")

suppressPackageStartupMessages({
  library("cowplot")
  library("EnrichedHeatmap")
  library("GenomicRanges")
  library("ggplot2")
  library("IRanges")
  library("rtracklayer")
  library("scales")
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

# Fixed universe: WT FOXL2 peaks
FOXL2_WT_peak_file <- snakemake@input[["FOXL2_WT"]]

# Normalised FOXL2 counts on the WT FOXL2 peak universe
FOXL2_count_file <- snakemake@input[["FOXL2_counts"]]

# Other transcription factor peaks
TRIM28_peak_file <- snakemake@input[["TRIM28"]]
NR5A2_peak_file <- snakemake@input[["NR5A2"]]
ESR2_peak_file <- snakemake@input[["ESR2"]]
RUNX_peak_file <- snakemake@input[["RUNX"]]

# BigWigs
FOXL2_WT_bw <- snakemake@input[["FOXL2_WT_bw"]]
FOXL2_cKO_bw <- snakemake@input[["FOXL2_cKO_bw"]]

ATAC_WT_bw <- snakemake@input[["ATAC_WT_bw"]]
ATAC_cKO_bw <- snakemake@input[["ATAC_cKO_bw"]]

SUMO1_WT_bw <- snakemake@input[["SUMO1_WT_bw"]]
SUMO1_cKO_bw <- snakemake@input[["SUMO1_cKO_bw"]]

SUMO2_WT_bw <- snakemake@input[["SUMO2_WT_bw"]]
SUMO2_cKO_bw <- snakemake@input[["SUMO2_cKO_bw"]]

# Parameters
WT_column <- snakemake@params[["WT_column"]]
cKO_column <- snakemake@params[["cKO_column"]]

pseudocount <- if ("pseudocount" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["pseudocount"]])
} else {
  1
}

profile_window <- if ("profile_window" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["profile_window"]])
} else {
  1000
}

bin_size <- if ("bin_size" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["bin_size"]])
} else {
  50
}

loss_fraction <- if ("loss_fraction" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["loss_fraction"]])
} else {
  0.25
}

# Outputs
output_regions <- snakemake@output[["regions"]]
output_summary <- snakemake@output[["summary"]]
output_stats <- snakemake@output[["stats"]]
output_enrichment <- if ("enrichment" %in% names(snakemake@output)) {
  snakemake@output[["enrichment"]]
} else {
  sub("\\.tsv$", "_enrichment.tsv", output_stats)
}
output_profiles <- snakemake@output[["profiles"]]

output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]


###########################################
# Load genomic regions
###########################################

log_message("Loading FOXL2 WT peaks.")

FOXL2_WT <- rtracklayer::import(
  FOXL2_WT_peak_file
)

if (length(FOXL2_WT) == 0) {
  log_error("No FOXL2 WT peaks were loaded.")
}

log_message(
  "FOXL2 WT peak universe: ",
  length(FOXL2_WT),
  " regions."
)

TRIM28 <- rtracklayer::import(
  TRIM28_peak_file
)

TF_list <- list(
  NR5A2 = rtracklayer::import(NR5A2_peak_file),
  ESR2 = rtracklayer::import(ESR2_peak_file),
  RUNX = rtracklayer::import(RUNX_peak_file)
)


###########################################
# Functions
###########################################

make_region_id <- function(gr) {

  paste0(
    as.character(seqnames(gr)),
    ":",
    start(gr),
    "-",
    end(gr)
  )
}


count_overlap_binary <- function(query, subject) {

  GenomicRanges::countOverlaps(
    query,
    subject,
    ignore.strand = TRUE
  ) > 0
}


# Associate normalised FOXL2 counts with the fixed WT FOXL2 peak universe
# using genomic overlaps rather than exact coordinate matching.
prepare_FOXL2_counts <- function(
    FOXL2_WT,
    count_file,
    WT_column,
    cKO_column,
    pseudocount) {

  log_message("Loading normalised FOXL2 counts.")

  counts <- read.csv(
    count_file,
    row.names = 1,
    check.names = FALSE
  )

  if (!WT_column %in% colnames(counts)) {
    log_error(
      "WT FOXL2 column not found in count matrix: ",
      WT_column
    )
  }

  if (!cKO_column %in% colnames(counts)) {
    log_error(
      "cKO FOXL2 column not found in count matrix: ",
      cKO_column
    )
  }

  log_message(
    "Count matrix contains ",
    nrow(counts),
    " quantified regions."
  )

  count_regions <- tryCatch(
    {
      GenomicRanges::GRanges(
        rownames(counts)
      )
    },
    error = function(e) {
      log_error(
        "Could not convert count-matrix row names to GRanges. ",
        "Expected region names such as chr1:100-200."
      )
    }
  )

  hits <- GenomicRanges::findOverlaps(
    count_regions,
    FOXL2_WT,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    log_error(
      "No quantified regions overlap the FOXL2 WT peak universe."
    )
  }

  count_indices <- sort(
    unique(
      S4Vectors::queryHits(hits)
    )
  )

  selected_regions <- count_regions[
    count_indices
  ]

  selected_counts <- counts[
    count_indices,
    ,
    drop = FALSE
  ]

  log_message(
    "Quantified regions overlapping FOXL2 WT peaks: ",
    length(selected_regions),
    " (",
    round(
      100 * length(selected_regions) / nrow(counts),
      1
    ),
    "% of quantified regions)."
  )

  WT_peak_hits <- GenomicRanges::findOverlaps(
    FOXL2_WT,
    selected_regions,
    ignore.strand = TRUE
  )

  log_message(
    "FOXL2 WT peaks represented by the selected quantified regions: ",
    length(
      unique(
        S4Vectors::queryHits(WT_peak_hits)
      )
    ),
    " / ",
    length(FOXL2_WT)
  )

  output <- data.frame(
    chromosome = as.character(
      seqnames(selected_regions)
    ),
    start = start(selected_regions),
    end = end(selected_regions),
    region = make_region_id(
      selected_regions
    ),
    FOXL2_WT = selected_counts[[WT_column]],
    FOXL2_cKO = selected_counts[[cKO_column]],
    stringsAsFactors = FALSE
  )

  output$delta_FOXL2 <- log2(
    (output$FOXL2_cKO + pseudocount) /
      (output$FOXL2_WT + pseudocount)
  )

  return(
    list(
      regions = selected_regions,
      table = output
    )
  )
}


# Annotate each quantified region associated with a FOXL2 WT peak directly
# for TRIM28 and ovarian TF occupancy.
annotate_FOXL2_regions <- function(
    FOXL2_regions,
    FOXL2_table,
    TRIM28,
    TF_list) {

  log_message(
    "Annotating quantified FOXL2 WT regions directly with TRIM28 and ovarian TF occupancy."
  )

  FOXL2_table$TRIM28 <- count_overlap_binary(
    FOXL2_regions,
    TRIM28
  )

  for (TF_name in names(TF_list)) {
    FOXL2_table[[TF_name]] <- count_overlap_binary(
      FOXL2_regions,
      TF_list[[TF_name]]
    )
  }

  TF_matrix <- as.matrix(
    FOXL2_table[
      ,
      names(TF_list),
      drop = FALSE
    ]
  )

  FOXL2_table$n_additional_TFs <- rowSums(
    TF_matrix
  )

  FOXL2_table$TRIM28_class <- ifelse(
    FOXL2_table$TRIM28,
    "FOXL2 + TRIM28",
    "FOXL2 without TRIM28"
  )

  FOXL2_table$hub_class <- NA_character_

  FOXL2_table$hub_class[
    FOXL2_table$TRIM28
  ] <- paste0(
    "TRIM28 + FOXL2 + ",
    FOXL2_table$n_additional_TFs[
      FOXL2_table$TRIM28
    ],
    ifelse(
      FOXL2_table$n_additional_TFs[
        FOXL2_table$TRIM28
      ] == 1,
      " additional TF",
      " additional TFs"
    )
  )

  FOXL2_table$hub_class <- factor(
    FOXL2_table$hub_class,
    levels = c(
      "TRIM28 + FOXL2 + 0 additional TFs",
      "TRIM28 + FOXL2 + 1 additional TF",
      "TRIM28 + FOXL2 + 2 additional TFs",
      "TRIM28 + FOXL2 + 3 additional TFs"
    )
  )

  log_message(
    "Quantified FOXL2 WT regions overlapping TRIM28: ",
    sum(FOXL2_table$TRIM28)
  )

  log_message("TRIM28-positive FOXL2 WT regions by additional TF number:")

  print(
    table(
      FOXL2_table$n_additional_TFs[
        FOXL2_table$TRIM28
      ]
    )
  )

  return(FOXL2_table)
}


# Define the strongest FOXL2-loss regions from the continuous delta.
define_FOXL2_loss <- function(
    FOXL2_table,
    loss_fraction = 0.25) {

  if (
    length(loss_fraction) != 1 ||
      is.na(loss_fraction) ||
      loss_fraction <= 0 ||
      loss_fraction >= 1
  ) {
    log_error(
      "loss_fraction must be a single value between 0 and 1."
    )
  }

  loss_threshold <- stats::quantile(
    FOXL2_table$delta_FOXL2,
    probs = loss_fraction,
    na.rm = TRUE,
    names = FALSE
  )

  FOXL2_table$FOXL2_loss <- (
    FOXL2_table$delta_FOXL2 <= loss_threshold
  )

  FOXL2_table$loss_class <- ifelse(
    FOXL2_table$FOXL2_loss,
    "Top FOXL2 loss",
    "Remaining FOXL2 sites"
  )

  log_message(
    "FOXL2 loss threshold (bottom ",
    loss_fraction * 100,
    "%): delta_FOXL2 <= ",
    signif(loss_threshold, 4)
  )

  log_message(
    "Strong FOXL2-loss regions: ",
    sum(FOXL2_table$FOXL2_loss),
    " / ",
    nrow(FOXL2_table)
  )

  return(FOXL2_table)
}


# Test enrichment of regulatory features among the strongest FOXL2-loss sites.
run_loss_enrichment <- function(
    FOXL2_table) {

  log_message(
    "Testing enrichment of TRIM28 and ovarian TF occupancy among strongest FOXL2-loss regions."
  )

  feature_table <- data.frame(
    feature = c(
      "TRIM28",
      "NR5A2",
      "ESR2",
      "RUNX",
      "TRIM28 + >=1 TF",
      "TRIM28 + >=2 TFs",
      "TRIM28 + 3 TFs"
    ),
    stringsAsFactors = FALSE
  )

  feature_vectors <- list(
    "TRIM28" = FOXL2_table$TRIM28,
    "NR5A2" = FOXL2_table$NR5A2,
    "ESR2" = FOXL2_table$ESR2,
    "RUNX" = FOXL2_table$RUNX,
    "TRIM28 + >=1 TF" = (
      FOXL2_table$TRIM28 &
        FOXL2_table$n_additional_TFs >= 1
    ),
    "TRIM28 + >=2 TFs" = (
      FOXL2_table$TRIM28 &
        FOXL2_table$n_additional_TFs >= 2
    ),
    "TRIM28 + 3 TFs" = (
      FOXL2_table$TRIM28 &
        FOXL2_table$n_additional_TFs == 3
    )
  )

  results <- lapply(
    feature_table$feature,
    function(feature_name) {

      feature_present <- feature_vectors[[feature_name]]
      loss <- FOXL2_table$FOXL2_loss

      contingency <- matrix(
        c(
          sum(loss & feature_present),
          sum(loss & !feature_present),
          sum(!loss & feature_present),
          sum(!loss & !feature_present)
        ),
        nrow = 2,
        byrow = TRUE
      )

      fisher_result <- fisher.test(
        contingency,
        alternative = "greater"
      )

      n_loss <- sum(loss)
      n_remaining <- sum(!loss)

      loss_with_feature <- sum(
        loss & feature_present
      )

      remaining_with_feature <- sum(
        !loss & feature_present
      )

      data.frame(
        feature = feature_name,
        n_loss = n_loss,
        loss_with_feature = loss_with_feature,
        loss_percent = 100 * loss_with_feature / n_loss,
        n_remaining = n_remaining,
        remaining_with_feature = remaining_with_feature,
        remaining_percent = 100 * remaining_with_feature / n_remaining,
        odds_ratio = unname(fisher_result$estimate),
        conf_low = fisher_result$conf.int[1],
        conf_high = fisher_result$conf.int[2],
        p_value = fisher_result$p.value,
        stringsAsFactors = FALSE
      )
    }
  )

  results <- do.call(
    rbind,
    results
  )

  results$FDR <- p.adjust(
    results$p_value,
    method = "BH"
  )

  return(results)
}


# Continuous region-level tests of FOXL2 loss according to regulatory context.
run_continuous_statistics <- function(
    FOXL2_table) {

  log_message(
    "Running continuous region-level FOXL2 occupancy tests."
  )

  trim28_test <- wilcox.test(
    delta_FOXL2 ~ TRIM28,
    data = FOXL2_table,
    alternative = "less",
    exact = FALSE
  )

  hub_data <- FOXL2_table[
    FOXL2_table$TRIM28,
    ,
    drop = FALSE
  ]

  spearman_test <- cor.test(
    hub_data$n_additional_TFs,
    hub_data$delta_FOXL2,
    method = "spearman",
    alternative = "less",
    exact = FALSE
  )

  data.frame(
    test = c(
      "Delta FOXL2: TRIM28+ vs TRIM28-",
      "Delta FOXL2 vs number of additional TFs at TRIM28+ sites"
    ),
    statistic = c(
      unname(trim28_test$statistic),
      unname(spearman_test$estimate)
    ),
    statistic_name = c(
      "Wilcoxon W",
      "Spearman rho"
    ),
    p_value = c(
      trim28_test$p.value,
      spearman_test$p.value
    ),
    stringsAsFactors = FALSE
  )
}


summarise_delta <- function(FOXL2_table) {

  log_message("Summarising relative FOXL2 occupancy.")

  trim28_summary <- do.call(
    rbind,
    lapply(
      split(
        FOXL2_table,
        FOXL2_table$TRIM28_class
      ),
      function(data) {

        data.frame(
          comparison = "TRIM28 occupancy",
          class = unique(data$TRIM28_class),
          n_regions = nrow(data),
          median_delta = median(
            data$delta_FOXL2,
            na.rm = TRUE
          ),
          mean_delta = mean(
            data$delta_FOXL2,
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  hub_data <- FOXL2_table[
    FOXL2_table$TRIM28,
    ,
    drop = FALSE
  ]

  hub_summary <- do.call(
    rbind,
    lapply(
      split(
        hub_data,
        hub_data$n_additional_TFs
      ),
      function(data) {

        data.frame(
          comparison = "Hub complexity",
          class = paste0(
            unique(data$n_additional_TFs),
            " additional TFs"
          ),
          n_regions = nrow(data),
          median_delta = median(
            data$delta_FOXL2,
            na.rm = TRUE
          ),
          mean_delta = mean(
            data$delta_FOXL2,
            na.rm = TRUE
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rbind(
    trim28_summary,
    hub_summary
  )
}




###########################################
# Signal profiles
###########################################

extract_profile_matrix <- function(
    bigwig_file,
    regions,
    extend,
    bin_size) {

  log_message("Loading BigWig signal: ", bigwig_file)

  signal <- rtracklayer::import(
    bigwig_file
  )

  if (length(signal) == 0) {
    log_error(
      "No signal was imported from BigWig file: ",
      bigwig_file
    )
  }

  if (!"score" %in% colnames(mcols(signal))) {
    log_error(
      "The imported BigWig does not contain a 'score' column: ",
      bigwig_file
    )
  }

  centred_regions <- GenomicRanges::resize(
    regions,
    width = 1,
    fix = "center"
  )

  common_seqlevels <- intersect(
    GenomeInfoDb::seqlevels(centred_regions),
    GenomeInfoDb::seqlevels(signal)
  )

  if (length(common_seqlevels) == 0) {
    log_error(
      "No common chromosomes were found between regions and BigWig: ",
      bigwig_file
    )
  }

  keep_regions <- as.character(
    GenomicRanges::seqnames(centred_regions)
  ) %in% common_seqlevels

  if (!all(keep_regions)) {
    log_message(
      sum(!keep_regions),
      " region(s) are on chromosomes absent from ",
      basename(bigwig_file),
      "; their profile values will be set to NA."
    )
  }

  signal_matrix <- matrix(
    NA_real_,
    nrow = length(centred_regions),
    ncol = ceiling((2 * extend) / bin_size)
  )

  if (any(keep_regions)) {

    matrix_common <- EnrichedHeatmap::normalizeToMatrix(
      signal,
      centred_regions[keep_regions],
      value_column = "score",
      extend = extend,
      w = bin_size,
      mean_mode = "coverage",
      background = 0,
      smooth = TRUE
    )

    signal_matrix[
      keep_regions,
      seq_len(ncol(matrix_common))
    ] <- as.matrix(matrix_common)
  }

  return(signal_matrix)
}


make_profiles <- function(
    FOXL2_regions,
    FOXL2_table,
    TRIM28,
    TF_list,
    bigwig_files,
    extend,
    bin_size) {

  log_message("Generating WT/cKO signal profiles on TRIM28-centred FOXL2 hubs.")

  TRIM28_FOXL2 <- TRIM28[
    count_overlap_binary(
      TRIM28,
      FOXL2_regions
    )
  ]

  if (length(TRIM28_FOXL2) == 0) {
    log_error(
      "No TRIM28 regions overlap quantified regions associated with FOXL2 WT peaks."
    )
  }

  hub_table <- data.frame(
    NR5A2 = count_overlap_binary(
      TRIM28_FOXL2,
      TF_list$NR5A2
    ),
    ESR2 = count_overlap_binary(
      TRIM28_FOXL2,
      TF_list$ESR2
    ),
    RUNX = count_overlap_binary(
      TRIM28_FOXL2,
      TF_list$RUNX
    ),
    stringsAsFactors = FALSE
  )

  hub_table$n_additional_TFs <- rowSums(
    hub_table[
      ,
      c(
        "NR5A2",
        "ESR2",
        "RUNX"
      )
    ]
  )

  profile_tables <- list()

  for (mark in names(bigwig_files)) {

    log_message("Processing profiles for ", mark)

    for (condition in names(bigwig_files[[mark]])) {

      log_message(
        "Importing ",
        mark,
        " ",
        condition,
        " signal."
      )

      matrix <- extract_profile_matrix(
        bigwig_file = bigwig_files[[mark]][[condition]],
        regions = TRIM28_FOXL2,
        extend = extend,
        bin_size = bin_size
      )

      positions <- seq(
        -extend,
        extend,
        length.out = ncol(matrix)
      )

      for (n_TFs in 0:3) {

        selected <- hub_table$n_additional_TFs == n_TFs

        if (sum(selected) == 0) {
          next
        }

        profile <- colMeans(
          matrix[
            selected,
            ,
            drop = FALSE
          ],
          na.rm = TRUE
        )

        profile_tables[[
          paste(
            mark,
            condition,
            n_TFs,
            sep = "_"
          )
        ]] <- data.frame(
          mark = mark,
          condition = condition,
          n_additional_TFs = n_TFs,
          position = positions,
          signal = profile,
          n_regions = sum(selected),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  return(
    do.call(
      rbind,
      profile_tables
    )
  )
}



###########################################
# Plot functions
###########################################

plot_TRIM28_delta <- function(
    FOXL2_table) {

  ggplot(
    FOXL2_table,
    aes(
      x = TRIM28_class,
      y = delta_FOXL2
    )
  ) +
    geom_violin(
      fill = "#D9D9D9",
      colour = "#555555",
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      fill = "white"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "#777777"
    ) +
    xlab(NULL) +
    ylab("Relative FOXL2 occupancy\nlog2(cKO / WT)") +
    ggtitle(
      "FOXL2 occupancy according to TRIM28 binding"
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      )
    )
}


plot_hub_delta <- function(
    FOXL2_table,
    stats_table) {

  data <- FOXL2_table[
    FOXL2_table$TRIM28,
    ,
    drop = FALSE
  ]

  spearman <- stats_table[
    stats_table$test ==
      "Delta FOXL2 vs number of additional TFs at TRIM28+ sites",
    ,
    drop = FALSE
  ]

  subtitle <- paste0(
    "Region-level Spearman rho = ",
    round(
      spearman$statistic,
      2
    ),
    ", P = ",
    formatC(
      spearman$p_value,
      format = "e",
      digits = 2
    )
  )

  ggplot(
    data,
    aes(
      x = factor(
        n_additional_TFs,
        levels = 0:3
      ),
      y = delta_FOXL2
    )
  ) +
    geom_violin(
      fill = "#D9D9D9",
      colour = "#555555",
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      fill = "white"
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "#777777"
    ) +
    xlab(
      "Number of additional ovarian TFs at TRIM28+ sites"
    ) +
    ylab(
      "Relative FOXL2 occupancy\nlog2(cKO / WT)"
    ) +
    ggtitle(
      "FOXL2 occupancy according to regulatory co-occupancy",
      subtitle = subtitle
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    )
}


plot_loss_enrichment <- function(
    enrichment_table) {

  data <- enrichment_table

  data$feature <- factor(
    data$feature,
    levels = rev(
      c(
        "TRIM28",
        "NR5A2",
        "ESR2",
        "RUNX",
        "TRIM28 + >=1 TF",
        "TRIM28 + >=2 TFs",
        "TRIM28 + 3 TFs"
      )
    )
  )

  ggplot(
    data,
    aes(
      x = odds_ratio,
      y = feature
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "#777777"
    ) +
    geom_errorbarh(
      aes(
        xmin = conf_low,
        xmax = conf_high
      ),
      height = 0.2
    ) +
    geom_point(
      size = 2.5
    ) +
    scale_x_log10() +
    xlab(
      "Odds ratio for enrichment in strongest FOXL2-loss regions"
    ) +
    ylab(NULL) +
    ggtitle(
      paste0(
        "Regulatory features enriched in the strongest ",
        loss_fraction * 100,
        "% of FOXL2-loss regions"
      )
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


plot_profiles <- function(
    profiles,
    mark) {

  data <- profiles[
    profiles$mark == mark,
    ,
    drop = FALSE
  ]

  class_labels <- unique(
    data[
      ,
      c(
        "n_additional_TFs",
        "n_regions"
      )
    ]
  )

  class_labels$facet <- paste0(
    class_labels$n_additional_TFs,
    " additional TF",
    ifelse(
      class_labels$n_additional_TFs == 1,
      "",
      "s"
    ),
    "\n(n = ",
    class_labels$n_regions,
    ")"
  )

  data <- merge(
    data,
    class_labels,
    by = c(
      "n_additional_TFs",
      "n_regions"
    ),
    all.x = TRUE
  )

  data$facet <- factor(
    data$facet,
    levels = class_labels$facet[
      order(
        class_labels$n_additional_TFs
      )
    ]
  )

  ggplot(
    data,
    aes(
      x = position,
      y = signal,
      linetype = condition
    )
  ) +
    geom_line(
      linewidth = 0.8,
      colour = "#320A28"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "#777777"
    ) +
    facet_wrap(
      ~facet,
      nrow = 1,
      scales = "fixed"
    ) +
    scale_linetype_manual(
      values = c(
        "WT" = "solid",
        "cKO" = "dashed"
      )
    ) +
    scale_x_continuous(
      breaks = c(
        -profile_window,
        0,
        profile_window
      ),
      labels = c(
        paste0(
          "-",
          profile_window / 1000,
          " kb"
        ),
        "centre",
        paste0(
          "+",
          profile_window / 1000,
          " kb"
        )
      )
    ) +
    xlab(NULL) +
    ylab("Mean coverage") +
    ggtitle(
      paste0(
        mark,
        " signal at TRIM28+FOXL2 regulatory regions"
      )
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      strip.text = element_text(
        face = "bold"
      ),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}


###########################################
# Run FOXL2 occupancy analysis
###########################################

FOXL2_data <- prepare_FOXL2_counts(
  FOXL2_WT = FOXL2_WT,
  count_file = FOXL2_count_file,
  WT_column = WT_column,
  cKO_column = cKO_column,
  pseudocount = pseudocount
)

FOXL2_regions <- FOXL2_data$regions
FOXL2_table <- FOXL2_data$table

FOXL2_table <- annotate_FOXL2_regions(
  FOXL2_regions = FOXL2_regions,
  FOXL2_table = FOXL2_table,
  TRIM28 = TRIM28,
  TF_list = TF_list
)

FOXL2_table <- define_FOXL2_loss(
  FOXL2_table,
  loss_fraction = loss_fraction
)

summary_table <- summarise_delta(
  FOXL2_table
)

stats_table <- run_continuous_statistics(
  FOXL2_table
)

enrichment_table <- run_loss_enrichment(
  FOXL2_table
)


###########################################
# Generate signal profiles
###########################################

bigwig_files <- list(
  FOXL2 = list(
    WT = FOXL2_WT_bw,
    cKO = FOXL2_cKO_bw
  ),
  ATAC = list(
    WT = ATAC_WT_bw,
    cKO = ATAC_cKO_bw
  ),
  SUMO1 = list(
    WT = SUMO1_WT_bw,
    cKO = SUMO1_cKO_bw
  ),
  SUMO2 = list(
    WT = SUMO2_WT_bw,
    cKO = SUMO2_cKO_bw
  )
)

profiles <- make_profiles(
  FOXL2_regions = FOXL2_regions,
  FOXL2_table = FOXL2_table,
  TRIM28 = TRIM28,
  TF_list = TF_list,
  bigwig_files = bigwig_files,
  extend = profile_window,
  bin_size = bin_size
)


###########################################
# Generate figures
###########################################

log_message("Generating FOXL2 occupancy figures.")

plot_A <- plot_TRIM28_delta(
  FOXL2_table
)

plot_B <- plot_hub_delta(
  FOXL2_table,
  stats_table
)

plot_enrichment <- plot_loss_enrichment(
  enrichment_table
)

plot_C <- plot_profiles(
  profiles,
  "FOXL2"
)

plot_D <- plot_profiles(
  profiles,
  "ATAC"
)

plot_E <- plot_profiles(
  profiles,
  "SUMO1"
)

plot_F <- plot_profiles(
  profiles,
  "SUMO2"
)

figure_top <- cowplot::plot_grid(
  plot_A,
  plot_B,
  plot_enrichment,
  ncol = 3,
  labels = c(
    "A",
    "B",
    "C"
  ),
  label_size = 12
)

figure_profiles <- cowplot::plot_grid(
  plot_C,
  plot_D,
  plot_E,
  plot_F,
  ncol = 1,
  labels = c(
    "D",
    "E",
    "F",
    "G"
  ),
  label_size = 12
)

figure <- cowplot::plot_grid(
  figure_top,
  figure_profiles,
  ncol = 1,
  rel_heights = c(
    1,
    3
  )
)


###########################################
# Export tables
###########################################

log_message("Writing FOXL2 region table: ", output_regions)

write.table(
  FOXL2_table,
  output_regions,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing FOXL2 summary table: ", output_summary)

write.table(
  summary_table,
  output_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing exploratory statistics: ", output_stats)

write.table(
  stats_table,
  output_stats,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing FOXL2-loss enrichment statistics: ", output_enrichment)

write.table(
  enrichment_table,
  output_enrichment,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing signal profiles: ", output_profiles)

write.table(
  profiles,
  output_profiles,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


###########################################
# Export figure
###########################################

log_message("Saving PDF figure: ", output_pdf)

cowplot::save_plot(
  output_pdf,
  figure,
  base_width = 25,
  base_height = 35,
  units = "cm",
  dpi = 300
)

log_message("Saving PNG figure: ", output_png)

cowplot::save_plot(
  output_png,
  figure,
  base_width = 25,
  base_height = 35,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message(
  "Exploratory FOXL2 occupancy analysis completed."
)