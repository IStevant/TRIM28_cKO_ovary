source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting TRIM28-positive DAR-a profile and trend analysis")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("IRanges")
  library("rtracklayer")
  library("EnrichedHeatmap")
  library("ggplot2")
  library("cowplot")
})

###########################################
# Parameters
###########################################

ATAC_bw_folder <- snakemake@params[["ATAC_bw"]]
H3K27ac_bw_folder <- snakemake@params[["H3K27ac_bw"]]
ChIP_bw_folder <- snakemake@params[["ChIP_bw"]]
SUMO <- snakemake@params[["SUMO"]]

ATAC_rep_bw_folder <- dirname(normalizePath(ATAC_bw_folder, mustWork = FALSE))
H3K27ac_rep_bw_folder <- dirname(normalizePath(H3K27ac_bw_folder, mustWork = FALSE))
ChIP_rep_bw_folder <- dirname(normalizePath(ChIP_bw_folder, mustWork = FALSE))

distance_to_feature <- if ("distance_to_feature" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["distance_to_feature"]])
} else {
  1000
}

central_window <- if ("central_window" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["central_window"]])
} else {
  250
}

bin_size <- if ("bin_size" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["bin_size"]])
} else {
  50
}

###########################################
# Load inputs
###########################################

samplesheet <- read.csv(
  snakemake@input[["samplesheet"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

conditions <- unique(gsub("_R.*", "", samplesheet$sample))

control_condition <- conditions[
  !grepl("ko|cko|mut", conditions, ignore.case = TRUE)
][1]

cko_condition <- conditions[
  grepl("ko|cko|mut", conditions, ignore.case = TRUE)
][1]

if (is.na(control_condition) || is.na(cko_condition)) {
  control_condition <- conditions[1]
  cko_condition <- conditions[2]
}

log_message("Control condition: ", control_condition)
log_message("cKO condition: ", cko_condition)

ATAC_clustering <- read.csv(
  snakemake@input[["ATAC_clustering"]],
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("region", "cluster") %in% colnames(ATAC_clustering))) {
  log_error("ATAC_clustering must contain columns 'region' and 'cluster'")
}

ATAC_clustering <- ATAC_clustering[
  grepl("^a", tolower(as.character(ATAC_clustering$cluster))),
  ,
  drop = FALSE
]

if (nrow(ATAC_clustering) == 0) {
  log_error("No DAR cluster a regions found")
}

TRIM28_peaks <- rtracklayer::import(snakemake@input[["TRIM28"]])

TF_list <- list(
  FOXL2 = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2 = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2  = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX  = rtracklayer::import(snakemake@input[["RUNX"]])
)

###########################################
# BigWig helpers
###########################################

find_bigwig_by_condition <- function(
  folder,
  condition,
  mark_name = NULL
) {
  matched_files <- list.files(
    folder,
    pattern = condition,
    full.names = TRUE,
    recursive = TRUE
  )

  matched_files <- matched_files[
    grepl("\\.bw$|\\.bigWig$", matched_files, ignore.case = TRUE)
  ]

  if (!is.null(mark_name)) {
    mark_matches <- matched_files[
      grepl(mark_name, basename(matched_files), ignore.case = TRUE)
    ]

    if (length(mark_matches) > 0) {
      matched_files <- mark_matches
    }
  }

  if (length(matched_files) != 1) {
    log_error(
      "Expected exactly one merged bigWig for ",
      condition,
      " / ",
      ifelse(is.null(mark_name), "signal", mark_name),
      ", found: ",
      paste(basename(matched_files), collapse = ", ")
    )
  }

  matched_files
}

get_merged_bigwigs <- function(folder, mark_name = NULL) {
  files <- sapply(
    conditions,
    function(condition) {
      find_bigwig_by_condition(
        folder = folder,
        condition = condition,
        mark_name = mark_name
      )
    }
  )

  names(files) <- conditions
  files
}

get_replicate_bigwigs <- function(
  folder,
  mark_name = NULL,
  replicate_pattern
) {
  all_files <- list()

  for (condition in conditions) {
    matched_files <- list.files(
      folder,
      pattern = condition,
      full.names = TRUE,
      recursive = TRUE
    )

    matched_files <- matched_files[
      grepl("\\.bw$|\\.bigWig$", matched_files, ignore.case = TRUE)
    ]

    if (!is.null(mark_name)) {
      mark_matches <- matched_files[
        grepl(mark_name, basename(matched_files), ignore.case = TRUE)
      ]

      if (length(mark_matches) > 0) {
        matched_files <- mark_matches
      }
    }

    matched_files <- matched_files[
      grepl(replicate_pattern, basename(matched_files), ignore.case = FALSE)
    ]

    if (length(matched_files) == 0) {
      log_error(
        "No replicate bigWigs found for ",
        condition,
        " / ",
        ifelse(is.null(mark_name), "signal", mark_name),
        " using pattern ",
        replicate_pattern
      )
    }

    replicate_names <- sub(
      "\\.(bw|bigWig)$",
      "",
      basename(matched_files),
      ignore.case = TRUE
    )

    names(matched_files) <- replicate_names

    all_files <- c(all_files, as.list(matched_files))
  }

  unlist(all_files)
}

merged_bw_files <- c(
  list(
    ATAC = get_merged_bigwigs(ATAC_bw_folder, "ATAC"),
    H3K27ac = get_merged_bigwigs(H3K27ac_bw_folder, "H3K27ac")
  ),
  setNames(
    lapply(
      SUMO,
      function(mark_name) {
        get_merged_bigwigs(
          folder = ChIP_bw_folder,
          mark_name = mark_name
        )
      }
    ),
    SUMO
  )
)

replicate_bw_files <- c(
  list(
    ATAC = get_replicate_bigwigs(
      folder = ATAC_rep_bw_folder,
      mark_name = "ATAC",
      replicate_pattern = "REP[0-9]+"
    ),
    H3K27ac = get_replicate_bigwigs(
      folder = H3K27ac_rep_bw_folder,
      mark_name = "H3K27ac",
      replicate_pattern = "R[0-9]+"
    )
  ),
  setNames(
    lapply(
      SUMO,
      function(mark_name) {
        get_replicate_bigwigs(
          folder = ChIP_rep_bw_folder,
          mark_name = mark_name,
          replicate_pattern = "Rep[0-9]+"
        )
      }
    ),
    SUMO
  )
)

###########################################
# Functions
###########################################

round_down <- function(x) {
  trunc(x * 100) / 100
}

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

annotate_TRIM28_regions <- function(
  TRIM28_peaks,
  TF_list
) {
  log_message("Annotating TRIM28 peaks with TF overlaps")

  TF_overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      count_overlap_binary(
        query = TRIM28_peaks,
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

  TRIM28_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(TRIM28_peaks)),
    start = GenomicRanges::start(TRIM28_peaks),
    end = GenomicRanges::end(TRIM28_peaks),
    region = make_region_id(TRIM28_peaks),
    TF_overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = TF_combination,
    stringsAsFactors = FALSE
  )

  TRIM28_table$TRIM28_TF_class <- ifelse(
    TRIM28_table$n_TFs == 0,
    "TRIM28 only",
    ifelse(
      TRIM28_table$n_TFs == 1,
      "TRIM28 + 1 TF",
      ifelse(
        TRIM28_table$n_TFs == 2,
        "TRIM28 + 2 TFs",
        ifelse(
          TRIM28_table$n_TFs == 3,
          "TRIM28 + 3 TFs",
          "TRIM28 + 4 TFs"
        )
      )
    )
  )

  TRIM28_table$TRIM28_TF_class <- factor(
    TRIM28_table$TRIM28_TF_class,
    levels = c(
      "TRIM28 only",
      "TRIM28 + 1 TF",
      "TRIM28 + 2 TFs",
      "TRIM28 + 3 TFs",
      "TRIM28 + 4 TFs"
    )
  )

  TRIM28_table
}

prepare_DAR_TRIM28_positive_regions <- function(
  ATAC_clustering,
  TRIM28_peaks,
  TRIM28_table
) {
  log_message("Preparing TRIM28-positive DAR-a regions")

  DAR_regions <- GenomicRanges::GRanges(
    ATAC_clustering$region
  )

  DAR_regions$region <- ATAC_clustering$region
  DAR_regions$DAR_cluster <- ATAC_clustering$cluster

  hits <- GenomicRanges::findOverlaps(
    DAR_regions,
    TRIM28_peaks,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    log_error("No DAR-a regions overlap TRIM28 peaks")
  }

  hit_table <- data.frame(
    DAR_index = S4Vectors::queryHits(hits),
    TRIM28_index = S4Vectors::subjectHits(hits),
    n_TFs = TRIM28_table$n_TFs[S4Vectors::subjectHits(hits)],
    TF_combination = TRIM28_table$TF_combination[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )

  best_hit <- do.call(
    rbind,
    lapply(
      split(hit_table, hit_table$DAR_index),
      function(x) {
        x[
          which.max(x$n_TFs),
          ,
          drop = FALSE
        ]
      }
    )
  )

  selected_regions <- DAR_regions[
    as.numeric(best_hit$DAR_index)
  ]

  selected_table <- data.frame(
    region = selected_regions$region,
    DAR_cluster = selected_regions$DAR_cluster,
    n_TFs = best_hit$n_TFs,
    TF_combination = best_hit$TF_combination,
    stringsAsFactors = FALSE
  )

  selected_table$TRIM28_TF_class <- ifelse(
    selected_table$n_TFs == 0,
    "TRIM28 only",
    ifelse(
      selected_table$n_TFs == 1,
      "TRIM28 + 1 TF",
      ifelse(
        selected_table$n_TFs == 2,
        "TRIM28 + 2 TFs",
        ifelse(
          selected_table$n_TFs == 3,
          "TRIM28 + 3 TFs",
          "TRIM28 + 4 TFs"
        )
      )
    )
  )

  selected_table$TRIM28_TF_class <- factor(
    selected_table$TRIM28_TF_class,
    levels = c(
      "TRIM28 only",
      "TRIM28 + 1 TF",
      "TRIM28 + 2 TFs",
      "TRIM28 + 3 TFs",
      "TRIM28 + 4 TFs"
    )
  )

  selected_regions$TRIM28_TF_class <- selected_table$TRIM28_TF_class

  log_message("TRIM28-positive DAR-a regions: ", length(selected_regions))

  for (class_name in levels(selected_table$TRIM28_TF_class)) {
    log_message(
      class_name,
      ": ",
      sum(selected_table$TRIM28_TF_class == class_name)
    )
  }

  list(
    regions = selected_regions,
    table = selected_table
  )
}

normalise_signal_matrix <- function(
  bigwig_file,
  regions,
  distance_to_feature,
  bin_size
) {
  log_message("Importing ", basename(bigwig_file))

  reads <- rtracklayer::import(
    bigwig_file
  )

  EnrichedHeatmap::normalizeToMatrix(
    reads,
    regions,
    background = 0,
    extend = distance_to_feature,
    w = bin_size,
    mean_mode = "coverage",
    value_column = "score",
    smooth = TRUE
  )
}

make_profile_table_for_mark <- function(
  mark_name,
  bigwig_files,
  regions,
  classes
) {
  all_profiles <- list()

  for (condition in conditions) {
    mat <- normalise_signal_matrix(
      bigwig_file = bigwig_files[[condition]],
      regions = regions,
      distance_to_feature = distance_to_feature,
      bin_size = bin_size
    )

    mat <- as.matrix(mat)

    position <- seq(
      from = -distance_to_feature,
      to = distance_to_feature,
      length.out = ncol(mat)
    )

    for (class_name in levels(classes)) {
      class_index <- classes == class_name

      if (sum(class_index, na.rm = TRUE) == 0) {
        next
      }

      profile <- colMeans(
        mat[class_index, , drop = FALSE],
        na.rm = TRUE
      )

      all_profiles[[paste(mark_name, condition, class_name, sep = "_")]] <- data.frame(
        mark = mark_name,
        condition = condition,
        TRIM28_TF_class = class_name,
        n_TFs = unique(selected_table$n_TFs[selected_table$TRIM28_TF_class == class_name])[1],
        position = position,
        signal = profile,
        n_regions = sum(class_index, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, all_profiles)
}

make_replicate_signal_table_for_mark <- function(
  mark_name,
  bigwig_files,
  regions,
  classes
) {
  all_signals <- list()

  for (replicate_name in names(bigwig_files)) {
    condition <- conditions[
      vapply(
        conditions,
        function(x) grepl(x, replicate_name, ignore.case = TRUE),
        logical(1)
      )
    ][1]

    if (is.na(condition)) {
      log_error("Could not infer condition for replicate: ", replicate_name)
    }

    mat <- normalise_signal_matrix(
      bigwig_file = bigwig_files[[replicate_name]],
      regions = regions,
      distance_to_feature = distance_to_feature,
      bin_size = bin_size
    )

    mat <- as.matrix(mat)

    position <- seq(
      from = -distance_to_feature,
      to = distance_to_feature,
      length.out = ncol(mat)
    )

    central_columns <- abs(position) <= central_window

    central_signal_per_region <- rowMeans(
      mat[, central_columns, drop = FALSE],
      na.rm = TRUE
    )

    for (class_name in levels(classes)) {
      class_index <- classes == class_name

      if (sum(class_index, na.rm = TRUE) == 0) {
        next
      }

      all_signals[[paste(mark_name, replicate_name, class_name, sep = "_")]] <- data.frame(
        mark = mark_name,
        condition = condition,
        replicate = replicate_name,
        TRIM28_TF_class = class_name,
        n_TFs = unique(selected_table$n_TFs[selected_table$TRIM28_TF_class == class_name])[1],
        signal = mean(
          central_signal_per_region[class_index],
          na.rm = TRUE
        ),
        n_regions = sum(class_index, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, all_signals)
}

summarise_replicate_delta <- function(
  replicate_signal_table
) {
  stats_table <- do.call(
    rbind,
    lapply(
      split(
        replicate_signal_table,
        paste(
          replicate_signal_table$mark,
          replicate_signal_table$TRIM28_TF_class,
          sep = "__"
        )
      ),
      function(x) {
        ctrl_values <- x$signal[x$condition == control_condition]
        cko_values <- x$signal[x$condition == cko_condition]

        mean_ctrl <- mean(ctrl_values, na.rm = TRUE)
        mean_cko <- mean(cko_values, na.rm = TRUE)

        data.frame(
          mark = unique(x$mark),
          TRIM28_TF_class = unique(x$TRIM28_TF_class),
          n_TFs = unique(x$n_TFs),
          n_regions = unique(x$n_regions),
          n_ctrl_replicates = length(ctrl_values),
          n_cko_replicates = length(cko_values),
          mean_ctrl = mean_ctrl,
          mean_cko = mean_cko,
          delta = mean_cko - mean_ctrl,
          delta_pct = 100 * (mean_cko - mean_ctrl) / mean_ctrl,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(stats_table) <- NULL
  stats_table
}

run_trend_correlations <- function(
  stats_table
) {
  marks_to_test <- c(
    "H3K27ac",
    SUMO
  )

  correlation_table <- do.call(
    rbind,
    lapply(
      marks_to_test,
      function(mark_name) {
        x <- stats_table[
          stats_table$mark == mark_name,
          ,
          drop = FALSE
        ]

        x <- x[
          order(x$n_TFs),
          ,
          drop = FALSE
        ]

        test <- stats::cor.test(
          x$n_TFs,
          x$delta_pct,
          method = "spearman",
          alternative = "less",
          exact = FALSE
        )

        data.frame(
          mark = mark_name,
          rho = unname(test$estimate),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  correlation_table
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

make_condition_linetypes <- function(conditions) {
  linetypes <- ifelse(
    conditions == cko_condition,
    "dashed",
    "solid"
  )

  names(linetypes) <- conditions
  linetypes
}

make_condition_colours <- function(conditions) {
  colours <- ifelse(
    conditions == cko_condition,
    "#7D1128",
    "#333333"
  )

  names(colours) <- conditions
  colours
}

plot_mark_profiles <- function(
  profile_table,
  stats_table,
  correlation_table,
  mark_name
) {
  plot_table <- profile_table[
    profile_table$mark == mark_name,
    ,
    drop = FALSE
  ]

  stats_subset <- stats_table[
    stats_table$mark == mark_name,
    ,
    drop = FALSE
  ]

  class_labels <- unique(
    plot_table[
      ,
      c("TRIM28_TF_class", "n_regions")
    ]
  )

  class_labels$facet_label <- paste0(
    class_labels$TRIM28_TF_class,
    "\n(n = ",
    class_labels$n_regions,
    ")"
  )

  plot_table <- merge(
    plot_table,
    class_labels,
    by = c(
      "TRIM28_TF_class",
      "n_regions"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  stats_subset <- merge(
    stats_subset,
    class_labels,
    by = c(
      "TRIM28_TF_class",
      "n_regions"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  class_order <- c(
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )

  facet_levels <- class_labels$facet_label[
    match(
      class_order,
      class_labels$TRIM28_TF_class
    )
  ]

  facet_levels <- facet_levels[!is.na(facet_levels)]

  plot_table$facet_label <- factor(
    plot_table$facet_label,
    levels = facet_levels
  )

  stats_subset$facet_label <- factor(
    stats_subset$facet_label,
    levels = facet_levels
  )

  y_max <- aggregate(
    signal ~ facet_label,
    data = plot_table,
    FUN = max
  )

  colnames(y_max) <- c(
    "facet_label",
    "y"
  )

  stats_subset <- merge(
    stats_subset,
    y_max,
    by = "facet_label",
    all.x = TRUE
  )

  stats_subset$x <- 0
  stats_subset$y <- stats_subset$y * 0.92

  stats_subset$label <- paste0(
    "delta = ",
    round(stats_subset$delta_pct, 1),
    "%"
  )

  subtitle <- ""

  if (mark_name %in% correlation_table$mark) {
    corr <- correlation_table[
      correlation_table$mark == mark_name,
      ,
      drop = FALSE
    ]

    subtitle <- paste0(
      "Spearman trend across TRIM28-TF classes: rho = ",
      round_down(corr$rho),
      ", ",
      format_pvalue(corr$p_value)
    )
  }

  ggplot(
    plot_table,
    aes(
      x = position,
      y = signal,
      colour = condition,
      linetype = condition
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "#555555",
      linewidth = 0.3
    ) +
    geom_text(
      data = stats_subset,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      hjust = 0.5
    ) +
    facet_wrap(
      ~facet_label,
      ncol = 5
    ) +
    scale_colour_manual(
      values = make_condition_colours(conditions)
    ) +
    scale_linetype_manual(
      values = make_condition_linetypes(conditions)
    ) +
    scale_x_continuous(
      breaks = c(
        -distance_to_feature,
        0,
        distance_to_feature
      ),
      labels = c(
        paste0("-", distance_to_feature / 1000, " kb"),
        "centre",
        paste0("+", distance_to_feature / 1000, " kb")
      )
    ) +
    xlab(NULL) +
    ylab("Mean coverage") +
    ggtitle(
      paste0(
        mark_name,
        " coverage at TRIM28-positive DAR-a regions"
      ),
      subtitle = subtitle
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 12
      ),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
    )
}

###########################################
# Prepare regions
###########################################

TRIM28_table <- annotate_TRIM28_regions(
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

selected <- prepare_DAR_TRIM28_positive_regions(
  ATAC_clustering = ATAC_clustering,
  TRIM28_peaks = TRIM28_peaks,
  TRIM28_table = TRIM28_table
)

selected_regions <- selected$regions
selected_table <- selected$table

classes <- factor(
  selected_table$TRIM28_TF_class,
  levels = c(
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )
)

keep_levels <- levels(classes)[
  levels(classes) %in% unique(as.character(classes))
]

classes <- factor(
  as.character(classes),
  levels = keep_levels
)

selected_regions$TRIM28_TF_class <- classes

###########################################
# Compute merged profiles
###########################################

profile_tables <- lapply(
  names(merged_bw_files),
  function(mark_name) {
    make_profile_table_for_mark(
      mark_name = mark_name,
      bigwig_files = merged_bw_files[[mark_name]],
      regions = selected_regions,
      classes = classes
    )
  }
)

profile_table <- do.call(rbind, profile_tables)

###########################################
# Compute replicate-based deltas
###########################################

replicate_signal_tables <- lapply(
  names(replicate_bw_files),
  function(mark_name) {
    make_replicate_signal_table_for_mark(
      mark_name = mark_name,
      bigwig_files = replicate_bw_files[[mark_name]],
      regions = selected_regions,
      classes = classes
    )
  }
)

replicate_signal_table <- do.call(rbind, replicate_signal_tables)

stats_table <- summarise_replicate_delta(
  replicate_signal_table = replicate_signal_table
)

correlation_table <- run_trend_correlations(
  stats_table = stats_table
)

###########################################
# Save tables
###########################################

# write.table(
#   selected_table,
#   file = snakemake@output[["regions"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

# write.table(
#   profile_table,
#   file = snakemake@output[["profiles"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

# write.table(
#   replicate_signal_table,
#   file = snakemake@output[["replicate_signals"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

# write.table(
#   stats_table,
#   file = snakemake@output[["stats"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

# write.table(
#   correlation_table,
#   file = snakemake@output[["correlations"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

###########################################
# Save figures
###########################################

mark_plots <- lapply(
  names(merged_bw_files),
  function(mark_name) {
    plot_mark_profiles(
      profile_table = profile_table,
      stats_table = stats_table,
      correlation_table = correlation_table,
      mark_name = mark_name
    )
  }
)

names(mark_plots) <- names(merged_bw_files)

pdf(
  file = snakemake@output[["pdf"]],
  width = 22 / 2.54,
  height = 8 / 2.54
)

for (mark_name in names(mark_plots)) {
  print(mark_plots[[mark_name]])
}

dev.off()

summary_plot <- cowplot::plot_grid(
  plotlist = mark_plots,
  ncol = 1,
  labels = names(mark_plots),
  label_size = 12
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = summary_plot,
  base_width = 22,
  base_height = 8 * length(mark_plots),
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")