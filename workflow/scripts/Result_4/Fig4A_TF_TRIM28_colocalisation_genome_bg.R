source(".Rprofile")

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}

log_message("Starting pairwise TF colocalisation analysis with precomputed genome background")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomeInfoDb")
  library("S4Vectors")
  library("rtracklayer")
  library("regioneR")
  library("ggplot2")
  library("cowplot")
})

###########################################
# Parameters
###########################################

n_permutations <- if ("n_permutations" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["n_permutations"]])
} else {
  1000
}

random_seed <- if ("random_seed" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["random_seed"]])
} else {
  1234
}

if (is.na(n_permutations) || n_permutations < 1) {
  log_error("Parameter 'n_permutations' must be a positive numeric value")
}

if (is.na(random_seed)) {
  log_error("Parameter 'random_seed' must be numeric")
}

###########################################
# Load input data
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

log_message("Importing TF / cofactor peak files")

peak_list <- list(
  TRIM28 = rtracklayer::import(snakemake@input[["TRIM28"]]),
  FOXL2  = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2  = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2   = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX  = rtracklayer::import(snakemake@input[["RUNX"]])
)

for (peak_name in names(peak_list)) {
  log_message("Number of ", peak_name, " peaks: ", length(peak_list[[peak_name]]))
}

###########################################
# Functions
###########################################

prepare_genome_background <- function(
  genome_gtf,
  peak_list
) {
  log_message("Preparing genome background")

  seq_lengths <- GenomeInfoDb::seqlengths(
    GenomeInfoDb::seqinfo(genome_gtf)
  )

  seq_lengths <- seq_lengths[
    !is.na(seq_lengths) &
      seq_lengths > 0
  ]

  if (length(seq_lengths) == 0) {
    log_message("No chromosome lengths found in GTF seqinfo; using maximum observed coordinates")

    all_peak_regions <- unlist(
      GenomicRanges::GRangesList(peak_list),
      use.names = FALSE
    )

    all_regions <- c(
      GenomicRanges::GRanges(genome_gtf),
      all_peak_regions
    )

    seq_lengths <- tapply(
      GenomicRanges::end(all_regions),
      as.character(GenomicRanges::seqnames(all_regions)),
      max,
      na.rm = TRUE
    )
  }

  seq_lengths <- seq_lengths[
    grepl("^chr([0-9]+|X|Y)$", names(seq_lengths))
  ]

  if (length(seq_lengths) == 0) {
    log_error("No valid standard chromosomes available for genome background")
  }

  genome <- data.frame(
    chr = names(seq_lengths),
    start = 1,
    end = as.numeric(seq_lengths),
    stringsAsFactors = FALSE
  )

  log_message(
    "Chromosomes used as genome background: ",
    paste(genome$chr, collapse = ", ")
  )

  genome
}

keep_standard_chromosomes <- function(
  regions,
  genome
) {
  regions <- regions[
    as.character(GenomicRanges::seqnames(regions)) %in% genome$chr
  ]

  GenomeInfoDb::seqlevels(regions) <- intersect(
    GenomeInfoDb::seqlevels(regions),
    genome$chr
  )

  regions
}

count_overlapping_regions <- function(
  query_regions,
  subject_regions
) {
  sum(
    GenomicRanges::countOverlaps(
      query_regions,
      subject_regions,
      ignore.strand = TRUE
    ) > 0
  )
}

precompute_random_backgrounds <- function(
  peak_list,
  genome,
  n_permutations,
  random_seed
) {
  log_message("Precomputing random genomic backgrounds")

  random_backgrounds <- lapply(names(peak_list), function(peak_name) {
    log_message("Randomising peaks for ", peak_name)

    set.seed(random_seed + match(peak_name, names(peak_list)))

    random_sets <- lapply(seq_len(n_permutations), function(i) {
      regioneR::randomizeRegions(
        A = peak_list[[peak_name]],
        genome = genome,
        per.chromosome = TRUE,
        allow.overlaps = TRUE
      )
    })

    names(random_sets) <- paste0("perm_", seq_len(n_permutations))

    random_sets
  })

  names(random_backgrounds) <- names(peak_list)

  random_backgrounds
}

run_pairwise_test_from_precomputed_background <- function(
  query_regions,
  target_regions,
  query_name,
  target_name,
  random_query_regions,
  n_permutations
) {
  log_message("Testing pair: ", query_name, " -> ", target_name)

  observed_overlap <- count_overlapping_regions(
    query_regions = query_regions,
    subject_regions = target_regions
  )

  random_overlaps <- vapply(
    random_query_regions,
    function(random_regions) {
      count_overlapping_regions(
        query_regions = random_regions,
        subject_regions = target_regions
      )
    },
    numeric(1)
  )

  expected_overlap <- mean(random_overlaps)
  expected_sd <- stats::sd(random_overlaps)

  empirical_pvalue <- (
    sum(random_overlaps >= observed_overlap) + 1
  ) / (
    n_permutations + 1
  )

  fold_enrichment <- ifelse(
    expected_overlap > 0,
    observed_overlap / expected_overlap,
    NA_real_
  )

  z_score <- ifelse(
    expected_sd > 0,
    (observed_overlap - expected_overlap) / expected_sd,
    NA_real_
  )

  data.frame(
    query_factor = query_name,
    target_factor = target_name,
    pair = paste0(query_name, "_vs_", target_name),
    observed_overlap = observed_overlap,
    query_total = length(query_regions),
    observed_percent = 100 * observed_overlap / length(query_regions),
    expected_overlap = expected_overlap,
    expected_sd = expected_sd,
    expected_percent = 100 * expected_overlap / length(query_regions),
    fold_enrichment = fold_enrichment,
    z_score = z_score,
    empirical_pvalue = empirical_pvalue,
    n_permutations = n_permutations,
    background = "precomputed_genome_randomisation_regioneR",
    stringsAsFactors = FALSE
  )
}

plot_pairwise_enrichment <- function(
  result_table
) {
  result_table$pair_label <- paste0(
    result_table$query_factor,
    " → ",
    result_table$target_factor
  )

  result_table$significance <- ifelse(
    result_table$empirical_pvalue < 0.001,
    "p < 0.001",
    ifelse(
      result_table$empirical_pvalue < 0.01,
      "p < 0.01",
      ifelse(
        result_table$empirical_pvalue < 0.05,
        "p < 0.05",
        "n.s."
      )
    )
  )

  result_table$pair_label <- factor(
    result_table$pair_label,
    levels = result_table$pair_label[
      order(result_table$fold_enrichment)
    ]
  )

  ggplot(
    result_table,
    aes(
      x = pair_label,
      y = fold_enrichment
    )
  ) +
    geom_col(
      fill = "#7D1128",
      colour = "#7D1128",
      linewidth = 0.2
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    geom_text(
      aes(
        label = paste0(
          round(fold_enrichment, 2),
          "x\n",
          significance
        )
      ),
      hjust = -0.05,
      size = 4
    ) +
    coord_flip() +
    ylab("Fold enrichment over genome background") +
    xlab(NULL) +
    ggtitle("Pairwise colocalisation enrichment") +
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
# Prepare genome and peaks
###########################################

genome_background <- prepare_genome_background(
  genome_gtf = genome_gtf,
  peak_list = peak_list
)

peak_list <- lapply(
  peak_list,
  keep_standard_chromosomes,
  genome = genome_background
)

for (peak_name in names(peak_list)) {
  log_message(
    "Number of ",
    peak_name,
    " peaks on standard chromosomes: ",
    length(peak_list[[peak_name]])
  )
}

###########################################
# Precompute randomised peak sets
###########################################

random_backgrounds <- precompute_random_backgrounds(
  peak_list = peak_list,
  genome = genome_background,
  n_permutations = n_permutations,
  random_seed = random_seed
)

###########################################
# Run pairwise tests
###########################################

pair_grid <- expand.grid(
  query_factor = names(peak_list),
  target_factor = names(peak_list),
  stringsAsFactors = FALSE
)

pair_grid <- pair_grid[
  pair_grid$query_factor != pair_grid$target_factor,
  ,
  drop = FALSE
]

pairwise_results <- lapply(seq_len(nrow(pair_grid)), function(i) {
  query_name <- pair_grid$query_factor[i]
  target_name <- pair_grid$target_factor[i]

  run_pairwise_test_from_precomputed_background(
    query_regions = peak_list[[query_name]],
    target_regions = peak_list[[target_name]],
    query_name = query_name,
    target_name = target_name,
    random_query_regions = random_backgrounds[[query_name]],
    n_permutations = n_permutations
  )
})

pairwise_results <- do.call(
  rbind,
  pairwise_results
)

pairwise_results <- pairwise_results[
  order(
    pairwise_results$empirical_pvalue,
    -pairwise_results$fold_enrichment
  ),
  ,
  drop = FALSE
]

###########################################
# Save output table
###########################################

log_message("Writing pairwise enrichment table")

write.table(
  pairwise_results,
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Save plot
###########################################

log_message("Plotting pairwise enrichment results")

figure <- plot_pairwise_enrichment(
  result_table = pairwise_results
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