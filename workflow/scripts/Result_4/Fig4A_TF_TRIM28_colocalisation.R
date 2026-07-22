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

log_message("Starting pairwise TF colocalisation enrichment analysis")

###########################################
# Libraries
###########################################

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomeInfoDb")
  library("S4Vectors")
  library("rtracklayer")
  library("ggplot2")
  library("cowplot")
  library("scales")
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

log_message("Importing ATAC control peaks as background")

ATAC_control_table <- read.csv(
  snakemake@input[["ATAC_control"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

if (!"region" %in% colnames(ATAC_control_table)) {
  log_error("ATAC control table must contain a 'region' column")
}

ATAC_control_peaks <- GenomicRanges::GRanges(
  ATAC_control_table$region
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

log_message("Number of ATAC control peaks: ", length(ATAC_control_peaks))

###########################################
# Functions
###########################################

get_standard_chromosomes <- function(...) {
  region_list <- list(...)

  chromosomes <- unique(unlist(lapply(
    region_list,
    function(regions) {
      as.character(GenomicRanges::seqnames(regions))
    }
  )))

  chromosomes[
    grepl("^chr([0-9]+|X|Y)$", chromosomes)
  ]
}

keep_standard_chromosomes <- function(
  regions,
  standard_chromosomes
) {
  regions <- regions[
    as.character(GenomicRanges::seqnames(regions)) %in% standard_chromosomes
  ]

  GenomeInfoDb::seqlevels(regions) <- intersect(
    GenomeInfoDb::seqlevels(regions),
    standard_chromosomes
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

run_pairwise_permutation_test <- function(
  query_regions,
  target_regions,
  ATAC_control_regions,
  query_name,
  target_name,
  n_permutations,
  random_seed
) {
  set.seed(random_seed)

  log_message("Testing pair: ", query_name, " -> ", target_name)

  observed_overlap <- count_overlapping_regions(
    query_regions = query_regions,
    subject_regions = target_regions
  )

  replace_sampling <- length(ATAC_control_regions) < length(query_regions)

  random_overlaps <- vapply(
    seq_len(n_permutations),
    function(i) {
      random_regions <- ATAC_control_regions[
        sample(
          seq_along(ATAC_control_regions),
          size = length(query_regions),
          replace = replace_sampling
        )
      ]

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
    background = "ATAC_control_resampling",
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
    ylab("Fold enrichment over ATAC background") +
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
# Prepare regions
###########################################

standard_chromosomes <- get_standard_chromosomes(
  ATAC_control_peaks,
  unlist(
    GenomicRanges::GRangesList(peak_list),
    use.names = FALSE
  )
)

ATAC_control_peaks <- keep_standard_chromosomes(
  regions = ATAC_control_peaks,
  standard_chromosomes = standard_chromosomes
)

peak_list <- lapply(
  peak_list,
  keep_standard_chromosomes,
  standard_chromosomes = standard_chromosomes
)

log_message(
  "Standard chromosomes retained: ",
  paste(standard_chromosomes, collapse = ", ")
)

###########################################
# Run pairwise enrichment tests
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

  run_pairwise_permutation_test(
    query_regions = peak_list[[query_name]],
    target_regions = peak_list[[target_name]],
    ATAC_control_regions = ATAC_control_peaks,
    query_name = query_name,
    target_name = target_name,
    n_permutations = n_permutations,
    random_seed = random_seed + i
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