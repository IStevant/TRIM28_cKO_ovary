source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting TF / TRIM28 colocalisation analysis with genome background")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomeInfoDb")
  library("S4Vectors")
  library("rtracklayer")
  library("regioneR")
  library("ChIPseeker")
  library("txdbmaker")
  library("cowplot")
  library("ggplot2")
  library("scales")
})

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

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

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28"]]
)

log_message("Importing TF peak files")

TF_list <- list(
  FOXL2 = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2 = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2  = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX = rtracklayer::import(snakemake@input[["RUNX"]])
)

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

dir.create(
  snakemake@output[["overlaps"]],
  showWarnings = FALSE,
  recursive = TRUE
)

get_tf_colours <- function(TF_name) {
  colours <- list(
    FOXL2 = c("#efc6ee", "#b93ec1"),
    NR5A2 = c("#bceff0", "#329ca3"),
    RUNX = c("#d5e99a", "#86b524"),
    ESR2  = c("#f7d7c9", "#d8581c")
  )

  colours[[TF_name]]
}

###########################################
# Genome annotation
###########################################

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

txdb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

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

prepare_genome_background <- function(genome_gtf, TF_list, TRIM28_peaks) {
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

    all_TF_regions <- unlist(
      GenomicRanges::GRangesList(TF_list),
      use.names = FALSE
    )

    all_regions <- c(
      GenomicRanges::GRanges(genome_gtf),
      TRIM28_peaks,
      all_TF_regions
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

  data.frame(
    chr = names(seq_lengths),
    start = 1,
    end = as.numeric(seq_lengths),
    stringsAsFactors = FALSE
  )
}

keep_standard_chromosomes <- function(regions, genome) {
  regions <- regions[
    as.character(GenomicRanges::seqnames(regions)) %in% genome$chr
  ]

  GenomeInfoDb::seqlevels(regions) <- intersect(
    GenomeInfoDb::seqlevels(regions),
    genome$chr
  )

  regions
}

count_overlapping_regions <- function(query_regions, subject_regions) {
  sum(
    GenomicRanges::countOverlaps(
      query_regions,
      subject_regions,
      ignore.strand = TRUE
    ) > 0
  )
}

run_genome_permutation_test <- function(
  TF_regions,
  TRIM28_regions,
  genome,
  TF_name,
  n_permutations,
  random_seed
) {
  set.seed(random_seed)

  log_message("Running genome permutation test for ", TF_name)

  observed_overlap <- count_overlapping_regions(
    query_regions = TF_regions,
    subject_regions = TRIM28_regions
  )

  test_result <- regioneR::overlapPermTest(
    A = TF_regions,
    B = TRIM28_regions,
    ntimes = n_permutations,
    genome = genome,
    alternative = "greater",
    count.once = TRUE,
    mc.cores = 1
  )

  expected_overlap <- mean(test_result$numOverlaps$permuted)
  expected_sd <- stats::sd(test_result$numOverlaps$permuted)

  fold_enrichment <- ifelse(
    expected_overlap > 0,
    observed_overlap / expected_overlap,
    NA_real_
  )

  data.frame(
    TF = TF_name,
    observed_overlap = observed_overlap,
    TF_total = length(TF_regions),
    observed_percent = 100 * observed_overlap / length(TF_regions),
    expected_overlap = expected_overlap,
    expected_sd = expected_sd,
    expected_percent = 100 * expected_overlap / length(TF_regions),
    fold_enrichment = fold_enrichment,
    z_score = test_result$numOverlaps$zscore,
    empirical_pvalue = test_result$numOverlaps$pval,
    n_permutations = n_permutations,
    background = "genome_randomisation_regioneR",
    stringsAsFactors = FALSE
  )
}

format_pvalue <- function(pvalue, n_permutations) {
  if (is.na(pvalue)) {
    return("empirical p = NA")
  }

  min_pvalue <- 1 / (n_permutations + 1)

  if (pvalue <= min_pvalue) {
    return(paste0("empirical p <= ", signif(min_pvalue, 3)))
  }

  paste0("empirical p = ", signif(pvalue, 3))
}

make_pie_plot <- function(TF_name, test_table) {
  n_overlap <- test_table$observed_overlap
  n_total <- test_table$TF_total
  n_non_overlap <- n_total - n_overlap

  plot_data <- data.frame(
    group = c(
      paste0(TF_name, "-only"),
      paste0(TF_name, " + TRIM28")
    ),
    percent = c(
      100 * n_non_overlap / n_total,
      100 * n_overlap / n_total
    ),
    count = c(
      n_non_overlap,
      n_overlap
    ),
    stringsAsFactors = FALSE
  )

  plot_data$label <- paste0(
    scales::percent(plot_data$percent / 100),
    "\n(",
    scales::comma(plot_data$count),
    ")"
  )

  plot_data$group <- factor(
    plot_data$group,
    levels = plot_data$group
  )

  test_label <- paste0(
    "Genome permutation test\n",
    "Observed = ",
    round(test_table$observed_percent, 1),
    "%; expected = ",
    round(test_table$expected_percent, 1),
    "% ± ",
    round(100 * test_table$expected_sd / n_total, 1),
    "%\n",
    "Enrichment = ",
    round(test_table$fold_enrichment, 2),
    "x; z = ",
    round(test_table$z_score, 2),
    "; ",
    format_pvalue(
      test_table$empirical_pvalue,
      test_table$n_permutations
    )
  )

  pie <- ggplot(
    plot_data,
    aes(
      x = "",
      y = percent,
      fill = group
    )
  ) +
    geom_bar(
      width = 1,
      stat = "identity",
      colour = "white",
      linewidth = 0.5
    ) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 5,
      colour = "black"
    ) +
    scale_fill_manual(
      values = get_tf_colours(TF_name)
    ) +
    ggtitle(
      paste0(TF_name, " peaks overlapping TRIM28")
    ) +
    theme_void() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 11)
    )

  label <- cowplot::ggdraw() +
    cowplot::draw_label(
      test_label,
      size = 10,
      hjust = 0.5
    )

  cowplot::plot_grid(
    pie,
    label,
    ncol = 1,
    rel_heights = c(1, 0.28)
  )
}

annotate_regions <- function(regions, txdb, gene2symbol) {
  if (length(regions) == 0) {
    return(data.frame(
      annotation = character(),
      nearest_gene = character(),
      distanceToTSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      regions,
      tssRegion = c(-2000, 0),
      TxDb = txdb,
      level = "gene",
      overlap = "all"
    )
  )

  annotation$geneId <- gene2symbol[
    annotation$geneId,
    "gene_name"
  ]

  data.frame(
    annotation = annotation$annotation,
    nearest_gene = annotation$geneId,
    distanceToTSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

write_overlap_table <- function(
  TF_name,
  common_regions,
  test_table,
  output_file
) {
  annotation_table <- annotate_regions(
    regions = common_regions,
    txdb = txdb,
    gene2symbol = gene2symbol
  )

  output_table <- data.frame(
    TF = TF_name,
    chromosome = as.character(GenomicRanges::seqnames(common_regions)),
    start = GenomicRanges::start(common_regions),
    end = GenomicRanges::end(common_regions),
    region = make_region_id(common_regions),
    annotation_table,
    observed_overlap = test_table$observed_overlap,
    TF_total = test_table$TF_total,
    expected_overlap = test_table$expected_overlap,
    expected_sd = test_table$expected_sd,
    fold_enrichment = test_table$fold_enrichment,
    z_score = test_table$z_score,
    empirical_pvalue = test_table$empirical_pvalue,
    n_permutations = test_table$n_permutations,
    background = test_table$background,
    stringsAsFactors = FALSE
  )

  write.table(
    output_table,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

###########################################
# Run analysis
###########################################

genome_background <- prepare_genome_background(
  genome_gtf = genome_gtf,
  TF_list = TF_list,
  TRIM28_peaks = TRIM28_peaks
)

TRIM28_peaks <- keep_standard_chromosomes(
  TRIM28_peaks,
  genome_background
)

TF_list <- lapply(
  TF_list,
  keep_standard_chromosomes,
  genome = genome_background
)

results <- lapply(names(TF_list), function(TF_name) {
  log_message("Processing ", TF_name)

  TF_regions <- TF_list[[TF_name]]

  hits <- GenomicRanges::findOverlaps(
    TF_regions,
    TRIM28_peaks,
    ignore.strand = TRUE
  )

  common_TF <- TF_regions[
    unique(S4Vectors::queryHits(hits))
  ]

  test_table <- run_genome_permutation_test(
    TF_regions = TF_regions,
    TRIM28_regions = TRIM28_peaks,
    genome = genome_background,
    TF_name = TF_name,
    n_permutations = n_permutations,
    random_seed = random_seed
  )

  output_file <- file.path(
    snakemake@output[["overlaps"]],
    paste0("TRIM28_", TF_name, ".tsv")
  )

  write_overlap_table(
    TF_name = TF_name,
    common_regions = common_TF,
    test_table = test_table,
    output_file = output_file
  )

  list(
    plot = make_pie_plot(
      TF_name = TF_name,
      test_table = test_table
    ),
    test = test_table
  )
})

names(results) <- names(TF_list)

plots <- lapply(results, `[[`, "plot")

test_table <- do.call(
  rbind,
  lapply(results, `[[`, "test")
)

write.table(
  test_table,
  file = file.path(
    snakemake@output[["overlaps"]],
    "TRIM28_TF_genome_background_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

figure <- cowplot::plot_grid(
  plotlist = plots,
  labels = names(plots),
  ncol = 1,
  align = "v"
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 15,
  base_height = 36,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 15,
  base_height = 36,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")