source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting gene-centric maximal TRIM28-TF class DEG analysis")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("IRanges")
  library("rtracklayer")
  library("txdbmaker")
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

DEG_table <- read.csv(
  snakemake@input[["DEG"]],
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

genome_file <- snakemake@input[["genome"]]

if (!all(c("Gene", "DEG") %in% colnames(DEG_table))) {
  log_error("DEG table must contain 'Gene' and 'DEG' columns")
}

has_expressed_genes <- "expressed_genes" %in% names(snakemake@input)

if (has_expressed_genes) {
  expressed_genes_table <- read.csv(
    snakemake@input[["expressed_genes"]],
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Gene" %in% colnames(expressed_genes_table)) {
    log_error("expressed_genes table must contain a 'Gene' column")
  }
}

###########################################
# Parameters
###########################################

distance_windows <- if ("distance_windows" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["distance_windows"]])
} else {
  c(50000)
}

distance_windows <- distance_windows[
  !is.na(distance_windows) &
    distance_windows > 0
]

if (length(distance_windows) == 0) {
  log_error("distance_windows must contain at least one positive numeric value")
}

###########################################
# Genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  genome_file
)

txdb <- txdbmaker::makeTxDbFromGFF(
  genome_file
)

genes_gr <- GenomicFeatures::genes(
  txdb
)

gene2symbol <- GenomicRanges::mcols(genome_gtf)[
  ,
  c("gene_id", "gene_name")
]

gene2symbol <- unique(gene2symbol)

gene2symbol <- gene2symbol[
  !is.na(gene2symbol$gene_id) &
    !is.na(gene2symbol$gene_name),
  ,
  drop = FALSE
]

rownames(gene2symbol) <- gene2symbol$gene_id

genes_gr$gene_id <- names(genes_gr)

genes_gr$gene_name <- gene2symbol[
  genes_gr$gene_id,
  "gene_name"
]

genes_gr$gene_name[
  is.na(genes_gr$gene_name)
] <- genes_gr$gene_id[
  is.na(genes_gr$gene_name)
]

gene_tss <- GenomicRanges::resize(
  genes_gr,
  width = 1,
  fix = "start"
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

count_overlap_binary <- function(query, subject) {
  GenomicRanges::countOverlaps(
    query,
    subject,
    ignore.strand = TRUE
  ) > 0
}

prepare_DEG_table <- function(DEG_table) {
  DEG_table <- DEG_table[
    !is.na(DEG_table$Gene) &
      DEG_table$Gene != "",
    ,
    drop = FALSE
  ]

  DEG_table$Gene <- as.character(DEG_table$Gene)
  DEG_table$DEG <- toupper(as.character(DEG_table$DEG))

  DEG_table <- DEG_table[
    DEG_table$DEG %in% c("UP", "DOWN"),
    ,
    drop = FALSE
  ]

  DEG_table <- DEG_table[
    !duplicated(DEG_table$Gene),
    c("Gene", "DEG"),
    drop = FALSE
  ]

  colnames(DEG_table) <- c("gene", "DEG")

  DEG_table
}

prepare_gene_universe <- function(
  genes_gr,
  DEG_table,
  has_expressed_genes,
  expressed_genes_table = NULL
) {
  if (has_expressed_genes) {
    log_message("Using expressed genes as statistical background")

    gene_universe <- unique(
      as.character(expressed_genes_table$Gene)
    )
  } else {
    log_message("No expressed_genes input provided; using all annotated genes as background")

    gene_universe <- unique(
      genes_gr$gene_name[
        !is.na(genes_gr$gene_name)
      ]
    )
  }

  gene_universe <- gene_universe[
    !is.na(gene_universe) &
      gene_universe != ""
  ]

  missing_DEGs <- setdiff(
    DEG_table$gene,
    gene_universe
  )

  if (length(missing_DEGs) > 0) {
    log_message(
      "Warning: ",
      length(missing_DEGs),
      " DEG genes are not present in the background gene universe"
    )
  }

  gene_universe
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

  region_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(TRIM28_peaks)),
    start = GenomicRanges::start(TRIM28_peaks),
    end = GenomicRanges::end(TRIM28_peaks),
    region = make_region_id(TRIM28_peaks),
    TF_overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = TF_combination,
    stringsAsFactors = FALSE
  )

  region_table$TRIM28_TF_class <- ifelse(
    region_table$n_TFs == 0,
    "TRIM28 only",
    ifelse(
      region_table$n_TFs == 1,
      "TRIM28 + 1 TF",
      ifelse(
        region_table$n_TFs == 2,
        "TRIM28 + 2 TFs",
        ifelse(
          region_table$n_TFs == 3,
          "TRIM28 + 3 TFs",
          "TRIM28 + 4 TFs"
        )
      )
    )
  )

  region_table$TRIM28_TF_class <- factor(
    region_table$TRIM28_TF_class,
    levels = c(
      "TRIM28 only",
      "TRIM28 + 1 TF",
      "TRIM28 + 2 TFs",
      "TRIM28 + 3 TFs",
      "TRIM28 + 4 TFs"
    )
  )

  region_table
}

get_gene_max_TRIM28_class <- function(
  TRIM28_peaks,
  TRIM28_table,
  gene_tss,
  window,
  gene_universe
) {
  log_message(
    "Assigning each gene to maximal nearby TRIM28-TF class within ±",
    window / 1000,
    " kb"
  )

  expanded_regions <- GenomicRanges::resize(
    TRIM28_peaks,
    width = GenomicRanges::width(TRIM28_peaks) + 2 * window,
    fix = "center"
  )

  expanded_regions <- GenomicRanges::trim(
    expanded_regions
  )

  hits <- GenomicRanges::findOverlaps(
    gene_tss,
    expanded_regions,
    ignore.strand = TRUE
  )

  links <- data.frame(
    gene = gene_tss$gene_name[
      S4Vectors::queryHits(hits)
    ],
    region = TRIM28_table$region[
      S4Vectors::subjectHits(hits)
    ],
    n_TFs = TRIM28_table$n_TFs[
      S4Vectors::subjectHits(hits)
    ],
    TF_combination = TRIM28_table$TF_combination[
      S4Vectors::subjectHits(hits)
    ],
    TRIM28_TF_class = as.character(
      TRIM28_table$TRIM28_TF_class[
        S4Vectors::subjectHits(hits)
      ]
    ),
    distance_window = window,
    stringsAsFactors = FALSE
  )

  links <- links[
    links$gene %in% gene_universe,
    ,
    drop = FALSE
  ]

  if (nrow(links) == 0) {
    log_error("No gene-TRIM28 links found for window ±", window / 1000, " kb")
  }

  gene_max <- aggregate(
    n_TFs ~ gene,
    data = links,
    FUN = max
  )

  colnames(gene_max) <- c(
    "gene",
    "max_n_TFs"
  )

  gene_max$max_TRIM28_TF_class <- ifelse(
    gene_max$max_n_TFs == 0,
    "TRIM28 only",
    ifelse(
      gene_max$max_n_TFs == 1,
      "TRIM28 + 1 TF",
      ifelse(
        gene_max$max_n_TFs == 2,
        "TRIM28 + 2 TFs",
        ifelse(
          gene_max$max_n_TFs == 3,
          "TRIM28 + 3 TFs",
          "TRIM28 + 4 TFs"
        )
      )
    )
  )

  gene_table <- data.frame(
    gene = gene_universe,
    distance_window = window,
    stringsAsFactors = FALSE
  )

  gene_table <- merge(
    gene_table,
    gene_max,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  gene_table$max_n_TFs[
    is.na(gene_table$max_n_TFs)
  ] <- -1

  gene_table$max_TRIM28_TF_class[
    is.na(gene_table$max_TRIM28_TF_class)
  ] <- "No nearby TRIM28"

  gene_table$max_TRIM28_TF_class <- factor(
    gene_table$max_TRIM28_TF_class,
    levels = c(
      "No nearby TRIM28",
      "TRIM28 only",
      "TRIM28 + 1 TF",
      "TRIM28 + 2 TFs",
      "TRIM28 + 3 TFs",
      "TRIM28 + 4 TFs"
    )
  )

  list(
    gene_table = gene_table,
    gene_region_links = links
  )
}

add_DEG_status <- function(
  gene_table,
  DEG_table
) {
  gene_table <- merge(
    gene_table,
    DEG_table,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  gene_table$DEG[
    is.na(gene_table$DEG)
  ] <- "-"

  gene_table$DEG_category <- ifelse(
    gene_table$DEG == "DOWN",
    "DOWN genes",
    ifelse(
      gene_table$DEG == "UP",
      "UP genes",
      "No change"
    )
  )

  gene_table$DEG_category <- factor(
    gene_table$DEG_category,
    levels = c(
      "DOWN genes",
      "UP genes",
      "No change"
    )
  )

  gene_table
}

summarise_gene_max_classes <- function(
  gene_table,
  window
) {
  class_order <- c(
    "No nearby TRIM28",
    "TRIM28 only",
    "TRIM28 + 1 TF",
    "TRIM28 + 2 TFs",
    "TRIM28 + 3 TFs",
    "TRIM28 + 4 TFs"
  )

  summary_table <- as.data.frame(
    table(
      gene_table$DEG_category,
      gene_table$max_TRIM28_TF_class
    ),
    stringsAsFactors = FALSE
  )

  colnames(summary_table) <- c(
    "DEG_category",
    "max_TRIM28_TF_class",
    "n_genes"
  )

  total_by_category <- aggregate(
    n_genes ~ DEG_category,
    data = summary_table,
    FUN = sum
  )

  colnames(total_by_category) <- c(
    "DEG_category",
    "n_total_genes"
  )

  summary_table <- merge(
    summary_table,
    total_by_category,
    by = "DEG_category",
    all.x = TRUE,
    sort = FALSE
  )

  summary_table$pct_genes <- 100 *
    summary_table$n_genes /
    summary_table$n_total_genes

  summary_table$distance_window <- window

  summary_table$DEG_category <- factor(
    summary_table$DEG_category,
    levels = c(
      "DOWN genes",
      "UP genes",
      "No change"
    )
  )

  summary_table$max_TRIM28_TF_class <- factor(
    summary_table$max_TRIM28_TF_class,
    levels = class_order
  )

  summary_table$label <- paste0(
    round(summary_table$pct_genes, 1),
    "%"
    # "%\n(",
    # summary_table$n_genes,
    # "/",
    # summary_table$n_total_genes,
    # ")"
  )

  summary_table
}

run_fisher_max_class <- function(
  gene_table,
  window
) {
  tests <- list(
    "DOWN genes enriched in TRIM28 + >=3 TFs" = list(
      response = gene_table$DEG_category == "DOWN genes",
      exposure = gene_table$max_n_TFs >= 3,
      response_label = "DOWN genes",
      exposure_label = "TRIM28 + >=3 TFs"
    ),
    "DOWN genes enriched in TRIM28 + 4 TFs" = list(
      response = gene_table$DEG_category == "DOWN genes",
      exposure = gene_table$max_n_TFs >= 4,
      response_label = "DOWN genes",
      exposure_label = "TRIM28 + 4 TFs"
    ),
    "UP genes enriched in TRIM28 + >=3 TFs" = list(
      response = gene_table$DEG_category == "UP genes",
      exposure = gene_table$max_n_TFs >= 3,
      response_label = "UP genes",
      exposure_label = "TRIM28 + >=3 TFs"
    ),
    "UP genes enriched in TRIM28 + 4 TFs" = list(
      response = gene_table$DEG_category == "UP genes",
      exposure = gene_table$max_n_TFs >= 4,
      response_label = "UP genes",
      exposure_label = "TRIM28 + 4 TFs"
    )
  )

  fisher_tables <- lapply(
    names(tests),
    function(test_name) {
      response <- tests[[test_name]]$response
      exposure <- tests[[test_name]]$exposure

      contingency_table <- matrix(
        c(
          sum(response & exposure),
          sum(response & !exposure),
          sum(!response & exposure),
          sum(!response & !exposure)
        ),
        nrow = 2,
        byrow = TRUE
      )

      fisher_result <- stats::fisher.test(
        contingency_table,
        alternative = "greater"
      )

      data.frame(
        distance_window = window,
        test = test_name,
        response = tests[[test_name]]$response_label,
        exposure = tests[[test_name]]$exposure_label,
        n_response_with_exposure = contingency_table[1, 1],
        n_response_total = sum(contingency_table[1, ]),
        pct_response_with_exposure = 100 * contingency_table[1, 1] / sum(contingency_table[1, ]),
        n_background_with_exposure = contingency_table[2, 1],
        n_background_total = sum(contingency_table[2, ]),
        pct_background_with_exposure = 100 * contingency_table[2, 1] / sum(contingency_table[2, ]),
        odds_ratio = unname(fisher_result$estimate),
        p_value = fisher_result$p.value,
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    fisher_tables
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

make_fisher_caption <- function(
  fisher_summary_all
) {
  fisher_summary_all$window_label <- paste0(
    "±",
    fisher_summary_all$distance_window / 1000,
    " kb"
  )

  fisher_summary_all$line <- paste0(
    fisher_summary_all$window_label,
    " | ",
    fisher_summary_all$response,
    " enriched in ",
    fisher_summary_all$exposure,
    ": OR = ",
    round(fisher_summary_all$odds_ratio, 2),
    ", ",
    vapply(
      fisher_summary_all$p_value,
      format_pvalue,
      character(1)
    )
  )

  paste(
    fisher_summary_all$line,
    collapse = "\n"
  )
}

plot_gene_max_class_stacked <- function(
  gene_summary_all,
  fisher_summary_all
) {
  gene_summary_all$window_label <- paste0(
    "±",
    gene_summary_all$distance_window / 1000,
    " kb"
  )

  gene_summary_all$max_TRIM28_TF_class <- factor(
    gene_summary_all$max_TRIM28_TF_class,
    levels = c(
      "No nearby TRIM28",
      "TRIM28 only",
      "TRIM28 + 1 TF",
      "TRIM28 + 2 TFs",
      "TRIM28 + 3 TFs",
      "TRIM28 + 4 TFs"
    )
  )

  main_plot <- ggplot(
    gene_summary_all,
    aes(
      x = DEG_category,
      y = pct_genes,
      fill = max_TRIM28_TF_class
    )
  ) +
    geom_col(
      colour = "#333333",
      linewidth = 0.25,
      width = 0.75
    ) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4,
      lineheight = 0.9
    ) +
    # facet_wrap(
    #   ~window_label
    # ) +
    scale_fill_manual(
      values = c(
        "No nearby TRIM28" = "grey90",
        "TRIM28 only" = "#C9C9C9",
        "TRIM28 + 1 TF" = "#B7A3CC",
        "TRIM28 + 2 TFs" = "#9C78B8",
        "TRIM28 + 3 TFs" = "#7D4EA3",
        "TRIM28 + 4 TFs" = "#4B1D73"
      )
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%")
      # limits = c(0, 100),
      # expand = expansion(mult = c(0, 0))
    ) +
    xlab(NULL) +
    ylab("Gene distribution") +
    ggtitle(
      "Distribution of maximal nearby TRIM28-TF class by DEG category"
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        size = 12
      ),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.text = element_text(
        face = "bold",
        size = 12
      )
    )

  fisher_text <- make_fisher_caption(
    fisher_summary_all = fisher_summary_all
  )

  fisher_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 1,
      label = fisher_text,
      hjust = 0,
      vjust = 1,
      size = 3.6
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 12)

  cowplot::plot_grid(
    main_plot,
    fisher_plot,
    ncol = 1,
    rel_heights = c(1, 0.22)
  )
}

###########################################
# Run analysis
###########################################

DEG_table <- prepare_DEG_table(
  DEG_table = DEG_table
)

gene_universe <- prepare_gene_universe(
  genes_gr = genes_gr,
  DEG_table = DEG_table,
  has_expressed_genes = has_expressed_genes,
  expressed_genes_table = if (has_expressed_genes) expressed_genes_table else NULL
)

TRIM28_table <- annotate_TRIM28_regions(
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

all_gene_tables <- list()
all_gene_region_links <- list()
all_gene_summaries <- list()
all_fisher_summaries <- list()

for (window in distance_windows) {
  log_message("Running analysis with window ±", window / 1000, " kb")

  max_class_results <- get_gene_max_TRIM28_class(
    TRIM28_peaks = TRIM28_peaks,
    TRIM28_table = TRIM28_table,
    gene_tss = gene_tss,
    window = window,
    gene_universe = gene_universe
  )

  gene_table <- add_DEG_status(
    gene_table = max_class_results$gene_table,
    DEG_table = DEG_table
  )

  gene_summary <- summarise_gene_max_classes(
    gene_table = gene_table,
    window = window
  )

  fisher_summary <- run_fisher_max_class(
    gene_table = gene_table,
    window = window
  )

  all_gene_tables[[as.character(window)]] <- gene_table
  all_gene_region_links[[as.character(window)]] <- max_class_results$gene_region_links
  all_gene_summaries[[as.character(window)]] <- gene_summary
  all_fisher_summaries[[as.character(window)]] <- fisher_summary
}

gene_table_all <- do.call(
  rbind,
  all_gene_tables
)

gene_region_links_all <- do.call(
  rbind,
  all_gene_region_links
)

gene_summary_all <- do.call(
  rbind,
  all_gene_summaries
)

fisher_summary_all <- do.call(
  rbind,
  all_fisher_summaries
)

###########################################
# Save outputs
###########################################

write.table(
  TRIM28_table,
  file = snakemake@output[["regions"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  gene_table_all,
  file = snakemake@output[["genes"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  gene_region_links_all,
  file = snakemake@output[["gene_region_links"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  gene_summary_all,
  file = snakemake@output[["gene_summary"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  fisher_summary_all,
  file = snakemake@output[["fisher_summary"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

###########################################
# Save figure
###########################################

figure <- plot_gene_max_class_stacked(
  gene_summary_all = gene_summary_all,
  fisher_summary_all = fisher_summary_all
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 18,
  base_height = 14,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 18,
  base_height = 14,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")