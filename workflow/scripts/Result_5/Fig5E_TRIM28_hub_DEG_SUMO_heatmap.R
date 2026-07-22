source(".Rprofile")

log_message <- function(...) message("[INFO] ", ...)
log_error <- function(...) stop("[ERROR] ", ...)

log_message("Starting ATAC / H3K27ac / SUMO profile plots around TRIM28-TF regions near DOWN genes")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("S4Vectors")
  library("IRanges")
  library("rtracklayer")
  library("txdbmaker")
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

gene_window <- if ("gene_window" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["gene_window"]])
} else {
  50000
}

distance_to_feature <- if ("distance_to_feature" %in% names(snakemake@params)) {
  as.numeric(snakemake@params[["distance_to_feature"]])
} else {
  1000
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

conditions <- unique(
  gsub("_R.*", "", samplesheet$sample)
)

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
  stringsAsFactors = FALSE,
  check.names = FALSE
)

genome_file <- snakemake@input[["genome"]]

if (!all(c("Gene", "DEG") %in% colnames(DEG_table))) {
  log_error("DEG table must contain 'Gene' and 'DEG' columns")
}

###########################################
# Find bigWig files
###########################################

find_bigwig_by_condition <- function(
  folder,
  condition,
  mark_name = NULL
) {
  matched_files <- list.files(
    folder,
    pattern = condition,
    full.names = TRUE
  )

  matched_files <- matched_files[
    grepl("\\.bw$|\\.bigWig$", matched_files, ignore.case = TRUE)
  ]

  if (length(matched_files) > 1 && !is.null(mark_name)) {
    matched_files <- matched_files[
      grepl(mark_name, basename(matched_files), ignore.case = TRUE)
    ]
  }

  if (length(matched_files) != 1) {
    log_error(
      "Expected exactly one ",
      ifelse(is.null(mark_name), "bigWig", mark_name),
      " file for ",
      condition,
      ", found: ",
      paste(basename(matched_files), collapse = ", ")
    )
  }

  matched_files
}

ATAC_bw_files <- sapply(
  conditions,
  function(condition) {
    find_bigwig_by_condition(
      folder = ATAC_bw_folder,
      condition = condition,
      mark_name = "ATAC"
    )
  }
)

names(ATAC_bw_files) <- conditions

H3K27ac_bw_files <- sapply(
  conditions,
  function(condition) {
    find_bigwig_by_condition(
      folder = H3K27ac_bw_folder,
      condition = condition,
      mark_name = "H3K27ac"
    )
  }
)

names(H3K27ac_bw_files) <- conditions

SUMO1_bw_files <- sapply(
  conditions,
  function(condition) {
    find_bigwig_by_condition(
      folder = ChIP_bw_folder,
      condition = condition,
      mark_name = "SUMO1"
    )
  }
)

names(SUMO1_bw_files) <- conditions

SUMO2_bw_files <- sapply(
  conditions,
  function(condition) {
    find_bigwig_by_condition(
      folder = ChIP_bw_folder,
      condition = condition,
      mark_name = "SUMO2"
    )
  }
)

names(SUMO2_bw_files) <- conditions

bw_files <- list(
  ATAC = ATAC_bw_files,
  H3K27ac = H3K27ac_bw_files,
  SUMO1 = SUMO1_bw_files,
  SUMO2 = SUMO2_bw_files
)

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

get_TRIM28_regions_near_DOWN_genes <- function(
  TRIM28_peaks,
  TRIM28_table,
  gene_tss,
  DEG_table,
  gene_window
) {
  down_genes <- unique(
    DEG_table$gene[
      DEG_table$DEG == "DOWN"
    ]
  )

  down_tss <- gene_tss[
    gene_tss$gene_name %in% down_genes
  ]

  if (length(down_tss) == 0) {
    log_error("No DOWN gene TSS found in genome annotation")
  }

  expanded_down_tss <- GenomicRanges::resize(
    down_tss,
    width = 1 + 2 * gene_window,
    fix = "center"
  )

  expanded_down_tss <- GenomicRanges::trim(
    expanded_down_tss
  )

  is_near_down_gene <- count_overlap_binary(
    query = TRIM28_peaks,
    subject = expanded_down_tss
  )

  TRIM28_subset <- TRIM28_peaks[
    is_near_down_gene
  ]

  TRIM28_table_subset <- TRIM28_table[
    is_near_down_gene,
    ,
    drop = FALSE
  ]

  log_message("DOWN genes: ", length(down_genes))
  log_message("TRIM28 regions near DOWN genes: ", length(TRIM28_subset))

  if (length(TRIM28_subset) == 0) {
    log_error("No TRIM28 regions found near DOWN genes")
  }

  list(
    regions = TRIM28_subset,
    table = TRIM28_table_subset
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
  classes,
  conditions,
  distance_to_feature,
  bin_size
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
        position = position,
        signal = profile,
        n_regions = sum(class_index, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(
    rbind,
    all_profiles
  )
}

make_condition_linetypes <- function(conditions) {
  linetypes <- ifelse(
    grepl("ko|cko|mut", conditions, ignore.case = TRUE),
    "dashed",
    "solid"
  )

  if (length(unique(linetypes)) == 1 && length(conditions) >= 2) {
    linetypes <- c(
      "solid",
      rep("dashed", length(conditions) - 1)
    )
  }

  names(linetypes) <- conditions
  linetypes
}

make_condition_colours <- function(conditions) {
  colours <- ifelse(
    grepl("ko|cko|mut", conditions, ignore.case = TRUE),
    "#7D1128",
    "#333333"
  )

  if (length(unique(colours)) == 1 && length(conditions) >= 2) {
    colours <- c(
      "#333333",
      rep("#7D1128", length(conditions) - 1)
    )
  }

  names(colours) <- conditions
  colours
}

plot_mark_profiles <- function(
  profile_table,
  mark_name,
  conditions
) {
  plot_table <- profile_table[
    profile_table$mark == mark_name,
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

  plot_table$facet_label <- factor(
    plot_table$facet_label,
    levels = class_labels$facet_label[
      match(
        levels(factor(class_labels$TRIM28_TF_class)),
        class_labels$TRIM28_TF_class
      )
    ]
  )

  ggplot(
    plot_table,
    aes(
      x = position,
      y = signal,
      colour = condition,
      linetype = condition
    )
  ) +
    geom_line(
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      colour = "#555555",
      linewidth = 0.3
    ) +
    facet_wrap(
      ~facet_label,
      ncol = 5,
      scales = "free_y"
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
        " coverage at TRIM28-TF regions near DOWN genes"
      )
    ) +
    theme_light(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
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

DEG_table <- prepare_DEG_table(
  DEG_table = DEG_table
)

TRIM28_table <- annotate_TRIM28_regions(
  TRIM28_peaks = TRIM28_peaks,
  TF_list = TF_list
)

selected <- get_TRIM28_regions_near_DOWN_genes(
  TRIM28_peaks = TRIM28_peaks,
  TRIM28_table = TRIM28_table,
  gene_tss = gene_tss,
  DEG_table = DEG_table,
  gene_window = gene_window
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

log_message("Selected regions by TRIM28-TF class:")

for (class_name in levels(classes)) {
  log_message(
    class_name,
    ": ",
    sum(classes == class_name)
  )
}

###########################################
# Compute profiles
###########################################

profile_tables <- lapply(
  names(bw_files),
  function(mark_name) {
    make_profile_table_for_mark(
      mark_name = mark_name,
      bigwig_files = bw_files[[mark_name]],
      regions = selected_regions,
      classes = classes,
      conditions = conditions,
      distance_to_feature = distance_to_feature,
      bin_size = bin_size
    )
  }
)

profile_table <- do.call(
  rbind,
  profile_tables
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

###########################################
# Save figures
###########################################

mark_plots <- lapply(
  names(bw_files),
  function(mark_name) {
    plot_mark_profiles(
      profile_table = profile_table,
      mark_name = mark_name,
      conditions = conditions
    )
  }
)

names(mark_plots) <- names(bw_files)

pdf(
  file = snakemake@output[["pdf"]],
  width = 22 / 2.54,
  height = 9 / 2.54
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
  base_height = 9 * length(mark_plots),
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")