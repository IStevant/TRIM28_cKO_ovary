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

log_message("Starting TFBS motif enrichment analysis on ATAC-Sertoli-overlapping regions")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomeInfoDb")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("Biostrings")
  library("BiocParallel")
  library("JASPAR2020")
  library("JASPAR2024")
  library("TFBSTools")
  library("monaLisa")
  library("SummarizedExperiment")
  library("ComplexHeatmap")
  library("circlize")
  library("motifStack")
  library("universalmotif")
  library("matrixStats")
  library("cowplot")
  library("grid")
  library("ggplot2")
  library("stringr")
  library("ade4")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 12))
genome_version <- "mm10"
background <- "genome"

TFBS_log10FDR <- as.numeric(
  snakemake@params[["TFBS_log10FDR"]]
)

if (is.null(TFBS_log10FDR) || length(TFBS_log10FDR) == 0) {
  log_error("Parameter 'TFBS_log10FDR' is missing from the Snakemake rule")
}

TFBS_FDR <- 10^(-TFBS_log10FDR)

log_message("Using -log10(FDR) threshold: ", TFBS_log10FDR)
log_message("Equivalent FDR threshold: ", TFBS_FDR)

BiocParallel::register(
  BiocParallel::MulticoreParam(n_cores)
)

###########################################
# Load input data
###########################################

log_message("Loading ATAC / ATAC-Sertoli overlap table")

overlap_table <- read.csv(
  file = snakemake@input[["sig_DARs"]],
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

required_columns <- c(
  "region",
  "cluster",
  "overlap_ATAC_Sertoli"
)

missing_columns <- setdiff(
  required_columns,
  colnames(overlap_table)
)

if (length(missing_columns) > 0) {
  log_error(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

overlap_table <- overlap_table[
  overlap_table$overlap_ATAC_Sertoli == "Yes" &
    !is.na(overlap_table$region) &
    !is.na(overlap_table$cluster),
  ,
  drop = FALSE
]

if (nrow(overlap_table) == 0) {
  log_error("No ATAC regions overlapping ATAC-Sertoli peaks were found")
}

clusters <- sort(
  unique(as.character(overlap_table$cluster))
)

log_message("Clusters retained: ", paste(clusters, collapse = ", "))

res_table <- snakemake@output[["res_table"]]

###########################################
# Load JASPAR database
###########################################

log_message("Loading JASPAR 2024 database")

JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024()@db

###########################################
# Functions
###########################################

get_genome_object <- function(genome_version) {
  if (genome_version == "mm10") {
    return(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  }

  log_error("Unsupported genome version: ", genome_version)
}

get_midpoint_gr <- function(regions) {
  mids <- floor(
    (GenomicRanges::start(regions) + GenomicRanges::end(regions)) / 2
  )

  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(regions),
    ranges = IRanges::IRanges(
      start = mids,
      end = mids
    ),
    strand = "*"
  )
}

get_enriched_TFs <- function(
  regions,
  cluster_name,
  TFBS_log10FDR
) {
  log_message("Running motif enrichment for cluster: ", cluster_name)

  pwms <- TFBSTools::getMatrixSet(
    JASPAR,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  regions <- GenomeInfoDb::keepSeqlevels(
    regions,
    paste0("chr", c(1:19, "X")),
    pruning.mode = "coarse"
  )

  if (length(regions) == 0) {
    log_message("No regions retained after seqlevel filtering for cluster: ", cluster_name)
    return(NULL)
  }

  summit_regions <- GenomicRanges::resize(
    get_midpoint_gr(regions),
    width = 500,
    fix = "center"
  )

  genome <- get_genome_object(genome_version)

  sequences <- Biostrings::getSeq(
    genome,
    summit_regions
  )

  enrichment <- monaLisa::calcBinnedMotifEnrR(
    seqs = sequences,
    pwmL = pwms,
    background = "genome",
    genome = genome,
    genome.regions = NULL,
    genome.oversample = 2,
    BPPARAM = BiocParallel::MulticoreParam(n_cores)
  )

  if (length(as.vector(SummarizedExperiment::rowData(enrichment)$motif.name)) <= 1) {
    return(NULL)
  }

  selected_motifs <- apply(
    SummarizedExperiment::assay(enrichment, "negLog10Padj"),
    1,
    function(x) {
      max(abs(x), na.rm = TRUE)
    }
  ) >= TFBS_log10FDR

  enrichment_selected <- enrichment[selected_motifs, ]

  if (nrow(enrichment_selected) == 0) {
    return(NULL)
  }

  enrichment_selected
}

extract_tf_prefix <- function(gene) {
  prefix <- sub("([a-zA-Z]+).*", "\\1", gene)

  if (prefix == "NR") {
    return(gene)
  }

  if (grepl("^HOX", prefix)) {
    return(stringr::str_sub(prefix, end = -2))
  }

  prefix
}

collapse_tf_names <- function(tf_vector) {
  grouped_genes <- split(
    tf_vector,
    sapply(tf_vector, extract_tf_prefix)
  )

  new_tf_vector <- sapply(grouped_genes, function(gene_group) {
    common_prefix <- extract_tf_prefix(gene_group[1])

    if (common_prefix == "HOX") {
      return("HOXs")
    }

    suffixes <- gsub(
      paste0("^", common_prefix),
      "",
      gene_group
    )

    numeric_suffixes <- sort(
      as.numeric(suffixes[grepl("^\\d+$", suffixes)]),
      na.last = TRUE
    )

    non_numeric_suffixes <- sort(
      suffixes[!grepl("^\\d+$", suffixes)]
    )

    paste0(
      common_prefix,
      paste(
        c(non_numeric_suffixes, numeric_suffixes),
        collapse = "/"
      )
    )
  })

  paste(
    new_tf_vector[order(new_tf_vector)],
    collapse = ";"
  )
}

prepare_motif_matrix <- function(
  enrichment_selected,
  TFBS_FDR
) {
  TF_enrichment <- SummarizedExperiment::assay(
    enrichment_selected,
    "log2enr"
  )

  TF_FDR <- 10^(
    -SummarizedExperiment::assay(
      enrichment_selected,
      "negLog10Padj"
    )
  )

  colnames(TF_FDR) <- paste0(
    "FDR_",
    colnames(TF_FDR)
  )

  enrichment_for_clustering <- cbind(
    TF_enrichment,
    TF_FDR
  )

  hcl <- stats::hclust(
    stats::dist(enrichment_for_clustering),
    method = "ward.D2"
  )

  clustering <- stats::cutree(
    hcl,
    k = 1
  )

  motif_signatures <- lapply(seq_len(1), function(cluster_id) {
    TFs <- names(clustering[clustering == cluster_id])

    if (length(TFs) > 1) {
      matrices <- SummarizedExperiment::rowData(enrichment_selected)$motif.pfm[TFs]

      pfms <- universalmotif::convert_motifs(
        matrices,
        class = "motifStack-pfm"
      )

      TF_names <- make.unique(
        unlist(lapply(pfms, function(x) x$name)),
        sep = "XXX"
      )

      TF_names <- gsub("::", "YYY", TF_names)
      TF_names <- gsub("-", "ZZZ", TF_names)

      pfms <- lapply(seq_along(pfms), function(i) {
        pfms[[i]]$name <- TF_names[[i]]
        pfms[[i]]
      })

      names(pfms) <- TF_names

      hc <- motifStack::clusterMotifs(pfms)
      phylog <- ade4::hclust2phylog(hc)
      pfms <- pfms[names(phylog$leaves)]

      motif_signature <- motifStack::motifSignature(
        pfms,
        phylog,
        cutoffPval = 0.0001,
        min.freq = 1
      )

      motifStack::signatures(motif_signature)
    } else {
      universalmotif::convert_motifs(
        SummarizedExperiment::rowData(enrichment_selected)$motif.pfm[TFs],
        class = "motifStack-pcm"
      )
    }
  })

  motif_signature_enrichment <- lapply(motif_signatures, function(cluster) {
    lapply(cluster, function(motif) {
      TF_names <- motif@name

      lapply(TF_names, function(TF) {
        tf_vector <- toupper(
          unlist(strsplit(TF, ";"))
        )

        enrichment <- enrichment_for_clustering

        enrichment_names <- make.unique(
          toupper(SummarizedExperiment::rowData(enrichment_selected)$motif.name)
        )

        enrichment_names <- gsub("[.]", "XXX", enrichment_names)
        enrichment_names <- gsub("::", "YYY", enrichment_names)
        enrichment_names <- gsub("-", "ZZZ", enrichment_names)

        rownames(enrichment) <- enrichment_names

        tf_enrichment <- enrichment[
          tf_vector,
          ,
          drop = FALSE
        ]

        if (nrow(tf_enrichment) > 1) {
          tf_names <- collapse_tf_names(tf_vector)

          tf_enrichment <- as.data.frame(
            matrixStats::colMedians(
              as.matrix(tf_enrichment)
            )
          )

          colnames(tf_enrichment) <- tf_names
          tf_enrichment <- as.data.frame(t(tf_enrichment))
          tf_enrichment$mat_names <- TF
        } else {
          tf_enrichment <- as.data.frame(tf_enrichment)
          tf_enrichment$mat_names <- TF
        }

        tf_enrichment
      })
    })
  })

  enrichment <- do.call(
    "rbind",
    do.call(
      "rbind",
      do.call("rbind", motif_signature_enrichment)
    )
  )

  enrichment <- enrichment[
    !duplicated(enrichment),
    ,
    drop = FALSE
  ]

  merged_pfms <- do.call(
    "c",
    do.call(
      "c",
      do.call("c", motif_signatures)
    )
  )

  names(merged_pfms) <- unlist(
    lapply(merged_pfms, function(motif) motif@name)
  )

  motifs <- merged_pfms[
    enrichment$mat_names
  ]

  motifs_pfms <- universalmotif::convert_motifs(
    motifs,
    class = "TFBSTools-PFMatrix"
  )

  max_width <- max(
    unlist(
      lapply(
        motifs_pfms,
        function(x) ncol(x@profileMatrix)
      )
    )
  )

  logo_grobs <- lapply(
    motifs_pfms,
    monaLisa::seqLogoGrob,
    xmax = max_width,
    xjust = "right"
  )

  enrichment_names <- rownames(enrichment)

  enrichment_names <- gsub("XXX", ".", enrichment_names)
  enrichment_names <- gsub("YYY", "::", enrichment_names)
  enrichment_names <- gsub("ZZZ", "-", enrichment_names)

  rownames(enrichment) <- enrichment_names
  names(logo_grobs) <- enrichment_names

  matrix <- enrichment[
    ,
    grep(
      "FDR_|mat_names",
      colnames(enrichment),
      invert = TRUE
    ),
    drop = FALSE
  ]

  FDR <- enrichment[
    ,
    grep("FDR_", colnames(enrichment)),
    drop = FALSE
  ]

  matrix[matrix > 0.5] <- 0.5
  matrix[matrix < 0] <- 0
  matrix[FDR > TFBS_FDR] <- NA

  matrix_no_na <- matrix
  matrix_no_na[is.na(matrix_no_na)] <- 0

  row_order <- order(
    -rowMeans(matrix_no_na, na.rm = TRUE)
  )

  matrix <- matrix[row_order, , drop = FALSE]
  logo_grobs <- logo_grobs[rownames(matrix)]

  list(
    matrix = matrix,
    logo_grobs = logo_grobs,
    FDR = FDR
  )
}

make_heatmap <- function(
  matrix,
  logo_grobs,
  cluster_name
) {
  warm_palette <- grDevices::colorRampPalette(
    c(
      "#fffee8",
      "#ffd9cb",
      "#ffb1ad",
      "#f5808d",
      "#eb2d62"
    )
  )

  colour_function <- circlize::colorRamp2(
    breaks = seq(0, 0.5, length.out = 12),
    colors = warm_palette(12)
  )

  logo_annotation <- ComplexHeatmap::HeatmapAnnotation(
    logo = monaLisa::annoSeqlogo(
      grobL = logo_grobs,
      which = "row",
      space = grid::unit(0.3, "mm"),
      width = grid::unit(1.05, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "Log2 enrichment",
    left_annotation = logo_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(
      fontsize = 7
    ),
    column_title = paste0(
      "TF motif enrichment - ATAC-Sertoli-overlapping regions - cluster ",
      cluster_name
    ),
    column_title_gp = grid::gpar(
      fontsize = 8.4,
      fontface = "bold"
    ),
    width = ncol(matrix) * grid::unit(2, "mm"),
    height = nrow(matrix) * grid::unit(2, "mm"),
    col = colour_function,
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.4
    ),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 7),
      labels_gp = grid::gpar(fontsize = 7),
      direction = "horizontal"
    )
  )

  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "bottom"
    ),
    wrap.grobs = TRUE
  )
}

make_empty_plot <- function(cluster_name) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = paste0(
        "No significant TF motif enrichment\nfor ATAC-Sertoli-overlapping regions\ncluster ",
        cluster_name
      ),
      size = 5,
      fontface = "bold"
    ) +
    ggplot2::ggtitle(
      paste0(
        "TF motif enrichment - cluster ",
        cluster_name
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 11
      )
    )
}

###########################################
# Run enrichment by cluster
###########################################

log_message("Running motif enrichment analysis")

plots <- list()
all_results <- list()

for (cluster_name in clusters) {
  cluster_table <- overlap_table[
    overlap_table$cluster == cluster_name,
    ,
    drop = FALSE
  ]

  cluster_regions <- GenomicRanges::GRanges(
    cluster_table$region
  )

  enrichment <- get_enriched_TFs(
    regions = cluster_regions,
    cluster_name = cluster_name,
    TFBS_log10FDR = TFBS_log10FDR
  )

  if (is.null(enrichment) || nrow(enrichment) == 0) {
    plots[[cluster_name]] <- make_empty_plot(cluster_name)
    next
  }

  heatmap_data <- prepare_motif_matrix(
    enrichment_selected = enrichment,
    TFBS_FDR = TFBS_FDR
  )

  plots[[cluster_name]] <- make_heatmap(
    matrix = heatmap_data$matrix,
    logo_grobs = heatmap_data$logo_grobs,
    cluster_name = cluster_name
  )

  log2_enrichment <- SummarizedExperiment::assay(
    enrichment,
    "log2enr"
  )

  FDR <- 10^(
    -SummarizedExperiment::assay(
      enrichment,
      "negLog10Padj"
    )
  )

  colnames(FDR) <- paste0(
    "FDR_",
    colnames(FDR)
  )

  all_results[[cluster_name]] <- data.frame(
    cluster = cluster_name,
    TF.name = SummarizedExperiment::rowData(enrichment)$motif.name,
    TF.matrix = rownames(log2_enrichment),
    log2_enrichment,
    FDR,
    check.names = FALSE
  )
}

###########################################
# Save motif enrichment table
###########################################

all_results <- do.call(
  rbind,
  all_results
)

if (!is.null(all_results)) {
  write.table(
    all_results,
    file = res_table,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
} else {
  write.table(
    "No significantly enriched TFBS motif found",
    file = res_table,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

###########################################
# Save output figures
###########################################

log_message("Saving PDF")

pdf(
  file = snakemake@output[["pdf"]],
  width = 12.35,
  height = 14.3
)

for (cluster_name in names(plots)) {
  grid::grid.newpage()

  if (
    inherits(plots[[cluster_name]], "grob") ||
      inherits(plots[[cluster_name]], "gTree")
  ) {
    grid::grid.draw(plots[[cluster_name]])
  } else {
    print(plots[[cluster_name]])
  }
}

dev.off()

log_message("Saving PNG")

png(
  filename = snakemake@output[["png"]],
  width = 4550,
  height = 5850,
  res = 300
)

cowplot::plot_grid(
  plotlist = plots,
  ncol = 1
)

dev.off()

log_message("Analysis completed successfully")