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

log_message("Starting TFBS motif enrichment analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("S4Vectors")
  library("rtracklayer")
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
  library("seriation")
  library("grid")
  library("cowplot")
  library("ggplot2")
  library("stringr")
  library("ade4")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

n_cores <- 12
genome_version <- "mm10"
background <- "genome"

BiocParallel::register(
  BiocParallel::MulticoreParam(n_cores)
)

TFBS_log10FDR <- as.numeric(
  snakemake@params[["TFBS_log10FDR"]]
)

if (is.na(TFBS_log10FDR)) {
  log_error("Parameter 'TFBS_log10FDR' must be numeric")
}

TFBS_FDR <- 10^(-TFBS_log10FDR)

log_message("Using FDR threshold: ", TFBS_FDR)

###########################################
# Load peak data
###########################################

log_message("Importing peak files")

H3K9me3_peaks <- rtracklayer::import(
  snakemake@input[["H3K9me3_peaks"]]
)

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

log_message("Number of H3K9me3 peaks: ", length(H3K9me3_peaks))
log_message("Number of TRIM28 peaks: ", length(TRIM28_peaks))

###########################################
# Load output paths
###########################################

log_message("Loading output paths")

res_table_H3K9me3_only <- snakemake@output[["res_table_H3K9me3_only"]]
res_table_TRIM28_only <- snakemake@output[["res_table_TRIM28_only"]]
res_table_TRIM28_H3K9me3 <- snakemake@output[["res_table_TRIM28_H3K9me3"]]

###########################################
# Load JASPAR database
###########################################

log_message("Loading JASPAR 2024 database")

JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024()@db

###########################################
# Functions
###########################################

make_empty_panel <- function(title) {
  cowplot::plot_grid(
    grid::textGrob(
      title,
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      )
    ),
    grid::textGrob(
      "Pas d'enrichissement",
      gp = grid::gpar(fontsize = 8.5)
    ),
    ncol = 1,
    rel_heights = c(0.12, 1)
  )
}

make_titled_panel <- function(plot, title) {
  cowplot::plot_grid(
    grid::textGrob(
      title,
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      )
    ),
    plot,
    ncol = 1,
    rel_heights = c(0.08, 1)
  )
}

get_genome_object <- function(genome_version) {
  if (genome_version == "mm10") {
    return(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
  }

  log_error("Unsupported genome version: ", genome_version)
}

get_enriched_TFs <- function(
  regions,
  result_table,
  jaspar_database,
  genome_version,
  background,
  n_cores,
  TFBS_log10FDR
) {
  log_message("Preparing PWM motif set")

  pwms <- TFBSTools::getMatrixSet(
    jaspar_database,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  seqlevels_keep <- paste0("chr", c(1:19, "X"))

  regions <- GenomeInfoDb::keepSeqlevels(
    regions,
    seqlevels_keep,
    pruning.mode = "coarse"
  )

  log_message("Number of regions retained after seqlevel filtering: ", length(regions))

  if (length(regions) == 0) {
    log_message("No regions available for motif enrichment")

    write.table(
      "No regions available for motif enrichment",
      file = result_table,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    return(NULL)
  }

  genome <- get_genome_object(genome_version)

  log_message("Extracting genomic sequences")

  sequences <- Biostrings::getSeq(
    genome,
    regions
  )

  log_message("Running motif enrichment analysis")

  if (background == "genome") {
    enrichment <- monaLisa::calcBinnedMotifEnrR(
      seqs = sequences,
      pwmL = pwms,
      background = "genome",
      genome = genome,
      genome.regions = NULL,
      genome.oversample = 2,
      BPPARAM = BiocParallel::MulticoreParam(n_cores)
    )
  } else {
    enrichment <- monaLisa::calcBinnedMotifEnrR(
      seqs = sequences,
      pwmL = pwms,
      BPPARAM = BiocParallel::MulticoreParam(n_cores)
    )
  }

  if (length(as.vector(SummarizedExperiment::rowData(enrichment)$motif.name)) <= 1) {
    log_message("No TFBS enrichment found")

    write.table(
      "No TFBS enrichment found",
      file = result_table,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    return(enrichment)
  }

  log_message("Filtering significantly enriched motifs")

  selected_motifs <- apply(
    SummarizedExperiment::assay(enrichment, "negLog10Padj"),
    1,
    function(x) max(abs(x), 0, na.rm = TRUE)
  ) >= TFBS_log10FDR

  enrichment_selected <- enrichment[selected_motifs, ]

  if (nrow(enrichment_selected) == 0) {
    log_message("No significantly enriched TFBS motifs found")

    write.table(
      "No significantly enriched TFBS motif found",
      file = result_table,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    return(enrichment_selected)
  }

  log_message("Writing motif enrichment result table")

  TF_summary <- data.frame(
    TF.name = SummarizedExperiment::rowData(enrichment_selected)$motif.name,
    TF.matrix = rownames(SummarizedExperiment::assay(enrichment_selected, "log2enr")),
    SummarizedExperiment::assay(enrichment_selected, "log2enr"),
    10^(-SummarizedExperiment::assay(enrichment_selected, "negLog10Padj")),
    check.names = FALSE
  )

  fdr_columns <- seq(
    3 + ncol(SummarizedExperiment::assay(enrichment_selected, "log2enr")),
    ncol(TF_summary)
  )

  colnames(TF_summary)[fdr_columns] <- paste0(
    "FDR_",
    colnames(SummarizedExperiment::assay(enrichment_selected, "negLog10Padj"))
  )

  write.table(
    TF_summary,
    file = result_table,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

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
      return(paste0(common_prefix, "s"))
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

    all_suffixes <- c(
      non_numeric_suffixes,
      numeric_suffixes
    )

    paste0(
      common_prefix,
      paste(all_suffixes, collapse = "/")
    )
  })

  paste(
    new_tf_vector[order(new_tf_vector)],
    collapse = ";"
  )
}

merge_TF_motifs <- function(
  enrichment_selected,
  background,
  TFBS_FDR
) {
  if (is.null(enrichment_selected) || nrow(enrichment_selected) == 0) {
    return(list(
      plot = ggplot2::ggplot() + ggplot2::theme_void(),
      n_rows = 5,
      n_cols = 1
    ))
  }

  log_message("Clustering motifs by enrichment profile")

  TF_enrichment <- SummarizedExperiment::assay(
    enrichment_selected,
    "log2enr"
  )

  TF_fdr <- 10^(
    -SummarizedExperiment::assay(
      enrichment_selected,
      "negLog10Padj"
    )
  )

  colnames(TF_fdr) <- paste0(
    "FDR_",
    colnames(TF_fdr)
  )

  TF_enrichment_for_clustering <- cbind(
    TF_enrichment,
    TF_fdr
  )

  n_clusters <- 1

  hcl <- stats::hclust(
    stats::dist(TF_enrichment_for_clustering),
    method = "ward.D"
  )

  clustering <- stats::cutree(
    hcl,
    k = n_clusters
  )

  motif_signatures <- lapply(seq_len(n_clusters), function(cluster_id) {
    TFs <- names(clustering[clustering == cluster_id])

    if (length(TFs) > 5) {
      log_message("Clustering motif PWMs")

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

      log_message("Building motif similarity tree")

      hc <- motifStack::clusterMotifs(pfms)
      phylog <- ade4::hclust2phylog(hc)
      leaves <- names(phylog$leaves)
      pfms <- pfms[leaves]

      log_message("Building motif signatures")

      motif_signature <- motifStack::motifSignature(
        pfms,
        phylog,
        cutoffPval = 0.01,
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

  log_message("Merging motif signatures")

  motif_signature_enrichment <- lapply(motif_signatures, function(cluster) {
    lapply(cluster, function(motif) {
      TF_names <- motif@name

      lapply(TF_names, function(TF) {
        tf_vector <- toupper(
          unlist(strsplit(TF, ";"))
        )

        enrichment <- TF_enrichment_for_clustering

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

  log_message("Preparing motif logo annotation")

  grob_list <- lapply(
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
  names(grob_list) <- enrichment_names

  matrix <- enrichment[
    ,
    grep(
      "FDR_|mat_names",
      colnames(enrichment),
      invert = TRUE
    ),
    drop = FALSE
  ]

  matrix <- matrix[
    ,
    1,
    drop = FALSE
  ]

  matrix[matrix > 0.5] <- 0.5
  matrix[matrix < -0.5] <- -0.5

  heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(
    logo = monaLisa::annoSeqlogo(
      grobL = grob_list,
      which = "row",
      space = grid::unit(0.5, "mm"),
      width = grid::unit(1.05, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  warm_palette <- grDevices::colorRampPalette(
    c("#fffee8", "#ffd9cb", "#ffb1ad", "#f5808d", "#eb2d62")
  )

  colour_function <- circlize::colorRamp2(
    breaks = seq(0, 0.5, length.out = 12),
    colors = warm_palette(12)
  )

  log_message("Ordering heatmap rows from most to least enriched")

  enrichment_score <- matrix[, 1]
  enrichment_score[is.na(enrichment_score)] <- -Inf

  row_order <- order(
    enrichment_score,
    decreasing = TRUE
  )

  if (background == "genome") {
    log_message("Replacing non-significant enrichment values by NA")

    fdr <- enrichment[
      ,
      grep("FDR_", colnames(enrichment)),
      drop = FALSE
    ]

    fdr <- fdr[
      ,
      1,
      drop = FALSE
    ]

    matrix[fdr > TFBS_FDR] <- NA
  }

  log_message("Drawing motif enrichment heatmap")

  cell_size <- grid::unit(4, "mm")

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "Log2 enrichment",
    left_annotation = heatmap_annotation,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = row_order,
    column_title = NULL,
    row_title = NULL,
    col = colour_function,
    width = ncol(matrix) * cell_size,
    height = nrow(matrix) * cell_size,
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.5
    ),
    row_names_max_width = grid::unit(17.5, "cm"),
    row_names_gp = grid::gpar(fontsize = 7),
    heatmap_legend_param = list(
      direction = "horizontal",
      labels_gp = grid::gpar(fontsize = 7),
      title_gp = grid::gpar(fontsize = 7)
    )
  )

  heatmap_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      heatmap,
      heatmap_legend_side = "bottom"
    ),
    wrap.grobs = TRUE
  )

  list(
    plot = heatmap_grob,
    n_rows = nrow(matrix),
    n_cols = ncol(matrix)
  )
}

###########################################
# Identify peak groups
###########################################

log_message("Finding overlaps between H3K9me3 and TRIM28 peaks")

hits <- GenomicRanges::findOverlaps(
  H3K9me3_peaks,
  TRIM28_peaks
)

query_hits <- S4Vectors::queryHits(hits)
subject_hits <- S4Vectors::subjectHits(hits)

H3K9me3_only <- H3K9me3_peaks[-query_hits]
TRIM28_only <- TRIM28_peaks[-subject_hits]

both <- GenomicRanges::pintersect(
  H3K9me3_peaks[query_hits],
  TRIM28_peaks[subject_hits]
)

H3K9me3_ov_TRIM28 <- H3K9me3_peaks[unique(query_hits)]
TRIM28_ov_H3K9me3 <- TRIM28_peaks[unique(subject_hits)]

log_message("Number of H3K9me3-only peaks: ", length(H3K9me3_only))
log_message("Number of TRIM28-only peaks: ", length(TRIM28_only))
log_message("Number of intersected peak regions: ", length(both))
log_message("Number of H3K9me3 peaks overlapping TRIM28: ", length(H3K9me3_ov_TRIM28))
log_message("Number of TRIM28 peaks overlapping H3K9me3: ", length(TRIM28_ov_H3K9me3))

###########################################
# Run motif enrichment
###########################################

log_message("Running motif enrichment for H3K9me3-only peaks")

enrichments_H3K9me3_only <- get_enriched_TFs(
  regions = H3K9me3_only,
  result_table = res_table_H3K9me3_only,
  jaspar_database = JASPAR,
  genome_version = genome_version,
  background = background,
  n_cores = n_cores,
  TFBS_log10FDR = TFBS_log10FDR
)

log_message("Running motif enrichment for TRIM28-only peaks")

enrichments_TRIM28_only <- get_enriched_TFs(
  regions = TRIM28_only,
  result_table = res_table_TRIM28_only,
  jaspar_database = JASPAR,
  genome_version = genome_version,
  background = background,
  n_cores = n_cores,
  TFBS_log10FDR = TFBS_log10FDR
)

log_message("Running motif enrichment for TRIM28 peaks overlapping H3K9me3")

enrichments_TRIM28_ov_H3K9me3 <- get_enriched_TFs(
  regions = TRIM28_ov_H3K9me3,
  result_table = res_table_TRIM28_H3K9me3,
  jaspar_database = JASPAR,
  genome_version = genome_version,
  background = background,
  n_cores = n_cores,
  TFBS_log10FDR = TFBS_log10FDR
)

###########################################
# Plot results
###########################################

log_message("Plotting motif enrichment results")

if (!is.null(enrichments_H3K9me3_only) &&
    nrow(enrichments_H3K9me3_only) > 1) {
  result1 <- merge_TF_motifs(
    enrichment_selected = enrichments_H3K9me3_only,
    background = background,
    TFBS_FDR = TFBS_FDR
  )

  figure1 <- make_titled_panel(
    result1$plot,
    "H3K9me3-only peaks"
  )

  n_rows_figure1 <- result1$n_rows
} else {
  figure1 <- make_empty_panel("H3K9me3-only peaks")
  n_rows_figure1 <- 5
}

if (!is.null(enrichments_TRIM28_only) &&
    nrow(enrichments_TRIM28_only) > 1) {
  result2 <- merge_TF_motifs(
    enrichment_selected = enrichments_TRIM28_only,
    background = background,
    TFBS_FDR = TFBS_FDR
  )

  figure2 <- make_titled_panel(
    result2$plot,
    "TRIM28-only peaks"
  )

  n_rows_figure2 <- result2$n_rows
} else {
  figure2 <- make_empty_panel("TRIM28-only peaks")
  n_rows_figure2 <- 5
}

if (!is.null(enrichments_TRIM28_ov_H3K9me3) &&
    nrow(enrichments_TRIM28_ov_H3K9me3) > 1) {
  result3 <- merge_TF_motifs(
    enrichment_selected = enrichments_TRIM28_ov_H3K9me3,
    background = background,
    TFBS_FDR = TFBS_FDR
  )

  figure3 <- make_titled_panel(
    result3$plot,
    "TRIM28 peaks overlapping H3K9me3"
  )

  n_rows_figure3 <- result3$n_rows
} else {
  figure3 <- make_empty_panel("TRIM28 peaks overlapping H3K9me3")
  n_rows_figure3 <- 5
}

###########################################
# Save output files
###########################################

log_message("Saving figures")

max_n_rows <- max(
  n_rows_figure1,
  n_rows_figure2,
  n_rows_figure3
)

pdf_width <- 13
pdf_height <- max(
  5,
  (max_n_rows * 4 / 25.4) + 2.5
)

pdf(
  file = snakemake@output[["pdf"]],
  width = pdf_width,
  height = pdf_height
)

print(figure1)
print(figure2)
print(figure3)

dev.off()

png(
  filename = snakemake@output[["png"]],
  width = pdf_width,
  height = pdf_height * 3,
  units = "in",
  res = 300
)

print(
  cowplot::plot_grid(
    figure1,
    figure2,
    figure3,
    ncol = 1,
    rel_heights = c(1, 1, 1)
  )
)

dev.off()

log_message("Analysis completed successfully")