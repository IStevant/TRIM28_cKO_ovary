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
  library("randomcoloR")
  library("stringr")
  library("ade4")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 24))
background <- snakemake@params[["background"]]
logos <- snakemake@params[["logos"]]
genome_version <- "mm10"

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

log_message("Loading significant differential regions")

filtered_DARs <- read.csv(
  file = snakemake@input[["sig_DARs"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (nrow(filtered_DARs) == 0) {
  log_error("No significant differential regions found")
}

if (!"cluster" %in% colnames(filtered_DARs)) {
  log_error("Column 'cluster' containing the cluster/bin assignment was not found")
}

log_message("Number of differential regions: ", nrow(filtered_DARs))

log_message("Loading TF gene list")

TFs <- as.vector(
  read.csv(
    snakemake@input[["TF_genes"]],
    header = FALSE
  )[, 1]
)

log_message("Loading samplesheet")

samplesheet <- read.csv(
  file = snakemake@input[["samplesheet"]],
  row.names = 1,
  check.names = FALSE
)

conditions <- sort(
  unique(samplesheet$conditions)
)

conditions_color <- randomcoloR::distinctColorPalette(
  length(conditions)
)

names(conditions_color) <- conditions

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

get_enriched_TFs <- function(
  DARs,
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

  peak_gr_all <- GenomicRanges::GRanges(
    rownames(DARs)
  )

  seqlevels_keep <- paste0("chr", c(1:19, "X"))

  keep_regions <- as.character(
    GenomicRanges::seqnames(peak_gr_all)
  ) %in% seqlevels_keep

  peak_gr <- peak_gr_all[keep_regions]
  bins <- factor(DARs$cluster[keep_regions])

  if (length(peak_gr) == 0) {
    log_error("No regions retained after filtering canonical chromosomes")
  }

  log_message("Number of regions used for motif enrichment: ", length(peak_gr))

  genome <- get_genome_object(
    genome_version
  )

  log_message("Extracting genomic sequences")

  sequences <- Biostrings::getSeq(
    genome,
    peak_gr
  )

  log_message("Running motif enrichment")

  if (background == "genome") {
    enrichment <- monaLisa::calcBinnedMotifEnrR(
      seqs = sequences,
      bins = bins,
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
      bins = bins,
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

    return(NULL)
  }

  log_message("Filtering enriched motifs")

  selected_motifs <- apply(
    SummarizedExperiment::assay(enrichment, "negLog10Padj"),
    1,
    function(x) max(abs(x), 0, na.rm = TRUE)
  ) >= TFBS_log10FDR

  enrichment_selected <- enrichment[selected_motifs, ]

  if (nrow(enrichment_selected) == 0) {
    log_message("No significantly enriched motifs after filtering")

    write.table(
      "No significantly enriched TFBS motif found",
      file = result_table,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    return(enrichment_selected)
  }

  log2_enrichment <- SummarizedExperiment::assay(
    enrichment_selected,
    "log2enr"
  )

  FDR <- 10^(
    -SummarizedExperiment::assay(
      enrichment_selected,
      "negLog10Padj"
    )
  )

  colnames(FDR) <- paste0(
    "FDR_",
    colnames(FDR)
  )

  TF_summary <- data.frame(
    TF.name = SummarizedExperiment::rowData(enrichment_selected)$motif.name,
    TF.matrix = rownames(log2_enrichment),
    log2_enrichment,
    FDR,
    check.names = FALSE
  )

  log_message("Writing motif enrichment table")

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

prepare_TF_motif_matrix <- function(
  enrichment_selected,
  background,
  TFBS_FDR
) {
  log_message("Preparing enrichment matrix")

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

  log_message("Clustering motif profiles")

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

      log_message("Building motif similarity tree")

      hc <- motifStack::clusterMotifs(pfms)
      phylog <- ade4::hclust2phylog(hc)
      leaves <- names(phylog$leaves)
      pfms <- pfms[leaves]

      log_message("Building motif signatures")

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

  log_message("Merging similar motifs")

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
          tf_names <- collapse_tf_names(
            tf_vector
          )

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

  enrichment_names <- rownames(enrichment)

  enrichment_names <- gsub("XXX", ".", enrichment_names)
  enrichment_names <- gsub("YYY", "::", enrichment_names)
  enrichment_names <- gsub("ZZZ", "-", enrichment_names)

  rownames(enrichment) <- enrichment_names

  matrix <- enrichment[
    ,
    grep(
      "FDR_|mat_names",
      colnames(enrichment),
      invert = TRUE
    ),
    drop = FALSE
  ]

  matrix[matrix > 1] <- 1
  matrix[matrix < -1] <- -1

  if (background == "genome") {
    log_message("Masking non-significant enrichment values")

    FDR <- enrichment[
      ,
      grep("FDR_", colnames(enrichment)),
      drop = FALSE
    ]

    matrix[FDR > TFBS_FDR] <- NA
  }

  matrix_no_na <- matrix
  matrix_no_na[is.na(matrix_no_na)] <- -10

  matrix <- matrix[
    order(
      -matrix_no_na[, 1],
      matrix_no_na[, 2]
    ),
    ,
    drop = FALSE
  ]

  logo_grobs <- NULL

  if (isTRUE(logos)) {
    log_message("Preparing motif logos")

    logo_grobs <- lapply(
      motifs_pfms,
      monaLisa::seqLogoGrob,
      xmax = max_width,
      xjust = "right"
    )

    names(logo_grobs) <- rownames(enrichment)

    logo_grobs <- logo_grobs[
      rownames(matrix)
    ]
  }

  list(
    matrix = matrix,
    logo_grobs = logo_grobs
  )
}

make_combined_heatmap <- function(
  matrix,
  logo_grobs,
  logos
) {
  log_message("Drawing combined motif enrichment heatmap")

  cold <- grDevices::colorRampPalette(
    c("#04bbc6", "#52d0cf", "#8addd8", "#c1f0e0", "#fffee8")
  )

  warm <- grDevices::colorRampPalette(
    c("#fffee8", "#ffd9cb", "#ffb1ad", "#f5808d", "#eb2d62")
  )

  colour_function <- circlize::colorRamp2(
    breaks = seq(-1, 1, length.out = 24),
    colors = c(cold(12), warm(12))
  )

  condition_colours <- rep(
    "black",
    ncol(matrix)
  )

  names(condition_colours) <- colnames(matrix)

  condition_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Stages = ComplexHeatmap::anno_block(
      gp = grid::gpar(
        fill = condition_colours,
        col = 0
      ),
      labels = names(condition_colours),
      labels_gp = grid::gpar(
        col = "white",
        fontsize = 9.8
      ),
      height = grid::unit(7, "mm")
    )
  )

  logo_annotation <- NULL

  if (isTRUE(logos) && !is.null(logo_grobs)) {
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
  }

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "Log2 enrichment",
    left_annotation = logo_annotation,
    top_annotation = condition_annotation,
    column_split = colnames(matrix),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_title = "TFBS motif enrichment",
    column_title_gp = grid::gpar(
      fontsize = 9.8,
      fontface = "bold"
    ),
    row_title = NULL,
    col = colour_function,
    width = ncol(matrix) * grid::unit(2, "mm"),
    height = nrow(matrix) * grid::unit(2, "mm"),
    rect_gp = grid::gpar(
      col = "white",
      lwd = 0.4
    ),
    row_names_max_width = grid::unit(17.5, "cm"),
    row_names_gp = grid::gpar(
      fontsize = 7
    ),
    heatmap_legend_param = list(
      direction = "horizontal",
      title_gp = grid::gpar(fontsize = 7),
      labels_gp = grid::gpar(fontsize = 7)
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

###########################################
# Run analysis
###########################################

log_message("Running motif enrichment")

enrichments <- get_enriched_TFs(
  DARs = filtered_DARs,
  result_table = res_table,
  jaspar_database = JASPAR,
  genome_version = genome_version,
  background = background,
  n_cores = n_cores,
  TFBS_log10FDR = TFBS_log10FDR
)

if (!is.null(enrichments) &&
    nrow(enrichments) > 5) {
  heatmap_data <- prepare_TF_motif_matrix(
    enrichment_selected = enrichments,
    background = background,
    TFBS_FDR = TFBS_FDR
  )

  figure <- make_combined_heatmap(
    matrix = heatmap_data$matrix,
    logo_grobs = heatmap_data$logo_grobs,
    logos = logos
  )
} else {
  figure <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No significantly enriched TFBS motifs found",
      size = 5
    ) +
    ggplot2::theme_void()
}

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = figure,
  base_width = 52,
  base_height = 39,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = figure,
  base_width = 52,
  base_height = 39,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")