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

log_message("Starting DAR/TRIM28 motif enrichment analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("GenomeInfoDb")
  library("rtracklayer")
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
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 12))
genome_version <- "mm10"
background <- "genome"

###########################################
# Load input data
###########################################

log_message("Loading DAR peak table")

peaks <- read.csv(
  file = snakemake@input[["DAR_peaks"]],
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!"cluster" %in% colnames(peaks)) {
  log_error("Column 'cluster' not found in DAR peak table")
}

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  snakemake@input[["TRIM28_peaks"]]
)

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
  res_table,
  jaspar_database,
  genome_version,
  background,
  n_cores
) {
  log_message("Preparing PWM motif set")

  pwms <- TFBSTools::getMatrixSet(
    jaspar_database,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )

  if (!"cond" %in% colnames(S4Vectors::mcols(DARs))) {
    log_error("Column 'cond' not found in DAR GRanges metadata")
  }

  seqlevels_keep <- paste0("chr", c(1:19, "X"))

  keep_regions <- as.character(
    GenomicRanges::seqnames(DARs)
  ) %in% seqlevels_keep

  DARs <- DARs[keep_regions]

  if (length(DARs) == 0) {
    log_error("No DARs retained after filtering canonical chromosomes")
  }

  bins <- factor(DARs$cond)

  genome <- get_genome_object(
    genome_version = genome_version
  )

  log_message("Extracting genomic sequences")

  sequences <- Biostrings::getSeq(
    genome,
    DARs
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
      file = res_table,
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
  ) >= 10

  enrichment_selected <- enrichment[selected_motifs, ]

  if (nrow(enrichment_selected) == 0) {
    log_message("No significantly enriched TFBS motifs found")

    write.table(
      "No significantly enriched TFBS motif found",
      file = res_table,
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

  write.table(
    TF_summary,
    file = res_table,
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

merge_TF_motifs <- function(
  enrichment_selected,
  background
) {
  log_message("Preparing motif enrichment matrix")

  TF_enrichment <- SummarizedExperiment::assay(
    enrichment_selected,
    "log2enr"
  )

  TF_pval <- 10^(
    -SummarizedExperiment::assay(
      enrichment_selected,
      "negLog10Padj"
    )
  )

  colnames(TF_pval) <- paste0(
    "p-val_",
    colnames(TF_pval)
  )

  enrichment_for_clustering <- cbind(
    TF_enrichment,
    TF_pval
  )

  n_clusters <- min(
    7,
    nrow(enrichment_for_clustering)
  )

  hcl <- stats::hclust(
    stats::dist(enrichment_for_clustering),
    method = "ward.D"
  )

  clustering <- stats::cutree(
    hcl,
    k = n_clusters
  )

  log_message("Building motif signatures")

  motif_signatures <- lapply(seq_len(n_clusters), function(cluster_id) {
    TFs <- names(clustering[clustering == cluster_id])

    if (length(TFs) > 5) {
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
      leaves <- names(phylog$leaves)
      pfms <- pfms[leaves]

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

  logo_grobs <- lapply(
    motifs_pfms,
    monaLisa::seqLogoGrob,
    xmax = max_width,
    xjust = "center"
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
      "p-val|mat_names",
      colnames(enrichment),
      invert = TRUE
    ),
    drop = FALSE
  ]

  matrix[matrix > 1] <- 1
  matrix[matrix < -1] <- -1

  if (background == "genome") {
    log_message("Masking non-significant enrichment values")

    pval <- enrichment[
      ,
      grep("p-val", colnames(enrichment)),
      drop = FALSE
    ]

    matrix[pval > 0.01] <- NA
  }

  matrix_no_na <- matrix
  matrix_no_na[is.na(matrix_no_na)] <- -10

  matrix <- matrix[
    do.call(
      order,
      unname(as.data.frame(matrix_no_na))
    ),
    ,
    drop = FALSE
  ]

  logo_grobs <- logo_grobs[
    rownames(matrix)
  ]

  heatmap_annotation <- ComplexHeatmap::HeatmapAnnotation(
    logo = monaLisa::annoSeqlogo(
      grobL = logo_grobs,
      which = "row",
      space = grid::unit(0.5, "mm"),
      width = grid::unit(1.5, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  condition_colours <- rep(
    "black",
    ncol(matrix)
  )

  names(condition_colours) <- colnames(matrix)

  condition_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Conditions = ComplexHeatmap::anno_block(
      gp = grid::gpar(
        fill = condition_colours,
        col = 0
      ),
      labels = names(condition_colours),
      labels_gp = grid::gpar(
        col = "white",
        fontsize = 14
      ),
      height = grid::unit(7, "mm")
    )
  )

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

  cell_height <- ifelse(
    nrow(matrix) < 10,
    10,
    6
  )

  log_message("Drawing motif enrichment heatmap")

  heatmap <- ComplexHeatmap::Heatmap(
    matrix,
    name = "Log2 enrichment",
    right_annotation = heatmap_annotation,
    top_annotation = condition_annotation,
    column_split = colnames(matrix),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_title = NULL,
    row_title = NULL,
    col = colour_function,
    width = ncol(matrix) * grid::unit(15, "mm"),
    height = nrow(matrix) * grid::unit(cell_height, "mm"),
    row_names_max_width = grid::unit(25, "cm"),
    heatmap_legend_param = list(
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

###########################################
# Prepare DAR groups
###########################################

log_message("Preparing DAR groups")

DAR_peaks <- GenomicRanges::GRanges(
  rownames(peaks)
)

DAR_peaks$cluster <- peaks$cluster

log_message("Finding overlaps with TRIM28")

hits <- GenomicRanges::findOverlaps(
  DAR_peaks,
  TRIM28_peaks,
  ignore.strand = TRUE
)

query_hits <- unique(
  S4Vectors::queryHits(hits)
)

DAR_only <- DAR_peaks[
  setdiff(seq_along(DAR_peaks), query_hits)
]

common <- DAR_peaks[
  query_hits
]

DAR_only$cond <- paste0(
  DAR_only$cluster,
  "_TRIM28-"
)

common$cond <- paste0(
  common$cluster,
  "_TRIM28+"
)

DARs <- c(
  common,
  DAR_only
)

log_message("DARs overlapping TRIM28: ", length(common))
log_message("DARs not overlapping TRIM28: ", length(DAR_only))

###########################################
# Run motif enrichment
###########################################

log_message("Running motif enrichment")

enrichments <- get_enriched_TFs(
  DARs = DARs,
  res_table = res_table,
  jaspar_database = JASPAR,
  genome_version = genome_version,
  background = background,
  n_cores = n_cores
)

###########################################
# Plot results
###########################################

log_message("Plotting results")

if (!is.null(enrichments) &&
    nrow(enrichments) > 5) {
  figure <- merge_TF_motifs(
    enrichment_selected = enrichments,
    background = background
  )
} else {
  figure <- ggplot2::ggplot() +
    ggplot2::theme_void()
}

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 28,
  base_height = 25,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 28,
  base_height = 25,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")