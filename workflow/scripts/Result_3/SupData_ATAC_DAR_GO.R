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

log_message("Starting GO enrichment analysis for genes overlapping DARs")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("GenomicFeatures")
  library("S4Vectors")
  library("rtracklayer")
  library("txdbmaker")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("cowplot")
  library("ggplot2")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

adj_pval <- as.numeric(snakemake@params[["adjpval"]])
log2FC <- as.numeric(snakemake@params[["log2FC"]])
output_folder <- snakemake@params[["path"]]

if (is.na(adj_pval)) {
  log_error("Parameter 'adjpval' must be numeric")
}

if (is.na(log2FC)) {
  log_error("Parameter 'log2FC' must be numeric")
}

dir.create(
  output_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

###########################################
# Load input data
###########################################

log_message("Loading DAR table")

DAR <- read.table(
  snakemake@input[["DAR"]],
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

if (nrow(DAR) == 0) {
  log_error("DAR table is empty")
}

required_columns <- c(
  "padj",
  "log2FoldChange"
)

missing_columns <- setdiff(
  required_columns,
  colnames(DAR)
)

if (length(missing_columns) > 0) {
  log_error(
    "Missing required columns in DAR table: ",
    paste(missing_columns, collapse = ", ")
  )
}

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

log_message("Building TxDb object")

txdb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

###########################################
# Functions
###########################################

prepare_gene_ranges <- function(
  txdb,
  genome_gtf
) {
  log_message("Preparing gene ranges")

  gene_ranges <- GenomicFeatures::genes(
    txdb
  )

  gene_annotation <- genome_gtf[
    genome_gtf$type == "gene"
  ]

  gene_metadata <- as.data.frame(
    GenomicRanges::mcols(gene_annotation)
  )

  if (!all(c("gene_id", "gene_name") %in% colnames(gene_metadata))) {
    log_error("Genome annotation must contain 'gene_id' and 'gene_name' columns")
  }

  gene2symbol <- unique(
    gene_metadata[, c("gene_id", "gene_name")]
  )

  rownames(gene2symbol) <- gene2symbol$gene_id

  gene_ranges$gene_id <- names(gene_ranges)
  gene_ranges$gene_name <- gene2symbol[
    gene_ranges$gene_id,
    "gene_name"
  ]

  gene_ranges <- gene_ranges[
    !is.na(gene_ranges$gene_name) &
      gene_ranges$gene_name != ""
  ]

  gene_ranges
}

prepare_significant_DARs <- function(
  DAR,
  adj_pval,
  log2FC
) {
  log_message("Filtering significant DARs")

  DAR_filtered <- DAR[
    !is.na(DAR$padj) &
      DAR$padj < adj_pval &
      abs(DAR$log2FoldChange) >= log2FC,
    ,
    drop = FALSE
  ]

  if (nrow(DAR_filtered) == 0) {
    log_error("No significant DARs found after filtering")
  }

  DAR_filtered$cluster <- ifelse(
    DAR_filtered$log2FoldChange >= log2FC,
    "Up in WT",
    "Up in KO"
  )

  log_message("Significant DARs retained: ", nrow(DAR_filtered))

  DAR_filtered
}

get_overlapping_genes <- function(
  DAR_filtered,
  gene_ranges
) {
  log_message("Finding genes directly overlapping DARs")

  DAR_gr <- GenomicRanges::GRanges(
    rownames(DAR_filtered)
  )

  DAR_gr$cluster <- DAR_filtered$cluster

  hits <- GenomicRanges::findOverlaps(
    DAR_gr,
    gene_ranges,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    log_error("No genes directly overlap the significant DARs")
  }

  overlapping_genes <- data.frame(
    region = rownames(DAR_filtered)[S4Vectors::queryHits(hits)],
    gene_id = gene_ranges$gene_id[S4Vectors::subjectHits(hits)],
    gene = gene_ranges$gene_name[S4Vectors::subjectHits(hits)],
    cluster = DAR_gr$cluster[S4Vectors::queryHits(hits)],
    stringsAsFactors = FALSE
  )

  overlapping_genes <- overlapping_genes[
    !is.na(overlapping_genes$gene) &
      overlapping_genes$gene != "",
    ,
    drop = FALSE
  ]

  overlapping_genes <- unique(
    overlapping_genes
  )

  log_message("Region-gene overlaps: ", nrow(overlapping_genes))
  log_message("Unique overlapping genes: ", length(unique(overlapping_genes$gene)))

  overlapping_genes
}

run_GO_enrichment <- function(
  overlapping_genes,
  output_file
) {
  log_message("Running GO enrichment")

  gene_clusters <- unique(
    overlapping_genes[, c("gene", "cluster")]
  )

  gene_clusters <- gene_clusters[
    !is.na(gene_clusters$gene) &
      gene_clusters$gene != "",
    ,
    drop = FALSE
  ]

  if (nrow(gene_clusters) == 0) {
    log_error("No genes available for GO enrichment")
  }

  GO_terms <- clusterProfiler::compareCluster(
    gene ~ cluster,
    data = gene_clusters,
    fun = "enrichGO",
    keyType = "SYMBOL",
    OrgDb = org.Mm.eg.db::org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  write.table(
    as.data.frame(GO_terms),
    file = output_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )

  GO_terms
}

plot_GO_terms <- function(
  GO_terms,
  nb_terms = 10
) {
  GO_table <- as.data.frame(
    GO_terms
  )

  if (nrow(GO_table) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = "No enriched GO terms found",
          size = 5
        ) +
        ggplot2::theme_void()
    )
  }

  GO_terms@compareClusterResult$Description <- gsub(
    "^([a-z])",
    "\\U\\1",
    GO_terms@compareClusterResult$Description,
    perl = TRUE
  )

  options(
    enrichplot.colours = c(
      "#77BFA3",
      "#98C9A3",
      "#BFD8BD",
      "#DDE7C7",
      "#EDEEC9"
    )
  )

  clusterProfiler::dotplot(
    GO_terms,
    showCategory = nb_terms
  ) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::ggtitle(
      paste0(
        "GO enrichment of genes directly overlapping DARs"
      )
    ) +
    ggplot2::labs(
      colour = "Adj. p-value",
      size = "Gene ratio"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(reverse = TRUE),
      size = ggplot2::guide_legend(reverse = TRUE)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 12,
        hjust = 0.5,
        face = "bold"
      ),
      axis.title.x = ggplot2::element_blank()
    )
}

###########################################
# Run analysis
###########################################

gene_ranges <- prepare_gene_ranges(
  txdb = txdb,
  genome_gtf = genome_gtf
)

DAR_filtered <- prepare_significant_DARs(
  DAR = DAR,
  adj_pval = adj_pval,
  log2FC = log2FC
)

overlapping_genes <- get_overlapping_genes(
  DAR_filtered = DAR_filtered,
  gene_ranges = gene_ranges
)

gene_file <- file.path(
  output_folder,
  "ATAC_genes_directly_overlapping_DAR.tsv"
)

log_message("Writing overlapping gene table")

write.table(
  overlapping_genes,
  file = gene_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

GO_file <- file.path(
  output_folder,
  "ATAC_GO_genes_directly_overlapping_DAR.tsv"
)

GO_terms <- run_GO_enrichment(
  overlapping_genes = overlapping_genes,
  output_file = GO_file
)

GO_plot <- plot_GO_terms(
  GO_terms = GO_terms,
  nb_terms = 10
)

###########################################
# Save output files
###########################################

log_message("Saving figures")

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = GO_plot,
  base_width = 20,
  base_height = 25,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = GO_plot,
  base_width = 20,
  base_height = 25,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")