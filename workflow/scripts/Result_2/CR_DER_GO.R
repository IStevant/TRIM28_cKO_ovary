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

log_message("Starting GO enrichment analysis for genes overlapping differential H3K9me3 regions")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("cowplot")
  library("ggplot2")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

AB <- snakemake@params[["AB"]]
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

log_message("Loading differential H3K9me3 regions")

DER <- read.table(
  snakemake@input[["DER"]],
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

if (nrow(DER) == 0) {
  log_error("DER table is empty")
}

log_message("Number of tested H3K9me3 regions: ", nrow(DER))

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

###########################################
# Functions
###########################################

prepare_gene_annotation <- function(
  genome_gtf
) {
  log_message("Preparing gene annotation")

  genes <- genome_gtf[
    genome_gtf$type == "gene"
  ]

  if (length(genes) == 0) {
    log_error("No gene features found in genome annotation")
  }

  gene_metadata <- as.data.frame(
    GenomicRanges::mcols(genes)
  )

  if (!all(c("gene_id", "gene_name") %in% colnames(gene_metadata))) {
    log_error("Genome annotation must contain 'gene_id' and 'gene_name' columns")
  }

  genes$gene_id <- gene_metadata$gene_id
  genes$gene_name <- gene_metadata$gene_name

  genes
}

prepare_differential_regions <- function(
  DER,
  adj_pval,
  log2FC
) {
  log_message("Filtering significant differential H3K9me3 regions")

  required_columns <- c(
    "padj",
    "log2FoldChange"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(DER)
  )

  if (length(missing_columns) > 0) {
    log_error(
      "Missing required columns in DER table: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  DER <- DER[
    !is.na(DER$padj) &
      DER$padj < adj_pval &
      abs(DER$log2FoldChange) >= log2FC,
    ,
    drop = FALSE
  ]

  if (nrow(DER) == 0) {
    log_error("No significant differential H3K9me3 regions after filtering")
  }

  DER$cluster <- ifelse(
    DER$log2FoldChange >= log2FC,
    "Higher H3K9me3 in WT",
    "Higher H3K9me3 in KO"
  )

  log_message("Significant differential H3K9me3 regions: ", nrow(DER))

  DER
}

get_overlapping_genes <- function(
  DER,
  genes
) {
  log_message("Finding genes overlapping differential H3K9me3 regions")

  DER_gr <- GenomicRanges::GRanges(
    rownames(DER)
  )

  DER_gr$cluster <- DER$cluster

  hits <- GenomicRanges::findOverlaps(
    DER_gr,
    genes,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    log_error("No genes overlap the significant differential H3K9me3 regions")
  }

  overlapping_genes <- data.frame(
    region = rownames(DER)[S4Vectors::queryHits(hits)],
    gene_id = genes$gene_id[S4Vectors::subjectHits(hits)],
    gene = genes$gene_name[S4Vectors::subjectHits(hits)],
    cluster = DER_gr$cluster[S4Vectors::queryHits(hits)],
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

  log_message("Number of region-gene overlaps: ", nrow(overlapping_genes))
  log_message("Number of unique overlapping genes: ", length(unique(overlapping_genes$gene)))

  overlapping_genes
}

run_GO_enrichment <- function(
  overlapping_genes,
  output_file
) {
  log_message("Running GO term enrichment")

  gene_clusters <- unique(
    overlapping_genes[, c("gene", "cluster")]
  )

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
    log_message("No enriched GO terms found")

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
        "GO enrichment of genes overlapping differential H3K9me3 regions"
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

genes <- prepare_gene_annotation(
  genome_gtf = genome_gtf
)

DER_filtered <- prepare_differential_regions(
  DER = DER,
  adj_pval = adj_pval,
  log2FC = log2FC
)

overlapping_genes <- get_overlapping_genes(
  DER = DER_filtered,
  genes = genes
)

gene_file <- file.path(
  output_folder,
  paste0("CR_", AB, "_genes_overlapping_H3K9me3_DER.tsv")
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
  paste0("CR_", AB, "_GO_genes_overlapping_H3K9me3_DER.tsv")
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