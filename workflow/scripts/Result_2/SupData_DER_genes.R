source(".Rprofile")

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}

log_message("Starting H3K9me3 DER gene / DEG summary script")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
})

adj_pval <- as.numeric(snakemake@params[["adjpval"]])
log2FC <- as.numeric(snakemake@params[["log2FC"]])

DER <- read.table(
  snakemake@input[["DER"]],
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

bulk_DEG <- read.csv(
  snakemake@input[["bulk_DEG"]],
  header = TRUE,
  check.names = FALSE
)

scRNA_DEG <- read.csv(
  snakemake@input[["scRNA_DEG"]],
  header = TRUE,
  check.names = FALSE
)

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

if (is.na(adj_pval)) {
  log_error("Parameter 'adjpval' must be numeric")
}

if (is.na(log2FC)) {
  log_error("Parameter 'log2FC' must be numeric")
}

if (!all(c("padj", "log2FoldChange") %in% colnames(DER))) {
  log_error("DER table must contain 'padj' and 'log2FoldChange' columns")
}

prepare_gene_annotation <- function(genome_gtf) {
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

prepare_DER <- function(
  DER,
  adj_pval,
  log2FC
) {
  log_message("Filtering significant H3K9me3 DERs")

  DER <- DER[
    !is.na(DER$padj) &
      DER$padj < adj_pval &
      abs(DER$log2FoldChange) >= log2FC,
    ,
    drop = FALSE
  ]

  if (nrow(DER) == 0) {
    log_error("No significant H3K9me3 DERs after filtering")
  }

  DER$DER_cluster <- ifelse(
    DER$log2FoldChange >= log2FC,
    "a",
    "b"
  )

  log_message("Significant H3K9me3 DERs retained: ", nrow(DER))

  DER
}

get_overlapping_genes <- function(
  DER,
  genes
) {
  log_message("Finding genes overlapping H3K9me3 DERs")

  DER_regions <- GenomicRanges::GRanges(
    rownames(DER)
  )

  DER_regions$DER_cluster <- DER$DER_cluster

  hits <- GenomicRanges::findOverlaps(
    DER_regions,
    genes,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    log_error("No genes overlap H3K9me3 DERs")
  }

  overlapping_genes <- data.frame(
    gene_id = genes$gene_id[S4Vectors::subjectHits(hits)],
    gene = genes$gene_name[S4Vectors::subjectHits(hits)],
    DER_cluster = DER_regions$DER_cluster[S4Vectors::queryHits(hits)],
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

  overlapping_genes <- overlapping_genes[
    order(
      overlapping_genes$DER_cluster,
      overlapping_genes$gene
    ),
    ,
    drop = FALSE
  ]

  log_message("Unique gene / cluster associations: ", nrow(overlapping_genes))
  log_message("Unique genes: ", length(unique(overlapping_genes$gene)))

  overlapping_genes
}

prepare_DEG_status <- function(
  DEG,
  table_name
) {
  log_message("Preparing DEG status for ", table_name)

  if (!all(c("Gene", "DEG") %in% colnames(DEG))) {
    log_error(
      table_name,
      " DEG table must contain 'Gene' and 'DEG' columns"
    )
  }

  DEG_status <- data.frame(
    gene = as.character(DEG$Gene),
    status = as.character(DEG$DEG),
    stringsAsFactors = FALSE
  )

  DEG_status$status <- toupper(DEG_status$status)

  DEG_status$status[
    !DEG_status$status %in% c("UP", "DOWN")
  ] <- "-"

  colnames(DEG_status) <- c(
    "gene",
    paste0(table_name, "_status")
  )

  DEG_status <- DEG_status[
    !duplicated(DEG_status$gene),
    ,
    drop = FALSE
  ]

  DEG_status
}

add_DEG_status <- function(
  overlapping_genes,
  DEG_status,
  table_name
) {
  log_message("Adding ", table_name, " DEG status")

  output <- merge(
    overlapping_genes,
    DEG_status,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  status_column <- paste0(table_name, "_status")

  output[[status_column]][is.na(output[[status_column]])] <- "-"

  output
}

genes <- prepare_gene_annotation(
  genome_gtf = genome_gtf
)

DER_filtered <- prepare_DER(
  DER = DER,
  adj_pval = adj_pval,
  log2FC = log2FC
)

overlapping_genes <- get_overlapping_genes(
  DER = DER_filtered,
  genes = genes
)

bulk_status <- prepare_DEG_status(
  DEG = bulk_DEG,
  table_name = "bulk"
)

scRNA_status <- prepare_DEG_status(
  DEG = scRNA_DEG,
  table_name = "scRNAseq"
)

summary_table <- add_DEG_status(
  overlapping_genes = overlapping_genes,
  DEG_status = bulk_status,
  table_name = "bulk"
)

summary_table <- add_DEG_status(
  overlapping_genes = summary_table,
  DEG_status = scRNA_status,
  table_name = "scRNAseq"
)

summary_table <- summary_table[, c(
  "gene_id",
  "gene",
  "DER_cluster",
  "bulk_status",
  "scRNAseq_status"
)]

summary_table <- unique(
  summary_table
)

summary_table <- summary_table[
  order(
    summary_table$DER_cluster,
    summary_table$gene
  ),
  ,
  drop = FALSE
]

log_message("Final table rows: ", nrow(summary_table))
log_message("Writing output table")

write.table(
  summary_table,
  file = snakemake@output[["table"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

log_message("Analysis completed successfully")