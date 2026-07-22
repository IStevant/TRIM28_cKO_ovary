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

log_message("Starting differential peak analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("DESeq2")
  library("GenomicRanges")
  library("rtracklayer")
  library("ChIPseeker")
  library("txdbmaker")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

promoter <- as.numeric(snakemake@params[["promoter"]])
adj_pval <- as.numeric(snakemake@params[["adjpval"]])
log2FC <- as.numeric(snakemake@params[["log2FC"]])

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}

if (is.na(adj_pval)) {
  log_error("Parameter 'adjpval' must be numeric")
}

if (is.na(log2FC)) {
  log_error("Parameter 'log2FC' must be numeric")
}

###########################################
# Load data
###########################################

log_message("Loading count matrix")

raw_counts <- read.csv(
  file = snakemake@input[["counts"]],
  row.names = 1,
  check.names = FALSE
)

raw_counts <- as.matrix(raw_counts)
storage.mode(raw_counts) <- "numeric"
raw_counts <- round(raw_counts)

if (nrow(raw_counts) == 0) {
  log_error("Count matrix has zero rows")
}

if (ncol(raw_counts) == 0) {
  log_error("Count matrix has zero columns")
}

if (anyNA(raw_counts)) {
  log_error("Count matrix contains NA values")
}

log_message("Number of peaks: ", nrow(raw_counts))
log_message("Number of samples: ", ncol(raw_counts))

log_message("Loading samplesheet")

samplesheet <- read.csv(
  file = snakemake@input[["samplesheet"]],
  row.names = 1,
  check.names = FALSE
)

if (!"conditions" %in% colnames(samplesheet)) {
  log_error("Samplesheet must contain a 'conditions' column")
}

samplesheet <- samplesheet[colnames(raw_counts), , drop = FALSE]

if (any(is.na(rownames(samplesheet)))) {
  log_error("Samplesheet rows do not match count matrix columns")
}

samplesheet$conditions <- factor(samplesheet$conditions)
conditions <- levels(samplesheet$conditions)

log_message("Conditions detected: ", paste(conditions, collapse = ", "))

if (length(conditions) < 2) {
  log_error("Not enough conditions to perform differential analysis")
}

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_file <- snakemake@input[["genome"]]

genome_gtf <- rtracklayer::import(
  genome_file
)

log_message("Preparing gene ID to symbol mapping")

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(
  genome_file
)

###########################################
# Functions
###########################################

annotate_regions <- function(
  regions,
  TxDb,
  gene2symbol,
  promoter
) {
  if (length(regions) == 0) {
    return(data.frame(
      region = character(),
      annotation = character(),
      nearest_gene = character(),
      distance_to_TSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      regions,
      genomicAnnotationPriority = c(
        "Promoter", "5UTR", "Exon", "Intron",
        "3UTR", "Downstream", "Intergenic"
      ),
      tssRegion = c(-promoter, 0),
      TxDb = TxDb,
      level = "gene",
      overlap = "all"
    )
  )

  annotation$geneId <- gene2symbol[
    annotation$geneId,
    "gene_name"
  ]

  annotation$region <- paste0(
    annotation$seqnames,
    ":",
    annotation$start,
    "-",
    annotation$end
  )

  annotation <- annotation[, c(
    "region",
    "annotation",
    "geneId",
    "distanceToTSS"
  )]

  colnames(annotation) <- c(
    "region",
    "annotation",
    "nearest_gene",
    "distance_to_TSS"
  )

  annotation
}

run_deseq_analysis <- function(
  raw_counts,
  samplesheet
) {
  log_message("Creating DESeq2 object")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  conditions <- levels(samplesheet$conditions)

  if (length(conditions) != 2) {
    log_error("This script expects exactly two conditions")
  }

  log_message(
    "Running DESeq2 Wald test: ",
    conditions[2],
    " versus ",
    conditions[1]
  )

  dds <- DESeq2::DESeq(
    dds,
    test = "Wald"
  )

  res <- DESeq2::results(
    dds,
    contrast = c(
      "conditions",
      conditions[2],
      conditions[1]
    )
  )

  res
}

###########################################
# Run differential analysis
###########################################

resDERs <- run_deseq_analysis(
  raw_counts = raw_counts,
  samplesheet = samplesheet
)

unfiltered_DERs <- as.data.frame(resDERs)

log_message("Filtering significant differential regions")

sig_DE <- subset(
  resDERs,
  !is.na(padj) &
    padj < adj_pval &
    abs(log2FoldChange) > log2FC
)

log_message("Number of significant differential regions: ", nrow(sig_DE))

DER_GR <- GenomicRanges::GRanges(
  rownames(sig_DE)
)

###########################################
# Annotate significant regions
###########################################

log_message("Annotating significant differential regions")

DE_anno <- annotate_regions(
  regions = DER_GR,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)

filtered_DERs <- data.frame(
  region = rownames(sig_DE),
  as.data.frame(sig_DE),
  stringsAsFactors = FALSE
)

filtered_DERs <- merge(
  filtered_DERs,
  DE_anno,
  by = "region",
  all.x = TRUE,
  sort = FALSE
)

rownames(filtered_DERs) <- filtered_DERs$region

filtered_DERs <- filtered_DERs[order(filtered_DERs$padj), ]

###########################################
# Save output files
###########################################

log_message("Writing filtered differential region table")

write.table(
  filtered_DERs,
  file = snakemake@output[["tsv"]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

log_message("Writing unfiltered differential region table")

write.table(
  unfiltered_DERs,
  file = paste0(snakemake@output[["tsv"]], "_no_filter"),
  quote = FALSE,
  sep = "\t",
  row.names = TRUE
)

log_message("Saving significant differential region IDs")

filtered_DER_ids <- rownames(filtered_DERs)

save(
  filtered_DER_ids,
  file = snakemake@output[["sig_DERs"]]
)

log_message("Writing BED file")

DER_bed <- GenomicRanges::GRanges(
  filtered_DER_ids
)

if (length(DER_bed) > 0) {
  GenomicRanges::mcols(DER_bed) <- filtered_DERs[
    ,
    setdiff(colnames(filtered_DERs), "region"),
    drop = FALSE
  ]
}

rtracklayer::export(
  DER_bed,
  snakemake@output[["DER_bed"]],
  format = "bed"
)

log_message("Analysis completed successfully")