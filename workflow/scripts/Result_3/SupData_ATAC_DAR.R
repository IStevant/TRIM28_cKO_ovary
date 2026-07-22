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

log_message("Starting differential accessibility analysis script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("DESeq2")
  library("GenomicRanges")
  library("rtracklayer")
  library("txdbmaker")
  library("ChIPseeker")
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
# Load input data
###########################################

log_message("Loading count matrix")

raw_counts <- read.csv(
  file = snakemake@input[["counts"]],
  row.names = 1,
  check.names = FALSE
)

raw_counts <- as.matrix(raw_counts)
storage.mode(raw_counts) <- "integer"

log_message("Loading samplesheet")

samplesheet <- read.csv(
  file = snakemake@input[["samplesheet"]],
  row.names = 1,
  check.names = FALSE
)

if (!"conditions" %in% colnames(samplesheet)) {
  log_error("Column 'conditions' not found in samplesheet")
}

samplesheet$conditions <- factor(samplesheet$conditions)

common_samples <- intersect(
  colnames(raw_counts),
  rownames(samplesheet)
)

if (length(common_samples) == 0) {
  log_error("No matching samples between count matrix columns and samplesheet rows")
}

raw_counts <- raw_counts[, common_samples, drop = FALSE]
samplesheet <- samplesheet[common_samples, , drop = FALSE]

conditions <- levels(droplevels(samplesheet$conditions))

if (length(conditions) != 2) {
  log_error(
    "This script is configured for exactly 2 conditions. Current number of conditions: ",
    length(conditions)
  )
}

log_message("Conditions: ", paste(conditions, collapse = ", "))

###########################################
# Load genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  snakemake@input[["genome"]]
)

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(
  snakemake@input[["genome"]]
)

###########################################
# Functions
###########################################

annotate_DARs <- function(
  DAR_GR,
  TxDb,
  gene2symbol,
  promoter
) {
  if (length(DAR_GR) == 0) {
    return(data.frame(
      annotation = character(),
      nearest.gene = character(),
      distanceToTSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      DAR_GR,
      genomicAnnotationPriority = c(
        "Promoter",
        "5UTR",
        "Exon",
        "Intron",
        "3UTR",
        "Downstream",
        "Intergenic"
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

  data.frame(
    annotation = annotation$annotation,
    nearest.gene = annotation$geneId,
    distanceToTSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

###########################################
# Run DESeq2 analysis
###########################################

log_message("Creating DESeq2 object")

DESeqObj <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = samplesheet,
  design = ~conditions
)

log_message("Running DESeq2 Wald test for two conditions")

DARs <- DESeq2::DESeq(
  DESeqObj,
  test = "Wald"
)

contrast <- c(
  "conditions",
  conditions[2],
  conditions[1]
)

log_message(
  "Extracting results for contrast: ",
  conditions[2],
  " vs ",
  conditions[1]
)

resDARs <- DESeq2::results(
  DARs,
  contrast = contrast
)

resDARs <- as.data.frame(resDARs)

resDARs <- resDARs[
  !is.na(resDARs$padj),
  ,
  drop = FALSE
]

sig_DA <- resDARs[
  resDARs$padj < adj_pval &
    abs(resDARs$log2FoldChange) > log2FC,
  ,
  drop = FALSE
]

log_message("Number of significant DARs: ", nrow(sig_DA))

DAR_GR <- GenomicRanges::GRanges(
  rownames(sig_DA)
)

###########################################
# Export significant DAR BED
###########################################

log_message("Writing significant DAR BED file")

rtracklayer::export.bed(
  DAR_GR,
  con = snakemake@output[["sig_DARs_bed"]]
)

###########################################
# Annotate DARs
###########################################

log_message("Annotating significant DARs")

DA_anno <- annotate_DARs(
  DAR_GR = DAR_GR,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)

filtered_DARs <- data.frame(
  sig_DA,
  annotation = DA_anno$annotation,
  nearest.gene = DA_anno$nearest.gene,
  distanceToTSS = DA_anno$distanceToTSS,
  stringsAsFactors = FALSE
)

###########################################
# Save output files
###########################################

log_message("Writing DAR result table")

write.table(
  filtered_DARs[order(filtered_DARs$padj), ],
  file = snakemake@output[["tsv"]],
  quote = FALSE,
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)

log_message("Saving significant DAR IDs")

filtered_DARs <- rownames(filtered_DARs)

save(
  filtered_DARs,
  file = snakemake@output[["sig_DARs"]]
)

log_message("Analysis completed successfully")