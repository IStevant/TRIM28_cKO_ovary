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

log_message("Starting count matrix preparation script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("DESeq2")
  library("SummarizedExperiment")
  library("GenomicRanges")
  library("rtracklayer")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

count_file <- snakemake@input[["counts"]]
min_reads <- as.numeric(snakemake@params[["minReads"]])
mapping_samplesheet <- snakemake@input[["samplesheet"]]

if (is.na(min_reads)) {
  log_error("Parameter 'minReads' must be numeric")
}

###########################################
# Functions
###########################################

get_peak_matrix <- function(count_file) {
  log_message("Reading featureCounts table")

  raw_counts <- read.csv(
    file = count_file,
    header = TRUE,
    sep = "\t",
    skip = 1,
    check.names = FALSE
  )

  peak_names <- paste0(
    raw_counts$Chr,
    ":",
    raw_counts$Start,
    "-",
    raw_counts$End
  )

  rownames(raw_counts) <- peak_names

  raw_samplenames <- grep(
    "\\.mLb\\.clN\\.sorted\\.bam$",
    colnames(raw_counts),
    value = TRUE
  )

  if (length(raw_samplenames) == 0) {
    log_error("No sample columns ending with '.mLb.clN.sorted.bam' were found")
  }

  samplenames <- gsub(
    "\\.mLb\\.clN\\.sorted\\.bam$",
    "",
    basename(raw_samplenames)
  )

  raw_counts <- raw_counts[, raw_samplenames, drop = FALSE]
  colnames(raw_counts) <- samplenames

  raw_counts <- as.data.frame(
    lapply(raw_counts, as.numeric),
    row.names = rownames(raw_counts),
    check.names = FALSE
  )

  raw_counts <- round(raw_counts)

  raw_counts <- raw_counts[
    grep("^chr([0-9]+|X|Y):", rownames(raw_counts)),
    ,
    drop = FALSE
  ]

  if (nrow(raw_counts) == 0) {
    log_error("No peaks retained on canonical chromosomes")
  }

  if (ncol(raw_counts) == 0) {
    log_error("Count matrix has zero sample columns")
  }

  log_message("Number of peaks imported: ", nrow(raw_counts))
  log_message("Number of samples imported: ", ncol(raw_counts))
  log_message("Count matrix class after import: ", paste(class(raw_counts), collapse = ", "))
  log_message("Count matrix column classes: ", paste(unique(vapply(raw_counts, class, character(1))), collapse = ", "))

  raw_counts
}

run_filter_low_counts <- function(counts, min_reads) {
  log_message("Filtering peaks with max read count < ", min_reads)

  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  keep <- apply(
    counts,
    1,
    function(row) max(row, na.rm = TRUE) >= min_reads
  )

  filtered_counts <- counts[keep, , drop = FALSE]

  if (nrow(filtered_counts) == 0) {
    log_error("No peaks retained after low-count filtering")
  }

  log_message("Number of peaks retained after filtering: ", nrow(filtered_counts))
  log_message("Count matrix mode after filtering: ", mode(filtered_counts))

  filtered_counts
}

make_samplesheet <- function(counts) {
  log_message("Preparing samplesheet")

  sample_names <- colnames(counts)

  conditions <- sapply(
    strsplit(sample_names, "_"),
    `[`,
    1
  )

  replicate <- sapply(
    strsplit(sample_names, "_"),
    `[`,
    3
  )

  samplesheet <- data.frame(
    sample = sample_names,
    conditions = conditions,
    replicate = replicate,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )

  log_message("Number of samples in samplesheet: ", nrow(samplesheet))

  samplesheet
}

order_by_samplesheet <- function(
  raw_counts,
  norm_counts = NULL,
  samplesheet,
  mapping_samplesheet
) {
  if (
    is.null(mapping_samplesheet) ||
    length(mapping_samplesheet) == 0 ||
    is.na(mapping_samplesheet) ||
    mapping_samplesheet == "" ||
    !file.exists(mapping_samplesheet)
  ) {
    log_message("No valid mapping samplesheet found. Keeping current sample order")

    return(list(
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      samplesheet = samplesheet
    ))
  }

  log_message("Ordering matrices using mapping samplesheet")

  samplesheet_table <- read.csv(
    file = mapping_samplesheet,
    header = TRUE,
    stringsAsFactors = FALSE
  )

print(samplesheet_table)

  sample_order <- samplesheet_table$sample

  log_message("Requested sample order: ", paste(sample_order, collapse = ", "))
  log_message("Available count columns: ", paste(colnames(raw_counts), collapse = ", "))

  matched_samples <- match(sample_order, colnames(raw_counts))
  matched_samples <- matched_samples[!is.na(matched_samples)]

  if (length(matched_samples) == 0) {
    log_error("No samples from mapping samplesheet matched the count matrix columns")
  }

  raw_counts <- raw_counts[, matched_samples, drop = FALSE]

  if (ncol(raw_counts) == 0) {
    log_error("Reordered count matrix has zero columns")
  }

  samplesheet <- samplesheet[
    match(colnames(raw_counts), samplesheet$sample),
    ,
    drop = FALSE
  ]

  if (anyNA(samplesheet$sample)) {
    log_error("Samplesheet order does not match count matrix columns")
  }

  if (!is.null(norm_counts)) {
    norm_counts <- norm_counts[, colnames(raw_counts), drop = FALSE]
  }

  log_message("Final sample order: ", paste(colnames(raw_counts), collapse = ", "))

  list(
    raw_counts = raw_counts,
    norm_counts = norm_counts,
    samplesheet = samplesheet
  )
}

get_normalized_counts <- function(raw_counts, samplesheet) {
  log_message("Normalising counts with DESeq2 VST")

  raw_counts <- as.matrix(raw_counts)
  storage.mode(raw_counts) <- "numeric"
  raw_counts <- round(raw_counts)

  log_message("Count matrix mode before DESeq2: ", mode(raw_counts))
  log_message("Count matrix dimensions before DESeq2: ", paste(dim(raw_counts), collapse = " x "))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  SummarizedExperiment::assay(
    DESeq2::vst(dds, blind = FALSE)
  )
}

get_size_factors <- function(raw_counts, samplesheet) {
  log_message("Computing DESeq2 size factors")

  raw_counts <- as.matrix(raw_counts)
  storage.mode(raw_counts) <- "numeric"
  raw_counts <- round(raw_counts)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  size_factors <- DESeq2::sizeFactors(dds)

  data.frame(
    sample = names(size_factors),
    SizeFactor = size_factors,
    row.names = NULL
  )
}

###########################################
# Generate count matrices
###########################################

raw_counts <- get_peak_matrix(count_file)

raw_counts <- run_filter_low_counts(
  counts = raw_counts,
  min_reads = min_reads
)

###########################################
# Prepare samplesheet
###########################################

samplesheet <- make_samplesheet(raw_counts)

ordered_objects <- order_by_samplesheet(
  raw_counts = raw_counts,
  samplesheet = samplesheet,
  mapping_samplesheet = mapping_samplesheet
)

raw_counts <- ordered_objects$raw_counts
samplesheet <- ordered_objects$samplesheet

###########################################
# Normalise read counts
###########################################

norm_counts <- get_normalized_counts(
  raw_counts = raw_counts,
  samplesheet = samplesheet
)

size_factors <- get_size_factors(
  raw_counts = raw_counts,
  samplesheet = samplesheet
)

###########################################
# Prepare BED file
###########################################

log_message("Preparing BED file")

OCR_GR <- GenomicRanges::GRanges(rownames(norm_counts))

###########################################
# Save output files
###########################################

log_message("Writing output files")

write.csv(
  raw_counts,
  snakemake@output[["counts"]],
  quote = FALSE
)

write.csv(
  norm_counts,
  snakemake@output[["norm_counts"]],
  quote = FALSE
)

write.csv(
  samplesheet,
  snakemake@output[["samplesheet"]],
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  size_factors,
  snakemake@output[["size_factors"]],
  quote = FALSE,
  row.names = FALSE
)

rtracklayer::export.bed(
  OCR_GR,
  con = snakemake@output[["bed"]]
)

log_message("Analysis completed successfully")