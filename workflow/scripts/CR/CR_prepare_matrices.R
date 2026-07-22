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


###########################################
# Snakemake inputs, outputs and parameters
###########################################

log_message("Reading Snakemake inputs, outputs and parameters.")

# Inputs
count_file <- snakemake@input[["counts"]]
mapping_samplesheet <- snakemake@input[["samplesheet"]]

# Parameters
min_reads <- snakemake@params[["minReads"]]

# Outputs
output_counts <- snakemake@output[["counts"]]
output_norm_counts <- snakemake@output[["norm_counts"]]
output_samplesheet <- snakemake@output[["samplesheet"]]
output_size_factors <- snakemake@output[["size_factors"]]
output_bed <- snakemake@output[["bed"]]


###########################################
# Validate inputs
###########################################

if (!file.exists(count_file)) {
  log_error(
    "FeatureCounts file does not exist: ",
    count_file
  )
}

if (length(min_reads) != 1 || is.na(min_reads) || min_reads < 0) {
  log_error(
    "The minimum read-count threshold must be a single non-negative value."
  )
}


###########################################
# Functions
###########################################

#' Generate the read-count matrix
#'
#' @param csv_file Path to the FeatureCounts read-count file.
#'
#' @return A numeric read-count data frame.
get_peak_matrix <- function(csv_file) {

  log_message("Loading FeatureCounts matrix: ", csv_file)

  # FeatureCounts adds an extra line at the beginning of the output file.
  raw_counts <- read.csv(
    file = csv_file,
    header = TRUE,
    sep = "\t",
    skip = 1,
    check.names = FALSE
  )

  required_columns <- c(
    "Geneid",
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
  )

  if (!all(required_columns %in% colnames(raw_counts))) {
    log_error(
      "The FeatureCounts file does not contain all expected annotation columns."
    )
  }

  # Use genomic coordinates as peak identifiers.
  peak_names <- raw_counts$Geneid
  rownames(raw_counts) <- peak_names

  # Identify BAM-derived sample columns.
  raw_sample_names <- basename(
    grep(
      "\\.target\\.markdup\\.sorted\\.bam$",
      colnames(raw_counts),
      value = TRUE
    )
  )

  if (length(raw_sample_names) == 0) {
    log_error(
      "No columns ending with '.target.markdup.sorted.bam' were found."
    )
  }

  sample_names <- gsub(
    "\\.target\\.markdup\\.sorted\\.bam$",
    "",
    raw_sample_names
  )

  sample_names <- sub(
    "^.+\\.",
    "",
    sample_names
  )

  raw_counts <- raw_counts[
    ,
    raw_sample_names,
    drop = FALSE
  ]

  colnames(raw_counts) <- sample_names

  # DESeq2 requires integer counts.
  raw_counts <- round(
    raw_counts,
    digits = 0
  )

  # Retain peaks located on canonical chromosomes.
  canonical_peaks <- grep(
    "^chr([0-9]{1,2}|X|Y):",
    rownames(raw_counts)
  )

  raw_counts <- raw_counts[
    canonical_peaks,
    ,
    drop = FALSE
  ]

  log_message(
    "Read-count matrix generated with ",
    nrow(raw_counts),
    " peaks and ",
    ncol(raw_counts),
    " samples."
  )

  return(raw_counts)
}


#' Normalise read counts using DESeq2
#'
#' @param raw_counts Raw read-count matrix.
#' @param samplesheet Sample information for DESeq2.
#'
#' @return A variance-stabilised read-count matrix.
get_normalized_counts <- function(
    raw_counts,
    samplesheet) {

  log_message("Computing variance-stabilised read counts.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(
    dds
  )

  norm_counts <- SummarizedExperiment::assay(
    DESeq2::vst(
      dds,
      blind = FALSE
    )
  )

  return(norm_counts)
}


#' Extract DESeq2 size factors
#'
#' The size factors are subsequently used to normalise BigWig files.
#'
#' @param raw_counts Raw read-count matrix.
#' @param samplesheet Sample information for DESeq2.
#'
#' @return A data frame containing one size factor per sample.
get_size_factors <- function(
    raw_counts,
    samplesheet) {

  log_message("Computing DESeq2 size factors.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(
    dds
  )

  size_factors <- DESeq2::sizeFactors(
    dds
  )

  size_factors <- data.frame(
    sample = names(size_factors),
    SizeFactor = size_factors
  )

  return(size_factors)
}


#' Filter a single row of read counts
#'
#' @param row Vector of read counts.
#' @param column_names Sample names.
#' @param min_expression Minimum read-count threshold.
#'
#' @return The original row or a zero-filled row.
filter_low_counts <- function(
    row,
    column_names,
    min_expression) {

  if (max(row) < min_expression) {
    return(
      setNames(
        rep(0, length(row)),
        column_names
      )
    )
  }

  return(
    setNames(
      row,
      column_names
    )
  )
}


#' Remove peaks with low read counts
#'
#' Peaks whose maximum read count across all samples is below the
#' selected threshold are removed.
#'
#' @param data Read-count matrix.
#' @param min_expression Minimum read-count threshold.
#'
#' @return A filtered read-count data frame.
run_filter_low_counts <- function(
    data,
    min_expression = 5) {

  log_message(
    "Filtering peaks with a maximum read count below ",
    min_expression,
    "."
  )

  initial_peak_count <- nrow(data)
  column_names <- colnames(data)

  data <- t(
    apply(
      data,
      1,
      filter_low_counts,
      column_names = column_names,
      min_expression = min_expression
    )
  )

  data <- as.data.frame(
    data
  )

  data <- data[
    rowSums(data) > 0,
    ,
    drop = FALSE
  ]

  log_message(
    "Retained ",
    nrow(data),
    " of ",
    initial_peak_count,
    " peaks."
  )

  return(data)
}


###########################################
# Generate read-count matrix
###########################################

raw_counts <- get_peak_matrix(
  count_file
)

raw_counts <- run_filter_low_counts(
  raw_counts,
  min_expression = min_reads
)

if (nrow(raw_counts) == 0) {
  log_error(
    "No peaks remained after low-count filtering."
  )
}


###########################################
# Generate DESeq2 sample sheet
###########################################

log_message("Generating the DESeq2 sample sheet.")

sample_name_parts <- strsplit(
  colnames(raw_counts),
  "_"
)

conditions <- vapply(
  sample_name_parts,
  `[`,
  character(1),
  1
)

replicates <- vapply(
  sample_name_parts,
  `[`,
  character(1),
  3
)

samplesheet <- data.frame(
  sample = colnames(raw_counts),
  conditions = conditions,
  replicate = replicates,
  row.names = colnames(raw_counts),
  stringsAsFactors = FALSE
)


###########################################
# Normalise read counts
###########################################

norm_counts <- get_normalized_counts(
  raw_counts,
  samplesheet
)

size_factors <- get_size_factors(
  raw_counts,
  samplesheet
)


###########################################
# Prepare BED regions
###########################################

log_message("Converting peak identifiers to genomic regions.")

peak_granges <- GenomicRanges::GRanges(
  rownames(norm_counts)
)


###########################################
# Order samples
###########################################

has_mapping_samplesheet <- (
  length(mapping_samplesheet) > 0 &&
    !is.na(mapping_samplesheet[1]) &&
    nzchar(mapping_samplesheet[1])
)

if (has_mapping_samplesheet) {

  if (!file.exists(mapping_samplesheet)) {
    log_error(
      "Mapping sample sheet does not exist: ",
      mapping_samplesheet
    )
  }

  log_message(
    "Ordering samples according to: ",
    mapping_samplesheet
  )

  samplesheet_table <- read.csv(
    file = mapping_samplesheet,
    header = TRUE
  )

  required_columns <- c(
    "group",
    "replicate"
  )

  if (!all(required_columns %in% colnames(samplesheet_table))) {
    log_error(
      "The mapping sample sheet must contain 'group' and 'replicate' columns."
    )
  }

  sample_order <- paste0(
    samplesheet_table$group,
    "_R",
    samplesheet_table$replicate
  )

  matched_samples <- na.omit(
    match(
      sample_order,
      colnames(raw_counts)
    )
  )

  raw_counts <- raw_counts[
    ,
    matched_samples,
    drop = FALSE
  ]

  norm_counts <- norm_counts[
    ,
    matched_samples,
    drop = FALSE
  ]

  samplesheet <- samplesheet[
    match(
      colnames(raw_counts),
      samplesheet$sample
    ),
    ,
    drop = FALSE
  ]
}


###########################################
# Export results
###########################################

log_message("Writing raw read-count matrix: ", output_counts)

write.csv(
  raw_counts,
  file = output_counts,
  quote = FALSE
)

log_message(
  "Writing normalised read-count matrix: ",
  output_norm_counts
)

write.csv(
  norm_counts,
  file = output_norm_counts,
  quote = FALSE
)

log_message(
  "Writing DESeq2 sample sheet: ",
  output_samplesheet
)

write.csv(
  samplesheet,
  file = output_samplesheet,
  quote = FALSE,
  row.names = FALSE
)

log_message(
  "Writing DESeq2 size factors: ",
  output_size_factors
)

write.csv(
  size_factors,
  file = output_size_factors,
  quote = FALSE,
  row.names = FALSE
)

log_message(
  "Writing peak BED file: ",
  output_bed
)

rtracklayer::export.bed(
  peak_granges,
  con = output_bed
)

log_message("Read-count processing completed.")