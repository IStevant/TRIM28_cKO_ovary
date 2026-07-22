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

count_file <- snakemake@input[["counts"]]
mapping_samplesheet <- snakemake@input[["samplesheet"]]

minReads <- snakemake@params[["minReads"]]

output_counts <- snakemake@output[["counts"]]
output_norm_counts <- snakemake@output[["norm_counts"]]
output_samplesheet <- snakemake@output[["samplesheet"]]
output_size_factors <- snakemake@output[["size_factors"]]
output_bed <- snakemake@output[["bed"]]


###########################################
# Functions
###########################################

#' Generate the read count matrix
#'
#' @param csv_file Path to the featureCounts read count table.
#'
#' @return A data frame containing the raw read counts.
get_peak_matrix <- function(csv_file) {

  log_message("Loading featureCounts table: ", csv_file)

  if (!file.exists(csv_file)) {
    log_error("Read count file does not exist: ", csv_file)
  }

  # featureCounts adds an extra line at the beginning of the count file.
  raw_counts <- read.csv(
    file = csv_file,
    header = TRUE,
    sep = "\t",
    skip = 1
  )

  required_columns <- c(
    "Chr",
    "Start",
    "End"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(raw_counts)
  )

  if (length(missing_columns) > 0) {
    log_error(
      "The featureCounts table is missing the following columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  # Define peak names as "chr:start-end" to facilitate conversion to
  # GRanges objects.
  peak_names <- paste0(
    raw_counts$Chr,
    ":",
    raw_counts$Start,
    "-",
    raw_counts$End
  )

  rownames(raw_counts) <- peak_names

  # Extract sample columns from the cleaned and sorted BAM filenames.
  raw_samplenames <- grep(
    "\\.mLb\\.clN\\.sorted\\.bam$",
    colnames(raw_counts),
    value = TRUE
  )

  if (length(raw_samplenames) == 0) {
    log_error(
      "No sample columns ending with '.mLb.clN.sorted.bam' were found in: ",
      csv_file
    )
  }

  samplenames <- gsub(
    "\\.mLb\\.clN\\.sorted\\.bam$",
    "",
    raw_samplenames
  )

  # Remove the featureCounts annotation columns.
  raw_counts <- raw_counts[
    ,
    raw_samplenames,
    drop = FALSE
  ]

  colnames(raw_counts) <- samplenames

  # Ensure integer counts for DESeq2.
  raw_counts <- round(
    raw_counts,
    digits = 0
  )

  # Retain peaks located on canonical chromosomes.
  raw_counts <- raw_counts[
    grep(
      "^chr([0-9]{1,2}|X|Y):",
      rownames(raw_counts)
    ),
    ,
    drop = FALSE
  ]

  if (nrow(raw_counts) == 0) {
    log_error("No peaks remained on canonical chromosomes.")
  }

  log_message(
    "Loaded ",
    scales::comma(nrow(raw_counts)),
    " peaks across ",
    ncol(raw_counts),
    " samples."
  )

  return(raw_counts)
}


#' Normalise read counts using DESeq2 size-factor normalisation
#'
#' @param raw_counts Read count matrix.
#' @param samplesheet Sample information for DESeq2.
#'
#' @return A matrix containing variance-stabilised read counts.
get_normalized_counts <- function(raw_counts, samplesheet) {

  log_message("Normalising read counts with DESeq2 VST.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  norm_counts <- SummarizedExperiment::assay(
    DESeq2::vst(
      dds,
      blind = FALSE
    )
  )

  return(norm_counts)
}


#' Calculate DESeq2 size factors
#'
#' The size factors can be used to normalise bigWig files.
#'
#' @param raw_counts Read count matrix.
#' @param samplesheet Sample information for DESeq2.
#'
#' @return A data frame containing one size factor per sample.
get_size_factors <- function(raw_counts, samplesheet) {

  log_message("Calculating DESeq2 size factors.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )

  dds <- DESeq2::estimateSizeFactors(dds)

  size_factors <- DESeq2::sizeFactors(dds)

  size_factors <- data.frame(
    sample = names(size_factors),
    SizeFactor = size_factors
  )

  return(size_factors)
}


#' Filter peaks with low read counts
#'
#' Peaks whose maximum read count across samples is lower than `minExp`
#' are removed.
#'
#' @param data Read count matrix.
#' @param minExp Minimum read count. Default is 5.
#'
#' @return A data frame containing the retained peaks.
run_filter_low_counts <- function(data, minExp = 5) {

  log_message(
    "Filtering peaks with a maximum read count below ",
    minExp,
    "."
  )

  initial_peak_count <- nrow(data)
  col_names <- colnames(data)

  data <- t(
    apply(
      data,
      1,
      filter_low_counts,
      col_names = col_names,
      minExp = minExp
    )
  )

  data <- as.data.frame(data)

  data <- data[
    rowSums(data) > 0,
    ,
    drop = FALSE
  ]

  log_message(
    "Retained ",
    scales::comma(nrow(data)),
    " of ",
    scales::comma(initial_peak_count),
    " peaks after filtering."
  )

  if (nrow(data) == 0) {
    log_error(
      "No peaks remained after filtering with minReads = ",
      minExp,
      "."
    )
  }

  return(data)
}


# Set all values of a low-count peak to zero.
filter_low_counts <- function(row, col_names, minExp) {

  if (max(row) < minExp) {
    return(
      setNames(
        rep(0, length(row)),
        col_names
      )
    )
  }

  return(
    setNames(
      row,
      col_names
    )
  )
}


###########################################
# Generate raw count matrix
###########################################

log_message("Generating raw peak count matrix.")

raw_counts <- get_peak_matrix(
  count_file
)

raw_counts <- run_filter_low_counts(
  raw_counts,
  minExp = minReads
)


###########################################
# Generate sample information
###########################################

log_message("Generating DESeq2 sample information.")

conditions <- sapply(
  strsplit(
    colnames(raw_counts),
    "_"
  ),
  `[`,
  1
)

replicate <- sapply(
  strsplit(
    colnames(raw_counts),
    "_"
  ),
  `[`,
  2
)

samplesheet <- data.frame(
  sample = colnames(raw_counts),
  conditions = conditions,
  replicate = replicate,
  row.names = colnames(raw_counts)
)

log_message(
  "Detected conditions: ",
  paste(
    unique(samplesheet$conditions),
    collapse = ", "
  )
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

log_message("Preparing consensus peak regions for BED export.")

OCR_GR <- GenomicRanges::GRanges(
  rownames(norm_counts)
)


###########################################
# Order matrices
###########################################

# When a mapping samplesheet is available, order the matrices according to
# the sample order defined in that file.
if (length(mapping_samplesheet) > 0 && nzchar(mapping_samplesheet)) {

  log_message(
    "Ordering matrices according to mapping samplesheet: ",
    mapping_samplesheet
  )

  if (!file.exists(mapping_samplesheet)) {
    log_error(
      "Mapping samplesheet does not exist: ",
      mapping_samplesheet
    )
  }

  samplesheet_table <- read.csv(
    file = mapping_samplesheet,
    header = TRUE
  )

  if (!"sample" %in% colnames(samplesheet_table)) {
    log_error(
      "The mapping samplesheet does not contain a 'sample' column: ",
      mapping_samplesheet
    )
  }

  sample_order <- as.vector(
    samplesheet_table$sample
  )

  missing_samples <- setdiff(
    sample_order,
    colnames(raw_counts)
  )

  if (length(missing_samples) > 0) {
    log_error(
      "The following samples from the mapping samplesheet were not found ",
      "in the count matrix: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  raw_counts <- raw_counts[
    ,
    match(sample_order, colnames(raw_counts)),
    drop = FALSE
  ]

  norm_counts <- norm_counts[
    ,
    match(sample_order, colnames(norm_counts)),
    drop = FALSE
  ]

  samplesheet <- samplesheet[
    match(sample_order, samplesheet$sample),
    ,
    drop = FALSE
  ]

  size_factors <- size_factors[
    match(sample_order, size_factors$sample),
    ,
    drop = FALSE
  ]
}


###########################################
# Export results
###########################################

log_message("Writing raw read count matrix: ", output_counts)

write.csv(
  raw_counts,
  output_counts,
  quote = FALSE
)

log_message("Writing normalised read count matrix: ", output_norm_counts)

write.csv(
  norm_counts,
  output_norm_counts,
  quote = FALSE
)

log_message("Writing sample information: ", output_samplesheet)

write.csv(
  samplesheet,
  output_samplesheet,
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing DESeq2 size factors: ", output_size_factors)

write.csv(
  size_factors,
  output_size_factors,
  quote = FALSE,
  row.names = FALSE
)

log_message("Writing consensus peak BED file: ", output_bed)

rtracklayer::export.bed(
  OCR_GR,
  con = output_bed
)

log_message("Peak count matrix preparation completed.")