source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

count_file <- snakemake@input[["counts"]]
minReads <- snakemake@params[["minReads"]]
mapping_samplesheet <- snakemake@input[["samplesheet"]]


###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

#' Generate the read count matrix
#' @param csv_file Path to the read count matrix.s
#' @return Return a dataframe.
get_peak_matrix <- function(csv_file, peaks) {
  # load file and skip the first line
  # FeatureCounts add an extra line at th ebegining of the read count file so we skip it to directly obtain a dataframe
  raw_counts <- read.csv(file = csv_file, header = TRUE, sep = "\t", skip = 1)
  # Define peak names as "chr:start-end" instead of "interval1234" so it will be easier to import as GRanges objects later.
  peak_names <- raw_counts$Geneid
  # Apply new peak names as rownames
  rownames(raw_counts) <- peak_names
  # The samplenames are the cleaned and sorted bam file names. We remove the unwanted file extention to only keep the samplename given for the mapping.
  raw_samplenames <- basename(grep(".target.markdup.sorted.bam", colnames(raw_counts), value=TRUE))
  samplenames <- gsub(".target.markdup.sorted.bam", "", raw_samplenames)
  samplenames <- sub("^.+[.]", "", samplenames)
  # Remove the extra columns of the file (Geneid  Chr Start End Strand  Length)
  raw_counts <- raw_counts[, colnames(raw_counts) %in% raw_samplenames]
  # Rename colums with the samplenames
  colnames(raw_counts) <- samplenames
  # Transform values as integers for DESeq2 that does not support floats. Normally, ATAC-seq featureCounts are not floats but this is just in case...
  raw_counts <- round(raw_counts, digits = 0)
  # Select peaks on canonical chromosomes
  raw_counts <- raw_counts[grep("^chr[0-9|X|Y]{1,2}[:]", rownames(raw_counts)), ]
  return(raw_counts)
}

#' Normalize the read count using the size factor normalization from DESeq2
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_normalized_counts <- function(raw_counts, samplesheet) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = samplesheet,
    design = ~conditions
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  # norm_counts <- DESeq2::counts(dds, normalized=TRUE)
  norm_counts <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))

  return(norm_counts)
}

#' Get the size factor from Deseq2. It is used to normanize the bigwig files.
#' @param raw_counts Read count matrix.
#' @param samplesheet Samplesheet for DESeq2.
#' @return Return a dataframe.
get_size_factors <- function(raw_counts, samplesheet) {
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

#' When the maximum value (read count) of a peak between samples is under a certain threshold, we considere it is not relevant and the values are set to 0.
#' @param data Read count matrix.
#' @param minExp Minimum value. Default is 5.
#' @return Return a dataframe.
run_filter_low_counts <- function(data, minExp = 5) {
  col_names <- colnames(data)
  data <- t(apply(data, 1, filter_low_counts, col_names = col_names, minExp = minExp))
  data <- as.data.frame(data)
  data <- data[rowSums(data[]) > 0, ]
  return(data)
}

filter_low_counts <- function(row, col_names, minExp) {
  if (max(row) < minExp) {
    return(setNames(rep(0, length(row)), col_names))
  } else {
    return(setNames(row, col_names))
  }
}





###########################################
#                                         #
#              Get matrices               #
#                                         #
###########################################

raw_counts <- get_peak_matrix(
  count_file
)

# Remove peaks if max value < x reads
raw_counts <- run_filter_low_counts(raw_counts, minReads)

###########################################
#                                         #
#             Get samplesheet             #
#                                         #
###########################################

conditions <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 1)

replicate <- sapply(strsplit(colnames(raw_counts), "_"), `[`, 3)

samplesheet <- data.frame(
  sample = colnames(raw_counts),
  conditions = conditions,
  replicate = replicate
)

###########################################
#                                         #
#           Normalize read counts         #
#                                         #
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
#                                         #
#             Prepare bed file            #
#                                         #
###########################################

OCR_GR <- GenomicRanges::GRanges(rownames(norm_counts))


###########################################
#                                         #
#             Order matrices              #
#                                         #
###########################################

if (length(mapping_samplesheet)) {
  samplesheet_table <- read.csv(file=mapping_samplesheet, header = TRUE)
  sample_order <- as.vector(paste0(samplesheet_table$group, "_R", samplesheet_table$replicate))
  raw_counts <- raw_counts[ , na.omit(match(sample_order, colnames(raw_counts)))]
  norm_counts <- norm_counts[ , na.omit(match(sample_order, colnames(norm_counts)))]
  samplesheet <- samplesheet[na.omit(match(sample_order, samplesheet$sample)), ]
}

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

write.csv(raw_counts, snakemake@output[["counts"]], quote=FALSE)
write.csv(norm_counts, snakemake@output[["norm_counts"]], quote=FALSE)
write.csv(samplesheet, snakemake@output[["samplesheet"]], quote=FALSE, row.names=FALSE)
write.csv(size_factors, snakemake@output[["size_factors"]], quote=FALSE)
rtracklayer::export.bed(OCR_GR,con=snakemake@output[["bed"]])
