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

log_message("Starting bigWig normalisation script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("rtracklayer")
})

###########################################
# Load parameters
###########################################

log_message("Loading parameters")

size_factors_file <- snakemake@input[["size_factors"]]
AB <- snakemake@params[["AB"]]
bigwig_folder <- snakemake@params[["bigwig_folder"]]
new_bigwig_folder <- snakemake@params[["new_bigwig_folder"]]
output_file <- snakemake@output[["output_file"]]

###########################################
# Prepare output directory
###########################################

log_message("Preparing output directory")

dir.create(
  new_bigwig_folder,
  showWarnings = FALSE,
  recursive = TRUE
)

###########################################
# Load size factors
###########################################

log_message("Loading size factors")

size_factors <- read.csv(
  size_factors_file,
  row.names = 1,
  check.names = FALSE
)

if (!"SizeFactor" %in% colnames(size_factors)) {
  log_error("Column 'SizeFactor' not found in size factors file")
}

###########################################
# Find bigWig files
###########################################

log_message("Finding bigWig files")

bigwig_files <- list.files(
  path = bigwig_folder,
  pattern = "\\.bigWig$",
  full.names = FALSE
)

bigwig_files <- grep(
  AB,
  bigwig_files,
  value = TRUE
)

if (length(bigwig_files) == 0) {
  log_error("No bigWig files found for pattern: ", AB)
}

samples <- sub("\\.bigWig$", "", bigwig_files)

log_message("Number of bigWig files selected: ", length(bigwig_files))
log_message("Selected samples: ", paste(samples, collapse = ", "))

###########################################
# Functions
###########################################

get_size_factor <- function(sample, size_factors) {
  log_message("Looking for size factor for sample: ", sample)
  log_message("Available size factor samples: ", paste(rownames(size_factors), collapse = ", "))

  clean_sample <- sub(".*\\.", "", sample)

  matched_rows <- rownames(size_factors)[
    tolower(rownames(size_factors)) == tolower(clean_sample)
  ]

  if (length(matched_rows) == 0) {
    log_error("No size factor found for sample: ", sample, " after cleaning to: ", clean_sample)
  }

  size_factor <- as.numeric(size_factors[matched_rows, "SizeFactor"])

  if (is.na(size_factor)) {
    log_error("Size factor is NA for sample: ", matched_rows)
  }

  size_factor
}

normalise_bigwig <- function(
  sample,
  bigwig_folder,
  new_bigwig_folder,
  size_factors
) {
  log_message("Processing sample: ", sample)

  input_bigwig <- file.path(
    bigwig_folder,
    paste0(sample, ".bigWig")
  )

  output_bigwig <- file.path(
    new_bigwig_folder,
    paste0(sample, ".bw")
  )

  log_message("Importing bigWig: ", input_bigwig)

  bigwig_data <- rtracklayer::import.bw(input_bigwig)

  size_factor <- get_size_factor(
    sample = sample,
    size_factors = size_factors
  )

  log_message("Size factor for ", sample, ": ", size_factor)

  normalised_bigwig <- bigwig_data
  normalised_bigwig$score <- bigwig_data$score * (1 / size_factor)

  log_message("Writing normalised bigWig: ", output_bigwig)

  rtracklayer::export.bw(
    normalised_bigwig,
    output_bigwig
  )

  output_bigwig
}

###########################################
# Normalise bigWig files
###########################################

normalised_bigwig_files <- lapply(
  samples,
  normalise_bigwig,
  bigwig_folder = bigwig_folder,
  new_bigwig_folder = new_bigwig_folder,
  size_factors = size_factors
)

###########################################
# Save output file
###########################################

log_message("Writing output summary file")

output_summary <- data.frame(
  sample = samples,
  normalised_bigwig = unlist(normalised_bigwig_files),
  stringsAsFactors = FALSE
)

write.csv(
  output_summary,
  file = output_file,
  quote = FALSE,
  row.names = FALSE
)

log_message("Analysis completed successfully")