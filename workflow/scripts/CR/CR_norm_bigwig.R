source(".Rprofile")

suppressPackageStartupMessages({
  library("rtracklayer")
})


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
size_factor_file <- snakemake@input[["size_factors"]]

# Parameters
AB <- snakemake@params[["AB"]]
bigwig_folder <- snakemake@params[["bigwig_folder"]]
new_bigwig_folder <- snakemake@params[["new_bigwig_folder"]]

# Outputs
output_file <- snakemake@output[["output_file"]]


###########################################
# Validate inputs
###########################################

if (!file.exists(size_factor_file)) {
  log_error(
    "Size-factor file does not exist: ",
    size_factor_file
  )
}

if (!dir.exists(bigwig_folder)) {
  log_error(
    "BigWig folder does not exist: ",
    bigwig_folder
  )
}

if (length(AB) != 1 || is.na(AB) || !nzchar(AB)) {
  log_error("A unique antibody name must be provided.")
}

dir.create(
  new_bigwig_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


###########################################
# Load size factors
###########################################

log_message("Loading DESeq2 size factors.")

size_factors <- read.csv(
  size_factor_file,
  row.names = 1,
  check.names = FALSE
)

if (!"SizeFactor" %in% colnames(size_factors)) {
  log_error(
    "The size-factor table does not contain a 'SizeFactor' column."
  )
}


###########################################
# Identify BigWig files
###########################################

log_message(
  "Searching for BigWig files associated with antibody: ",
  AB
)

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
  log_error(
    "No BigWig files found for antibody ",
    AB,
    " in: ",
    bigwig_folder
  )
}

samples <- sub(
  "\\.bigWig$",
  "",
  bigwig_files
)

log_message(
  "Found ",
  length(samples),
  " BigWig file(s)."
)


###########################################
# Normalise BigWig files
###########################################

for (sample in samples) {

  log_message("Processing sample: ", sample)

  input_bigwig <- file.path(
    bigwig_folder,
    paste0(sample, ".bigWig")
  )

  log_message(
    "Importing BigWig file: ",
    input_bigwig
  )

  bigwig_data <- rtracklayer::import.bw(
    input_bigwig
  )

  size_factor_matches <- grep(
    tolower(sample),
    tolower(rownames(size_factors))
  )

  if (length(size_factor_matches) == 0) {
    log_error(
      "No DESeq2 size factor found for sample: ",
      sample
    )
  }

  if (length(size_factor_matches) > 1) {
    log_error(
      "Multiple DESeq2 size factors found for sample: ",
      sample
    )
  }

  size_factor <- size_factors[
    size_factor_matches,
    "SizeFactor"
  ]

  if (
    length(size_factor) != 1 ||
      is.na(size_factor) ||
      !is.finite(size_factor) ||
      size_factor <= 0
  ) {
    log_error(
      "Invalid DESeq2 size factor for sample ",
      sample,
      ": ",
      size_factor
    )
  }

  log_message(
    "DESeq2 size factor: ",
    size_factor
  )

  normalised_bigwig <- bigwig_data
  normalised_bigwig$score <- bigwig_data$score / size_factor

  output_bigwig <- file.path(
    new_bigwig_folder,
    paste0(sample, ".bw")
  )

  log_message(
    "Exporting normalised BigWig: ",
    output_bigwig
  )

  rtracklayer::export.bw(
    normalised_bigwig,
    output_bigwig
  )
}


###########################################
# Export size factors
###########################################

log_message(
  "Writing size-factor table: ",
  output_file
)

write.csv(
  size_factors,
  file = output_file,
  quote = FALSE
)

log_message("BigWig normalisation completed.")