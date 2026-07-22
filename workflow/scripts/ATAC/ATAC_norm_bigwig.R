source(".Rprofile")

suppressPackageStartupMessages({
  library("doParallel")
  library("foreach")
  library("rtracklayer")
})

doParallel::registerDoParallel(cores = 12)


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

size_factor_file <- snakemake@input[["size_factors"]]

bigwig_folder <- snakemake@params[["bigwig_folder"]]
new_bigwig_folder <- snakemake@params[["new_bigwig_folder"]]

output_file <- snakemake@output[["output_file"]]


###########################################
# Load data
###########################################

if (!file.exists(size_factor_file)) {
  log_error("Size factor file does not exist: ", size_factor_file)
}

log_message("Loading DESeq2 size factors.")

size_factors <- read.csv(
  size_factor_file,
  row.names = 1
)


###########################################
# Identify BigWig files
###########################################

pattern <- "\\.mLb\\.clN\\.bigWig$"

bigwig_files <- list.files(
  path = bigwig_folder,
  pattern = pattern
)

if (length(bigwig_files) == 0) {
  log_error(
    "No BigWig files found in: ",
    bigwig_folder
  )
}

samples <- gsub(
  "\\.mLb\\.clN\\.bigWig$",
  "",
  bigwig_files
)

scale_folder <- file.path(
  bigwig_folder,
  "scale"
)

if (!dir.exists(scale_folder)) {
  log_error(
    "Scale folder does not exist: ",
    scale_folder
  )
}

dir.create(
  new_bigwig_folder,
  recursive = TRUE,
  showWarnings = FALSE
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

  bw_file <- file.path(
    bigwig_folder,
    paste0(sample, ".mLb.clN.bigWig")
  )

  log_message("Importing BigWig file.")

  bw_data <- rtracklayer::import.bw(
    bw_file
  )

  scale_file <- list.files(
    path = scale_folder,
    pattern = paste0("^", sample)
  )

  if (length(scale_file) == 0) {
    log_error(
      "No scaling factor found for sample: ",
      sample
    )
  }

  if (length(scale_file) > 1) {
    log_error(
      "Multiple scaling files found for sample: ",
      sample
    )
  }

  scale <- read.csv(
    file.path(scale_folder, scale_file),
    header = FALSE
  )[1, 1]

  log_message(
    "Original scaling factor: ",
    scale
  )

  denorm_bw <- bw_data
  denorm_bw$score <- bw_data$score / scale

  size_factor <- size_factors[
    grep(
      tolower(sample),
      tolower(rownames(size_factors))
    ),
    "SizeFactor"
  ]

  if (length(size_factor) == 0) {
    log_error(
      "No DESeq2 size factor found for sample: ",
      sample
    )
  }

  if (length(size_factor) > 1) {
    log_error(
      "Multiple DESeq2 size factors found for sample: ",
      sample
    )
  }

  norm_bw <- denorm_bw
  norm_bw$score <- denorm_bw$score / size_factor

  output_bw <- file.path(
    new_bigwig_folder,
    paste0(sample, ".bw")
  )

  log_message(
    "Exporting normalised BigWig: ",
    output_bw
  )

  rtracklayer::export.bw(
    norm_bw,
    output_bw
  )
}


log_message(
  "Writing size factor table: ",
  output_file
)

write.csv(
  size_factors,
  file = output_file
)

log_message("BigWig normalisation completed.")