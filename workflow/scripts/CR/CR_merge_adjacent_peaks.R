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

peak_folder <- snakemake@params[["raw_peaks"]]
distance <- snakemake@params[["distance"]]

output_folder <- snakemake@output[["output_folder"]]
output_file <- snakemake@output[["peak_counts"]]


###########################################
# Validate inputs
###########################################

if (!dir.exists(peak_folder)) {
  log_error("Peak folder does not exist: ", peak_folder)
}

if (length(distance) == 0) {
  log_error("No peak-merging distances were provided.")
}

if (is.null(names(distance)) || any(names(distance) == "")) {
  log_error(
    "Each peak-merging distance must be associated with an antibody name."
  )
}

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


###########################################
# Merge peaks
###########################################

# For each antibody, merge peaks separated by no more than the distance
# defined in the analysis parameters.
for (antibody in names(distance)) {

  merge_distance <- as.numeric(
    distance[[antibody]]
  )

  log_message(
    "Processing antibody ",
    antibody,
    " with a merging distance of ",
    merge_distance,
    " bp."
  )

  bed_files <- list.files(
    path = peak_folder,
    pattern = antibody
  )

  if (length(bed_files) == 0) {
    log_error(
      "No peak files found for antibody: ",
      antibody
    )
  }

  for (bed_file in bed_files) {

    input_bed <- file.path(
      peak_folder,
      bed_file
    )

    output_bed <- file.path(
      output_folder,
      bed_file
    )

    log_message("Merging peaks from: ", bed_file)

    bedtoolsr::bt.merge(
      i = input_bed,
      d = merge_distance,
      output = output_bed
    )
  }
}


###########################################
# Count merged peaks
###########################################

log_message("Counting merged peaks in output BED files.")

new_bed_files <- list.files(
  path = output_folder,
  pattern = "\\.bed$",
  full.names = TRUE
)

if (length(new_bed_files) == 0) {
  log_error(
    "No merged BED files were generated in: ",
    output_folder
  )
}

peak_counts <- vapply(
  new_bed_files,
  R.utils::countLines,
  numeric(1)
)

peak_count_table <- data.frame(
  peak_count = peak_counts,
  row.names = basename(new_bed_files)
)


###########################################
# Export results
###########################################

log_message("Writing peak count table: ", output_file)

write.table(
  peak_count_table,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = FALSE
)

log_message("Peak merging and counting completed.")