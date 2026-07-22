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

log_message("Starting intergenic TE / ATAC cluster overlap script")

###########################################
# Libraries
###########################################

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("IRanges")
  library("S4Vectors")
})

###########################################
# Load inputs and outputs
###########################################

TE_up_file <- snakemake@input[["TE_up"]]
TE_down_file <- snakemake@input[["TE_down"]]
ATAC_clusters_file <- snakemake@input[["ATAC_clusters"]]
repeatMasker_file <- snakemake@input[["repeatMasker"]]

output_table <- snakemake@output[["table"]]
output_summary <- snakemake@output[["summary"]]

###########################################
# Functions
###########################################

make_region_id <- function(regions) {
  paste0(
    as.character(GenomicRanges::seqnames(regions)),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

import_TE_bed <- function(file, direction) {
  log_message("Importing ", direction, " TE BED file: ", file)

  bed <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = ""
  )

  if (ncol(bed) < 6) {
    log_error("TE BED file must contain at least 6 columns")
  }

  TE_gr <- GenomicRanges::GRanges(
    seqnames = bed[[1]],
    ranges = IRanges::IRanges(
      start = bed[[2]] + 1,
      end = bed[[3]]
    ),
    strand = bed[[6]]
  )

  TE_gr$TE_name <- bed[[4]]
  TE_gr$score <- bed[[5]]
  TE_gr$direction <- direction
  TE_gr$TE_region <- make_region_id(TE_gr)

  TE_info <- strsplit(
    as.character(TE_gr$TE_name),
    ":",
    fixed = TRUE
  )

  TE_gr$TE_family <- vapply(
    TE_info,
    function(x) ifelse(length(x) >= 1, x[1], NA_character_),
    character(1)
  )

  TE_gr$TE_subclass <- vapply(
    TE_info,
    function(x) ifelse(length(x) >= 2, x[2], NA_character_),
    character(1)
  )

  TE_gr$TE_class <- vapply(
    TE_info,
    function(x) ifelse(length(x) >= 3, x[3], NA_character_),
    character(1)
  )

  return(TE_gr)
}

import_ATAC_clusters <- function(file) {
  log_message("Importing ATAC cluster table: ", file)

  ATAC <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("region", "cluster") %in% colnames(ATAC))) {
    log_error("ATAC cluster table must contain 'region' and 'cluster' columns")
  }

  ATAC_gr <- GenomicRanges::GRanges(
    ATAC$region
  )

  ATAC_gr$cluster <- ATAC$cluster
  ATAC_gr$ATAC_region <- ATAC$region

  return(ATAC_gr)
}

import_repeatMasker_TEs <- function(file) {
  log_message("Importing RepeatMasker annotation: ", file)

  rmsk <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    check.names = FALSE
  )

  if (!all(c("repClass", "repFamily") %in% colnames(rmsk))) {
    log_error("RepeatMasker file must contain 'repClass' and 'repFamily' columns")
  }

  TE_classes <- c(
    "LINE",
    "SINE",
    "LTR",
    "DNA",
    "RC"
  )

  rmsk <- rmsk[
    rmsk$repClass %in% TE_classes,
    ,
    drop = FALSE
  ]

  if (nrow(rmsk) == 0) {
    log_error("No transposable elements retained after filtering RepeatMasker")
  }

  repeats <- GenomicRanges::makeGRangesFromDataFrame(
    rmsk,
    keep.extra.columns = TRUE
  )

  repeats$repeat_region <- make_region_id(repeats)

  log_message("RepeatMasker TE annotations retained: ", length(repeats))

  return(repeats)
}

make_overlap_table <- function(TE_gr, ATAC_gr) {
  log_message("Finding overlaps between DEG TE regions and ATAC clusters")

  hits <- GenomicRanges::findOverlaps(
    TE_gr,
    ATAC_gr,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    return(data.frame())
  }

  query_idx <- S4Vectors::queryHits(hits)
  subject_idx <- S4Vectors::subjectHits(hits)

  overlap_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(TE_gr))[query_idx],
    start = GenomicRanges::start(TE_gr)[query_idx],
    end = GenomicRanges::end(TE_gr)[query_idx],
    TE_name = TE_gr$TE_name[query_idx],
    TE_family = TE_gr$TE_family[query_idx],
    TE_subclass = TE_gr$TE_subclass[query_idx],
    TE_class = TE_gr$TE_class[query_idx],
    score = TE_gr$score[query_idx],
    strand = as.character(GenomicRanges::strand(TE_gr))[query_idx],
    direction = TE_gr$direction[query_idx],
    ATAC_cluster = ATAC_gr$cluster[subject_idx],
    ATAC_region = ATAC_gr$ATAC_region[subject_idx],
    stringsAsFactors = FALSE
  )

  overlap_table <- unique(overlap_table)

  overlap_table <- overlap_table[
    order(
      overlap_table$direction,
      overlap_table$ATAC_cluster,
      overlap_table$chromosome,
      overlap_table$start
    ),
    ,
    drop = FALSE
  ]

  log_message("Number of DEG TE / ATAC cluster overlaps: ", nrow(overlap_table))

  return(overlap_table)
}

count_all_repeatMasker_TEs_in_ATAC_clusters <- function(
  repeatMasker_TEs,
  ATAC_gr
) {
  log_message("Counting all RepeatMasker TEs overlapping each ATAC cluster")

  hits <- GenomicRanges::findOverlaps(
    repeatMasker_TEs,
    ATAC_gr,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {
    return(data.frame(
      ATAC_cluster = unique(ATAC_gr$cluster),
      n_total_repeatMasker_TEs_in_cluster = 0,
      stringsAsFactors = FALSE
    ))
  }

  repeat_idx <- S4Vectors::queryHits(hits)
  atac_idx <- S4Vectors::subjectHits(hits)

  overlap_df <- data.frame(
    repeat_region = repeatMasker_TEs$repeat_region[repeat_idx],
    ATAC_cluster = ATAC_gr$cluster[atac_idx],
    stringsAsFactors = FALSE
  )

  overlap_df <- unique(overlap_df)

  count_table <- aggregate(
    repeat_region ~ ATAC_cluster,
    data = overlap_df,
    FUN = function(x) length(unique(x))
  )

  colnames(count_table) <- c(
    "ATAC_cluster",
    "n_total_repeatMasker_TEs_in_cluster"
  )

  all_clusters <- data.frame(
    ATAC_cluster = unique(ATAC_gr$cluster),
    stringsAsFactors = FALSE
  )

  count_table <- merge(
    all_clusters,
    count_table,
    by = "ATAC_cluster",
    all.x = TRUE,
    sort = FALSE
  )

  count_table$n_total_repeatMasker_TEs_in_cluster[
    is.na(count_table$n_total_repeatMasker_TEs_in_cluster)
  ] <- 0

  return(count_table)
}

make_summary_table <- function(
  overlap_table,
  repeatMasker_TE_counts
) {
  if (nrow(overlap_table) == 0) {
    summary_table <- data.frame(
      direction = character(),
      ATAC_cluster = character(),
      n_DEG_TE_overlaps = integer(),
      n_unique_DEG_TE = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    split_key <- paste(
      overlap_table$direction,
      overlap_table$ATAC_cluster,
      sep = "___"
    )

    summary_list <- lapply(
      unique(split_key),
      function(key) {
        subset_table <- overlap_table[
          split_key == key,
          ,
          drop = FALSE
        ]

        data.frame(
          direction = unique(subset_table$direction),
          ATAC_cluster = unique(subset_table$ATAC_cluster),
          n_DEG_TE_overlaps = nrow(subset_table),
          n_unique_DEG_TE = length(unique(paste0(
            subset_table$chromosome,
            ":",
            subset_table$start,
            "-",
            subset_table$end,
            ":",
            subset_table$TE_name
          ))),
          stringsAsFactors = FALSE
        )
      }
    )

    summary_table <- do.call(
      rbind,
      summary_list
    )
  }

  summary_table <- merge(
    summary_table,
    repeatMasker_TE_counts,
    by = "ATAC_cluster",
    all.x = TRUE,
    sort = FALSE
  )

  summary_table <- summary_table[
    order(
      summary_table$direction,
      summary_table$ATAC_cluster
    ),
    ,
    drop = FALSE
  ]

  return(summary_table)
}

###########################################
# Run analysis
###########################################

TE_up <- import_TE_bed(
  file = TE_up_file,
  direction = "UP"
)

TE_down <- import_TE_bed(
  file = TE_down_file,
  direction = "DOWN"
)

TE_all <- c(
  TE_up,
  TE_down
)

log_message("Total UP DEG TE regions: ", length(TE_up))
log_message("Total DOWN DEG TE regions: ", length(TE_down))

ATAC_clusters <- import_ATAC_clusters(
  file = ATAC_clusters_file
)

repeatMasker_TEs <- import_repeatMasker_TEs(
  file = repeatMasker_file
)

overlap_table <- make_overlap_table(
  TE_gr = TE_all,
  ATAC_gr = ATAC_clusters
)

repeatMasker_TE_counts <- count_all_repeatMasker_TEs_in_ATAC_clusters(
  repeatMasker_TEs = repeatMasker_TEs,
  ATAC_gr = ATAC_clusters
)

summary_table <- make_summary_table(
  overlap_table = overlap_table,
  repeatMasker_TE_counts = repeatMasker_TE_counts
)

###########################################
# Save output files
###########################################

write.table(
  overlap_table,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  summary_table,
  file = output_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

log_message("Analysis completed successfully")