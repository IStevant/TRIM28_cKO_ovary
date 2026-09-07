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

log_message("Starting ChIP consensus peak annotation script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("IRanges")
  library("rtracklayer")
  library("txdbmaker")
  library("ChIPseeker")
})

###########################################
# Load inputs, outputs and parameters
###########################################

log_message("Loading Snakemake inputs, outputs and parameters")

raw_peak_folder <- snakemake@input[["raw_peaks"]]
genome <- snakemake@input[["genome"]]

condition <- snakemake@params[["condition"]]
promoter <- as.numeric(snakemake@params[["promoter"]])
minimum_rep <- as.numeric(snakemake@params[["minimum_rep"]])

consensus_peak_table <- snakemake@output[["consensus_peak_table"]]
consensus_peak_bed <- snakemake@output[["consensus_peak_bed"]]

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}

if (is.na(minimum_rep)) {
  log_error("Parameter 'minimum_rep' must be numeric")
}

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(
  genome
)

log_message("Preparing gene ID to symbol mapping")

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(
  genome
)

###########################################
# Functions
###########################################

export_table <- function(
  gr,
  file
) {
  log_message("Writing consensus peak annotation table: ", file)

  region <- paste0(
    GenomicRanges::seqnames(gr),
    ":",
    GenomicRanges::start(gr),
    "-",
    GenomicRanges::end(gr)
  )

  annotation_columns <- c(
    "annotation",
    "geneId",
    "distanceToTSS"
  )

  output_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    region = region,
    as.data.frame(
      GenomicRanges::mcols(gr)[, annotation_columns]
    ),
    stringsAsFactors = FALSE
  )

  colnames(output_table) <- c(
    "chromosome",
    "start",
    "end",
    "region",
    "annotation",
    "nearest_gene",
    "distanceToTSS"
  )

  write.table(
    output_table,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

read_peak_file <- function(
  file
) {
  log_message("Importing peak file: ", basename(file))

  peak_gr <- rtracklayer::import(
    file
  )

  peak_gr <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(peak_gr),
    ranges = GenomicRanges::ranges(peak_gr),
    strand = "*"
  )

  return(peak_gr)
}

get_consensus <- function(
  gr_list,
  min_overlap
) {
  log_message("Building consensus peaks")

  if (length(gr_list) == 0) {
    log_error("No peak files were provided")
  }

  for (i in seq_along(gr_list)) {
    GenomicRanges::mcols(gr_list[[i]])$source_file <- paste0("file_", i)
  }

  all_gr <- do.call(
    c,
    gr_list
  )

  blocks <- GenomicRanges::disjoin(
    all_gr
  )

  presence_matrix <- sapply(
    gr_list,
    function(gr) {
      !is.na(
        GenomicRanges::findOverlaps(
          blocks,
          gr,
          select = "first"
        )
      )
    }
  )

  support_count <- rowSums(
    presence_matrix
  )

  kept_blocks <- blocks[
    support_count >= min_overlap
  ]

  if (length(kept_blocks) == 0) {
    log_error("No consensus peak passed the minimum replicate threshold")
  }

  consensus_regions <- GenomicRanges::reduce(
    kept_blocks,
    min.gapwidth = 1
  )

  log_message("Consensus peaks retained: ", length(consensus_regions))

  log_message("Annotating consensus peaks")

  peak_anno <- ChIPseeker::as.GRanges(
    ChIPseeker::annotatePeak(
      consensus_regions,
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

  peak_anno$geneId <- gene2symbol[
    peak_anno$geneId,
    "gene_name"
  ]

  return(peak_anno)
}

###########################################
# Build ChIP consensus peaks
###########################################

log_message("Searching ChIP peak files for condition: ", condition)

if (any(grep("FOXL2", condition))){
  peak_pattern <- paste0(
    "^",
    condition,
    "_peaks\\.narrowPeak$"
  )
} else {
  peak_pattern <- paste0(
    "^",
    condition,
    "_Rep[0-9]+_peaks\\.narrowPeak$"
  )
}



bed_files <- list.files(
  path = raw_peak_folder,
  pattern = peak_pattern,
  full.names = TRUE
)

if (length(bed_files) == 0) {
  log_error(
    "No peak files found for condition: ",
    condition,
    " with expected pattern: ",
    peak_pattern
  )
}

log_message(
  "Peak files found: ",
  paste(basename(bed_files), collapse = ", ")
)

gr_list <- lapply(
  bed_files,
  read_peak_file
)

consensus <- get_consensus(
  gr_list = gr_list,
  min_overlap = minimum_rep
)

###########################################
# Save outputs
###########################################

log_message("Writing consensus peak BED: ", consensus_peak_bed)

rtracklayer::export(
  consensus,
  consensus_peak_bed,
  format = "BED"
)

export_table(
  gr = consensus,
  file = consensus_peak_table
)

log_message("Analysis completed successfully")