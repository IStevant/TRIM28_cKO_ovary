source(".Rprofile")

suppressPackageStartupMessages({
  library("GenomicRanges")
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
raw_peak_folder <- snakemake@input[["raw_peaks"]]
domains_folder <- snakemake@input[["domains"]]
genome <- snakemake@input[["genome"]]

# Parameters
AB <- snakemake@params[["AB"]]
promoter <- snakemake@params[["promoter"]]
distance <- snakemake@params[["distance"]]
minimum_rep <- snakemake@params[["minimum_rep"]]

# Outputs
consensus_peak_table <- snakemake@output[["consensus_peak_table"]]
consensus_domain_table <- snakemake@output[["consensus_domain_table"]]

consensus_peak_bed <- snakemake@output[["consensus_peak_bed"]]
consensus_domain_bed <- snakemake@output[["consensus_domain_bed"]]

gtf_files <- snakemake@output[["gtf_files"]]


###########################################
# Validate inputs
###########################################

if (!dir.exists(raw_peak_folder)) {
  log_error(
    "Raw peak folder does not exist: ",
    raw_peak_folder
  )
}

if (!dir.exists(domains_folder)) {
  log_error(
    "Domain folder does not exist: ",
    domains_folder
  )
}

if (!file.exists(genome)) {
  log_error(
    "Genome annotation file does not exist: ",
    genome
  )
}

if (length(AB) != 1 || is.na(AB) || AB == "") {
  log_error("A unique antibody name must be provided.")
}

if (!AB %in% names(distance)) {
  log_error(
    "No domain-merging distance was defined for antibody: ",
    AB
  )
}

distance_between_peaks <- distance[[AB]]


###########################################
# Prepare genome annotation
###########################################

log_message("Preparing genome annotation: ", genome)

genome_gtf <- rtracklayer::import(
  genome
)

gene2symbol <- GenomicRanges::mcols(
  genome_gtf
)[, c("gene_id", "gene_name")]

gene2symbol <- unique(
  gene2symbol
)

rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- txdbmaker::makeTxDbFromGFF(
  genome
)


###########################################
# Functions
###########################################

#' Export annotated genomic regions
#'
#' @param gr Annotated GRanges object.
#' @param file Output table path.
#'
#' @return No return value.
export_table <- function(gr, file) {

  region <- paste0(
    seqnames(gr),
    ":",
    start(gr),
    "-",
    end(gr)
  )

  annotation_columns <- c(
    "annotation",
    "geneId",
    "distanceToTSS"
  )

  data <- data.frame(
    chromosome = seqnames(gr),
    start = start(gr),
    end = end(gr),
    region = region,
    as.data.frame(
      mcols(gr[, annotation_columns])
    ),
    stringsAsFactors = FALSE
  )

  write.table(
    data,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}


#' Generate consensus genomic regions
#'
#' Regions are retained when they are supported by at least
#' `min_overlap` replicate files. Adjacent retained blocks are then
#' merged according to `distance`.
#'
#' @param gr_list List of GRanges objects.
#' @param min_overlap Minimum number of supporting replicates.
#' @param distance Maximum distance used to merge retained blocks.
#'
#' @return An annotated GRanges object.
get_consensus <- function(
    gr_list,
    min_overlap,
    distance) {

  if (length(gr_list) == 0) {
    log_error(
      "Cannot generate consensus regions from an empty GRanges list."
    )
  }

  for (i in seq_along(gr_list)) {
    mcols(gr_list[[i]])$source_file <- paste0(
      "file_",
      i
    )
  }

  all_regions <- do.call(
    c,
    gr_list
  )

  blocks <- GenomicRanges::disjoin(
    all_regions
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

  if (is.null(dim(presence_matrix))) {
    presence_matrix <- matrix(
      presence_matrix,
      ncol = 1
    )
  }

  support_count <- rowSums(
    presence_matrix
  )

  kept_blocks <- blocks[
    support_count >= min_overlap
  ]

  if (length(kept_blocks) == 0) {
    log_error(
      "No region was supported by at least ",
      min_overlap,
      " replicate(s)."
    )
  }

  fused_regions <- GenomicRanges::reduce(
    kept_blocks,
    min.gapwidth = as.numeric(distance)
  )

  peak_anno <- ChIPseeker::as.GRanges(
    ChIPseeker::annotatePeak(
      fused_regions,
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


#' Export genomic regions in SAF format
#'
#' @param gr GRanges object.
#'
#' @return No return value.
make_SAF <- function(gr) {

  data <- as.data.frame(
    gr
  )

  genes <- paste0(
    seqnames(gr),
    ":",
    start(gr),
    "-",
    end(gr)
  )

  saf <- data.frame(
    genes,
    data
  )

  write.table(
    saf,
    file = gtf_files,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )
}


###########################################
# Generate consensus peaks
###########################################

log_message(
  "Searching for peak files associated with antibody: ",
  AB
)

peak_files <- list.files(
  path = raw_peak_folder,
  pattern = AB,
  full.names = TRUE
)

if (length(peak_files) == 0) {
  log_error(
    "No peak files found for antibody ",
    AB,
    " in: ",
    raw_peak_folder
  )
}

log_message(
  "Found ",
  length(peak_files),
  " peak file(s)."
)

peak_granges <- lapply(
  peak_files,
  function(file) {

    bed <- read.csv(
      file,
      header = FALSE,
      sep = "\t"
    )

    GenomicRanges::GRanges(
      seqnames = bed[[1]],
      ranges = IRanges(
        start = bed[[2]],
        end = bed[[3]]
      ),
      strand = "*"
    )
  }
)

log_message("Generating consensus peaks.")

consensus_peaks <- get_consensus(
  peak_granges,
  minimum_rep,
  1
)

log_message(
  "Consensus peaks identified: ",
  length(consensus_peaks)
)

log_message(
  "Exporting consensus peak BED file: ",
  consensus_peak_bed
)

rtracklayer::export(
  consensus_peaks,
  consensus_peak_bed,
  format = "BED"
)

log_message(
  "Exporting consensus peak table: ",
  consensus_peak_table
)

export_table(
  consensus_peaks,
  consensus_peak_table
)


###########################################
# Generate consensus domains
###########################################

log_message(
  "Searching for domain files associated with antibody: ",
  AB
)

domain_files <- list.files(
  path = domains_folder,
  pattern = AB,
  full.names = TRUE
)

if (length(domain_files) == 0) {
  log_error(
    "No domain files found for antibody ",
    AB,
    " in: ",
    domains_folder
  )
}

log_message(
  "Found ",
  length(domain_files),
  " domain file(s)."
)

domain_granges <- lapply(
  domain_files,
  rtracklayer::import
)

log_message(
  "Generating consensus domains using a merging distance of ",
  distance_between_peaks,
  " bp."
)

consensus_domains <- get_consensus(
  domain_granges,
  minimum_rep,
  distance_between_peaks
)

log_message(
  "Consensus domains identified: ",
  length(consensus_domains)
)

log_message(
  "Exporting consensus domain BED file: ",
  consensus_domain_bed
)

rtracklayer::export(
  consensus_domains,
  consensus_domain_bed,
  format = "BED"
)

log_message(
  "Exporting consensus domain table: ",
  consensus_domain_table
)

export_table(
  consensus_domains,
  consensus_domain_table
)


###########################################
# Export SAF file
###########################################

log_message(
  "Exporting consensus domains in SAF format: ",
  gtf_files
)

make_SAF(
  consensus_domains
)

log_message(
  "Consensus peak and domain generation completed."
)