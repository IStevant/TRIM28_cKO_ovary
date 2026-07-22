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
condition <- snakemake@params[["condition"]]
promoter <- snakemake@params[["promoter"]]
distance <- snakemake@params[["distance"]]
minimum_rep <- snakemake@params[["minimum_rep"]]
tmp_folder <- snakemake@params[["tmp"]]

# Outputs
consensus_peak_table <- snakemake@output[["consensus_peak_table"]]
consensus_domain_table <- snakemake@output[["consensus_domain_table"]]

consensus_peak_bed <- snakemake@output[["consensus_peak_bed"]]
consensus_domain_bed <- snakemake@output[["consensus_domain_bed"]]


###########################################
# Validate inputs
###########################################

for (folder in c(raw_peak_folder, domains_folder)) {
  if (!dir.exists(folder)) {
    log_error("Directory does not exist: ", folder)
  }
}

if (!file.exists(genome)) {
  log_error("Genome annotation file does not exist: ", genome)
}


###########################################
# Prepare genome annotation
###########################################

log_message("Preparing genome annotation.")

# Get antibody from condition name
AB <- unique(
  sapply(
    strsplit(condition, "_"),
    `[`,
    2
  )
)

if (length(AB) != 1) {
  log_error(
    "Unable to determine a unique antibody from condition: ",
    condition
  )
}

distance_between_peaks <- distance[[AB]]

genome_gtf <- rtracklayer::import(genome)

gene2symbol <- unique(
  GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
)

rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- txdbmaker::makeTxDbFromGFF(genome)


###########################################
# Functions
###########################################

# Export annotated GRanges as a tab-separated table.
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


# Generate consensus regions supported by a minimum number of replicates.
get_consensus <- function(
    gr_list,
    min_overlap,
    distance) {

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

  blocks <- GenomicRanges::disjoin(all_regions)

  presence_matrix <- sapply(
    gr_list,
    function(gr) {
      !is.na(
        findOverlaps(
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

  fused_regions <- GenomicRanges::reduce(
    kept_blocks,
    min.gapwidth = as.numeric(distance)
  )

  annotations <- ChIPseeker::as.GRanges(
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

  annotations$geneId <- gene2symbol[
    annotations$geneId,
    "gene_name"
  ]

  return(annotations)
}


###########################################
# Consensus peaks
###########################################

log_message("Generating consensus peaks.")

peak_files <- list.files(
  path = raw_peak_folder,
  pattern = condition,
  full.names = TRUE
)

if (length(peak_files) == 0) {
  log_error(
    "No peak files found for condition: ",
    condition
  )
}

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

consensus_peaks <- get_consensus(
  peak_granges,
  minimum_rep,
  1
)

log_message(
  "Consensus peaks identified: ",
  length(consensus_peaks)
)

rtracklayer::export(
  consensus_peaks,
  consensus_peak_bed,
  format = "BED"
)

export_table(
  consensus_peaks,
  consensus_peak_table
)


###########################################
# Consensus domains
###########################################

log_message("Generating consensus domains.")

domain_files <- list.files(
  path = domains_folder,
  pattern = condition,
  full.names = TRUE
)

if (length(domain_files) == 0) {
  log_error(
    "No domain files found for condition: ",
    condition
  )
}

domain_granges <- lapply(
  domain_files,
  rtracklayer::import
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

rtracklayer::export(
  consensus_domains,
  consensus_domain_bed,
  format = "BED"
)

export_table(
  consensus_domains,
  consensus_domain_table
)

log_message("Consensus peak and domain generation completed.")