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

raw_peak_folder <- snakemake@input[["raw_peaks"]]

promoter <- snakemake@params[["promoter"]]
minimum_rep <- snakemake@params[["minimum_rep"]]

anno_folder <- snakemake@output[["anno_folder"]]

dir.create(
  anno_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


###########################################
# Input files and experimental conditions
###########################################

log_message("Searching for peak files.")

raw_bed_files <- list.files(
  path = raw_peak_folder,
  pattern = "\\.mLb\\.clN_peaks\\.broadPeak$"
)

if (length(raw_bed_files) == 0) {
  log_error("No broadPeak files found in: ", raw_peak_folder)
}

# Extract condition names from peak filenames.
conditions <- unique(
  gsub(
    "_R....mLb\\.clN_peaks\\.broadPeak$",
    "",
    raw_bed_files
  )
)

log_message(
  "Detected conditions: ",
  paste(conditions, collapse = ", ")
)


###########################################
# Functions
###########################################

# Export annotated genomic regions as a tab-separated table.
export_table <- function(gr, file) {

  log_message("Writing annotation table: ", basename(file))

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

  output_table <- data.frame(
    chromosome = seqnames(gr),
    start = start(gr),
    end = end(gr),
    region = region,
    as.data.frame(mcols(gr[, annotation_columns])),
    stringsAsFactors = FALSE
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


# Generate consensus regions supported by at least `min_overlap`
# replicate peak sets, fuse adjacent blocks, and annotate the regions.
get_consensus <- function(gr_list, min_overlap, distance) {

  log_message(
    "Generating consensus peaks (minimum overlap = ",
    min_overlap,
    ", merge distance = ",
    distance,
    ")."
  )

  for (i in seq_along(gr_list)) {
    mcols(gr_list[[i]])$source_file <- paste0("file_", i)
  }

  all_gr <- do.call(c, gr_list)

  blocks <- disjoin(all_gr)

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

  support_count <- rowSums(presence_matrix)

  kept_blocks <- blocks[support_count >= min_overlap]

  fused_regions <- reduce(
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

  log_message(
    "Annotated ",
    length(peak_anno),
    " consensus peaks."
  )

  return(peak_anno)
}


###########################################
# Genome annotation
###########################################

log_message("Loading GENCODE genome annotation.")

gtf_url <- paste0(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/",
  "Gencode_mouse/release_M25/",
  "gencode.vM25.annotation.gtf.gz"
)

alt_gtf_url <- paste0(
  "http://ftp.cbi.pku.edu.cn/pub/mirror/GENCODE/",
  "Gencode_mouse/release_M25/",
  "gencode.vM25.annotation.gtf.gz"
)

genome_gtf <- tryCatch(
  {
    rtracklayer::import(gtf_url)
  },
  error = function(e) {
    log_message("Primary GENCODE mirror unavailable. Trying alternative mirror.")
    rtracklayer::import(alt_gtf_url)
  }
)

gene2symbol <- GenomicRanges::mcols(genome_gtf)[
  ,
  c("gene_id", "gene_name")
]

gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- txdbmaker::makeTxDbFromGFF(gtf_url)


###########################################
# Generate consensus peaks
###########################################

for (condition in conditions) {

  log_message("Processing condition: ", condition)

  bed_files <- list.files(
    path = raw_peak_folder,
    pattern = paste0(
      condition,
      ".*\\.mLb\\.clN_peaks\\.broadPeak$"
    ),
    full.names = TRUE
  )

  if (length(bed_files) == 0) {
    log_error("No peak files found for condition: ", condition)
  }

  gr_list <- lapply(
    bed_files,
    function(file) {

      log_message("Reading ", basename(file))

      bed <- read.csv(
        file,
        header = FALSE,
        sep = "\t"
      )

      GRanges(
        seqnames = bed[[1]],
        ranges = IRanges(
          start = bed[[2]],
          end = bed[[3]]
        ),
        strand = "*"
      )
    }
  )

  consensus <- get_consensus(
    gr_list,
    minimum_rep,
    1
  )

  export_table(
    consensus,
    file.path(
      anno_folder,
      paste0(
        condition,
        "_consensus_peaks.txt"
      )
    )
  )
}

log_message("Consensus peak generation completed.")