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

log_message("Starting BED peak annotation script")

###########################################
# Libraries
###########################################

log_message("Loading libraries")

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("rtracklayer")
  library("txdbmaker")
  library("ChIPseeker")
})

###########################################
# Load inputs, outputs and parameters
###########################################

log_message("Loading Snakemake inputs, outputs and parameters")

bed_file <- snakemake@input[["bed"]]
genome <- snakemake@input[["genome"]]

promoter <- as.numeric(snakemake@params[["promoter"]])

output_table <- snakemake@output[["table"]]

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
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

make_region_id <- function(
  regions
) {
  paste0(
    GenomicRanges::seqnames(regions),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

annotate_bed_peaks <- function(
  bed_file,
  TxDb,
  gene2symbol,
  promoter
) {
  log_message("Importing BED file: ", bed_file)

  peaks <- rtracklayer::import(
    bed_file
  )

  if (length(peaks) == 0) {
    log_error("BED file contains no peaks")
  }

  log_message("Number of peaks to annotate: ", length(peaks))

  log_message("Annotating peaks with ChIPseeker")

  peak_anno <- ChIPseeker::as.GRanges(
    ChIPseeker::annotatePeak(
      peaks,
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

export_annotation_table <- function(
  annotated_peaks,
  output_file
) {
  log_message("Preparing annotation table")

  region <- make_region_id(
    annotated_peaks
  )

  annotation_columns <- c(
    "annotation",
    "geneId",
    "distanceToTSS"
  )

  output_table <- data.frame(
    chromosome = as.character(
      GenomicRanges::seqnames(annotated_peaks)
    ),
    start = GenomicRanges::start(annotated_peaks),
    end = GenomicRanges::end(annotated_peaks),
    region = region,
    as.data.frame(
      GenomicRanges::mcols(
        annotated_peaks
      )[, annotation_columns]
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

  log_message("Writing annotated peak table: ", output_file)

  write.table(
    output_table,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

###########################################
# Run annotation
###########################################

annotated_peaks <- annotate_bed_peaks(
  bed_file = bed_file,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)

export_annotation_table(
  annotated_peaks = annotated_peaks,
  output_file = output_table
)

log_message("Analysis completed successfully")