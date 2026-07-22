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

log_message("Starting ATAC ctrl/cKO annotation and overlap script")

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
# ChIPseeker options
###########################################

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)

###########################################
# Load inputs, outputs and parameters
###########################################

log_message("Loading Snakemake inputs, outputs and parameters")

ATAC_file <- snakemake@input[["ATAC"]]
TRIM28_file <- snakemake@input[["TRIM28"]]
FOXL2_file <- snakemake@input[["FOXL2"]]
H3K9me3_file <- snakemake@input[["H3K9me3"]]
genome <- snakemake@input[["genome"]]

promoter <- as.numeric(snakemake@params[["promoter"]])

output_table <- snakemake@output[["table"]]

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}

###########################################
# Load peak files
###########################################

log_message("Importing ATAC peaks")
ATAC <- read.csv(
  ATAC_file,
  header = TRUE,
  sep = "\t"
)

ATAC <- GenomicRanges::GRanges(ATAC$region)

log_message("Importing TRIM28 peaks")
TRIM28 <- rtracklayer::import(TRIM28_file)

log_message("Importing FOXL2 peaks")
FOXL2 <- rtracklayer::import(FOXL2_file)

log_message("Importing H3K9me3 peaks")
H3K9me3 <- rtracklayer::import(H3K9me3_file)

log_message("ATAC peaks: ", length(ATAC))
log_message("TRIM28 peaks: ", length(TRIM28))
log_message("FOXL2 peaks: ", length(FOXL2))
log_message("H3K9me3 peaks: ", length(H3K9me3))

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(genome)

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(genome)

###########################################
# Functions
###########################################

make_region_id <- function(regions) {
  paste0(
    GenomicRanges::seqnames(regions),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

get_overlap_status <- function(query_regions, subject_regions) {
  GenomicRanges::countOverlaps(
    query_regions,
    subject_regions,
    ignore.strand = TRUE
  ) > 0
}

annotate_ATAC_peaks <- function(
  regions,
  condition,
  TxDb,
  gene2symbol,
  promoter
) {
  log_message("Annotating ATAC peaks for condition: ", condition)

  if (length(regions) == 0) {
    return(data.frame(
      condition = character(),
      chromosome = character(),
      start = integer(),
      end = integer(),
      region = character(),
      annotation = character(),
      nearest_gene = character(),
      distance_to_TSS = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  annotation <- as.data.frame(
    ChIPseeker::annotatePeak(
      regions,
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

  annotation$geneId <- gene2symbol[
    annotation$geneId,
    "gene_name"
  ]

  data.frame(
    condition = condition,
    chromosome = as.character(GenomicRanges::seqnames(regions)),
    start = GenomicRanges::start(regions),
    end = GenomicRanges::end(regions),
    region = make_region_id(regions),
    annotation = annotation$annotation,
    nearest_gene = annotation$geneId,
    distance_to_TSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

make_ATAC_table <- function(
  regions,
  condition,
  TRIM28,
  FOXL2,
  H3K9me3,
  TxDb,
  gene2symbol,
  promoter
) {
  log_message("Computing overlaps for condition: ", condition)

  overlap_TRIM28 <- get_overlap_status(regions, TRIM28)
  overlap_FOXL2 <- get_overlap_status(regions, FOXL2)
  overlap_H3K9me3 <- get_overlap_status(regions, H3K9me3)

  log_message(condition, " peaks overlapping TRIM28: ", sum(overlap_TRIM28))
  log_message(condition, " peaks overlapping FOXL2: ", sum(overlap_FOXL2))
  log_message(condition, " peaks overlapping H3K9me3: ", sum(overlap_H3K9me3))

  annotation_table <- annotate_ATAC_peaks(
    regions = regions,
    condition = condition,
    TxDb = TxDb,
    gene2symbol = gene2symbol,
    promoter = promoter
  )

  overlap_table <- data.frame(
    condition = condition,
    region = make_region_id(regions),
    overlap_TRIM28 = ifelse(overlap_TRIM28, "YES", "NO"),
    overlap_FOXL2 = ifelse(overlap_FOXL2, "YES", "NO"),
    overlap_H3K9me3 = ifelse(overlap_H3K9me3, "YES", "NO"),
    stringsAsFactors = FALSE
  )

  final_table <- merge(
    overlap_table,
    annotation_table,
    by = c("condition", "region"),
    all.x = TRUE,
    sort = FALSE
  )

  final_table[, c(
    "condition",
    "chromosome",
    "start",
    "end",
    "region",
    "overlap_TRIM28",
    "overlap_FOXL2",
    "overlap_H3K9me3",
    "annotation",
    "nearest_gene",
    "distance_to_TSS"
  )]
}

###########################################
# Run analysis
###########################################

final_table <- make_ATAC_table(
  regions = ATAC,
  condition = "ctrl",
  TRIM28 = TRIM28,
  FOXL2 = FOXL2,
  H3K9me3 = H3K9me3,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)


###########################################
# Save output table
###########################################

log_message("Writing output table: ", output_table)

write.table(
  final_table,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

log_message("Analysis completed successfully")