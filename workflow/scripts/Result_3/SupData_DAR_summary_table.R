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

log_message("Starting ATAC DAR overlap summary script")

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

DAR_file <- snakemake@input[["ATAC"]]
TRIM28_file <- snakemake@input[["TRIM28"]]
FOXL2_file <- snakemake@input[["FOXL2"]]
ATAC_Sertoli_file <- snakemake@input[["ATAC_Sertoli"]]
ATAC_Granulosa_file <- snakemake@input[["ATAC_Granulosa"]]
H3K9me3_ctrl_file <- snakemake@input[["H3K9me3_ctrl"]]
H3K9me3_cKO_file <- snakemake@input[["H3K9me3_cKO"]]
genome <- snakemake@input[["genome"]]

promoter <- as.numeric(snakemake@params[["promoter"]])
output_table <- snakemake@output[["table"]]

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}

###########################################
# Load ATAC DAR regions
###########################################

log_message("Loading ATAC DAR table")

DAR <- read.csv(
  DAR_file,
  header = TRUE,
  check.names = FALSE,
  sep = "\t"
)

if (nrow(DAR) == 0) {
  log_error("ATAC DAR table contains no regions")
}

if (!"cluster" %in% colnames(DAR)) {
  log_error("Column 'cluster' not found in ATAC DAR table")
}

DAR_regions <- GenomicRanges::GRanges(
  DAR$region
)

names(DAR_regions) <- rownames(DAR)

log_message("Number of ATAC DAR regions: ", length(DAR_regions))

###########################################
# Load overlap datasets
###########################################

log_message("Importing TRIM28 peaks")
TRIM28_peaks <- rtracklayer::import(TRIM28_file)

log_message("Importing FOXL2 peaks")
FOXL2_peaks <- rtracklayer::import(FOXL2_file)

log_message("Importing Sertoli-biased ATAC peaks")
ATAC_Sertoli_peaks <- rtracklayer::import(ATAC_Sertoli_file)

log_message("Importing Granulosa-biased ATAC peaks")
ATAC_Granulosa_peaks <- rtracklayer::import(ATAC_Granulosa_file)

log_message("Importing H3K9me3 ctrl peaks")
H3K9me3_ctrl_peaks <- rtracklayer::import(H3K9me3_ctrl_file)

log_message("Importing H3K9me3 cKO peaks")
H3K9me3_cKO_peaks <- rtracklayer::import(H3K9me3_cKO_file)

log_message("TRIM28 peaks: ", length(TRIM28_peaks))
log_message("FOXL2 peaks: ", length(FOXL2_peaks))
log_message("Sertoli-biased ATAC peaks: ", length(ATAC_Sertoli_peaks))
log_message("Granulosa-biased ATAC peaks: ", length(ATAC_Granulosa_peaks))
log_message("H3K9me3 ctrl peaks: ", length(H3K9me3_ctrl_peaks))
log_message("H3K9me3 cKO peaks: ", length(H3K9me3_cKO_peaks))

###########################################
# Prepare genome annotation
###########################################

log_message("Loading genome annotation")

genome_gtf <- rtracklayer::import(genome)

log_message("Preparing gene ID to symbol mapping")

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

annotate_regions <- function(
  regions,
  TxDb,
  gene2symbol,
  promoter
) {
  log_message("Annotating ATAC DAR regions")

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
    region = make_region_id(regions),
    annotation = annotation$annotation,
    nearest_gene = annotation$geneId,
    distance_to_TSS = annotation$distanceToTSS,
    stringsAsFactors = FALSE
  )
}

make_ATAC_bias_category <- function(
  overlap_sertoli,
  overlap_granulosa
) {
  category <- rep(
    "None",
    length(overlap_sertoli)
  )

  category[overlap_sertoli] <- "Sertoli-biased ATAC"
  category[!overlap_sertoli & overlap_granulosa] <- "Granulosa-biased ATAC"

  category
}

###########################################
# Compute overlaps
###########################################

log_message("Computing overlaps")

overlap_TRIM28 <- get_overlap_status(DAR_regions, TRIM28_peaks)
overlap_FOXL2 <- get_overlap_status(DAR_regions, FOXL2_peaks)
overlap_ATAC_Sertoli <- get_overlap_status(DAR_regions, ATAC_Sertoli_peaks)
overlap_ATAC_Granulosa <- get_overlap_status(DAR_regions, ATAC_Granulosa_peaks)
overlap_H3K9me3_ctrl <- get_overlap_status(DAR_regions, H3K9me3_ctrl_peaks)
overlap_H3K9me3_cKO <- get_overlap_status(DAR_regions, H3K9me3_cKO_peaks)

ATAC_bias_category <- make_ATAC_bias_category(
  overlap_sertoli = overlap_ATAC_Sertoli,
  overlap_granulosa = overlap_ATAC_Granulosa
)

log_message("Regions overlapping TRIM28: ", sum(overlap_TRIM28))
log_message("Regions overlapping FOXL2: ", sum(overlap_FOXL2))
log_message("Regions overlapping Sertoli-biased ATAC: ", sum(overlap_ATAC_Sertoli))
log_message("Regions overlapping Granulosa-biased ATAC: ", sum(overlap_ATAC_Granulosa))
log_message("Regions overlapping H3K9me3 ctrl: ", sum(overlap_H3K9me3_ctrl))
log_message("Regions overlapping H3K9me3 cKO: ", sum(overlap_H3K9me3_cKO))

###########################################
# Annotate regions
###########################################

annotation_table <- annotate_regions(
  regions = DAR_regions,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)

###########################################
# Build output table
###########################################

log_message("Building final summary table")

region_id <- make_region_id(DAR_regions)

output_table_df <- data.frame(
  chromosome = as.character(GenomicRanges::seqnames(DAR_regions)),
  start = GenomicRanges::start(DAR_regions),
  end = GenomicRanges::end(DAR_regions),
  region = region_id,
  DAR_cluster = DAR$cluster,
  overlap_TRIM28 = ifelse(overlap_TRIM28, "YES", "NO"),
  overlap_FOXL2 = ifelse(overlap_FOXL2, "YES", "NO"),
  overlap_H3K9me3_ctrl = ifelse(overlap_H3K9me3_ctrl, "YES", "NO"),
  overlap_H3K9me3_cKO = ifelse(overlap_H3K9me3_cKO, "YES", "NO"),
  overlap_ATAC_Sertoli_biased = ifelse(overlap_ATAC_Sertoli, "YES", "NO"),
  overlap_ATAC_Granulosa_biased = ifelse(overlap_ATAC_Granulosa, "YES", "NO"),
  stringsAsFactors = FALSE
)

output_table_df <- merge(
  output_table_df,
  annotation_table,
  by = "region",
  all.x = TRUE,
  sort = FALSE
)

output_table_df <- output_table_df[, c(
  "chromosome",
  "start",
  "end",
  "region",
  "DAR_cluster",
  "overlap_TRIM28",
  "overlap_FOXL2",
  "overlap_H3K9me3_ctrl",
  "overlap_H3K9me3_cKO",
  "overlap_ATAC_Sertoli_biased",
  "overlap_ATAC_Granulosa_biased",
  "annotation",
  "nearest_gene",
  "distance_to_TSS"
)]

###########################################
# Save output table
###########################################

log_message("Writing output table: ", output_table)

write.table(
  output_table_df,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

log_message("Analysis completed successfully")