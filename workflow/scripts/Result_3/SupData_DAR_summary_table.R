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
  library("S4Vectors")
  library("rtracklayer")
  library("txdbmaker")
  library("GenomicFeatures")
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
SOX9_file <- snakemake@input[["SOX9"]]
DMRT1_file <- snakemake@input[["DMRT1"]]

ATAC_Sertoli_file <- snakemake@input[["ATAC_Sertoli"]]
ATAC_Granulosa_file <- snakemake@input[["ATAC_Granulosa"]]

H3K9me3_ctrl_file <- snakemake@input[["H3K9me3_ctrl"]]
H3K9me3_cKO_file <- snakemake@input[["H3K9me3_cKO"]]

expressed_genes_file <- snakemake@input[["expressed_genes"]]
DEG_8weeks_file <- snakemake@input[["DEG_8weeks"]]
DEG_7months_file <- snakemake@input[["DEG_7months"]]
sex_bias_file <- snakemake@input[["sex_bias"]]

genome <- snakemake@input[["genome"]]

promoter <- as.numeric(
  snakemake@params[["promoter"]]
)

nearest_gene_max_distance <- 50000

output_table <- snakemake@output[["table"]]

if (is.na(promoter)) {
  log_error("Parameter 'promoter' must be numeric")
}


###########################################
# Functions
###########################################

#' Generate genomic region IDs
#'
#' @param regions GRanges object.
#' @return Character vector formatted as chr:start-end.
make_region_id <- function(regions) {
  paste0(
    GenomicRanges::seqnames(regions),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}


#' Determine whether query regions overlap subject regions
#'
#' @param query_regions Query GRanges object.
#' @param subject_regions Subject GRanges object.
#' @return Logical vector.
get_overlap_status <- function(
  query_regions,
  subject_regions
) {
  GenomicRanges::countOverlaps(
    query_regions,
    subject_regions,
    ignore.strand = TRUE
  ) > 0
}


#' Load DEG direction table
#'
#' The input table must contain:
#'   - Gene
#'   - DEG: "up" or "down"
#'
#' @param file Path to DEG table.
#' @return Dataframe with columns gene and direction.
load_DEG_direction <- function(file) {
  data <- read.csv(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!all(c("Gene", "DEG") %in% colnames(data))) {
    log_error(
      "DEG table must contain columns 'Gene' and 'DEG': ",
      file
    )
  }

  output <- data.frame(
    gene = data$Gene,
    direction = tolower(data$DEG),
    stringsAsFactors = FALSE
  )

  output <- output[
    !is.na(output$gene) &
      output$gene != "",
    ,
    drop = FALSE
  ]

  invalid_direction <- !output$direction %in% c(
    "up",
    "down"
  )

  if (any(invalid_direction)) {
    log_error(
      "Column 'DEG' must contain only 'up' or 'down' in file: ",
      file
    )
  }

  output$direction <- ifelse(
    output$direction == "up",
    "Up",
    "Down"
  )

  output <- output[
    !duplicated(output$gene),
    ,
    drop = FALSE
  ]

  return(output)
}


#' Annotate regions using ChIPseeker
#'
#' @param regions GRanges object.
#' @param TxDb TxDb annotation.
#' @param gene2symbol Gene ID to gene symbol mapping.
#' @param promoter Promoter size.
#' @return Annotation dataframe.
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


#' Find nearest expressed gene
#'
#' If a region overlaps an expressed gene body, that gene is returned.
#' Otherwise, the nearest expressed gene based on TSS distance is returned
#' if it lies within max_distance.
#'
#' @param regions Query GRanges object.
#' @param gene_ranges GRanges object containing expressed gene bodies.
#' @param gene_TSS GRanges object containing expressed gene TSS.
#' @param max_distance Maximum accepted distance in bp.
#' @return Character vector of nearest expressed gene symbols.
get_nearest_expressed_gene <- function(
  regions,
  gene_ranges,
  gene_TSS,
  max_distance = 50000
) {
  log_message(
    "Finding overlapping or nearest expressed gene within +/- ",
    max_distance,
    " bp"
  )

  nearest_gene <- rep(
    NA_character_,
    length(regions)
  )

  if (length(gene_ranges) == 0) {
    log_error(
      "No expressed genes found in genome annotation"
    )
  }

  ###########################################
  # First assign overlapping gene bodies
  ###########################################

  overlap_hits <- GenomicRanges::findOverlaps(
    regions,
    gene_ranges,
    ignore.strand = TRUE
  )

  if (length(overlap_hits) > 0) {
    overlap_query <- S4Vectors::queryHits(
      overlap_hits
    )

    overlap_subject <- S4Vectors::subjectHits(
      overlap_hits
    )

    overlapping_queries <- unique(
      overlap_query
    )

    for (query_id in overlapping_queries) {
      subject_ids <- overlap_subject[
        overlap_query == query_id
      ]

      if (length(subject_ids) == 1) {
        selected_subject <- subject_ids
      } else {
        tss_distances <- GenomicRanges::distance(
          regions[query_id],
          gene_TSS[subject_ids],
          ignore.strand = TRUE
        )

        selected_subject <- subject_ids[
          which.min(tss_distances)
        ]
      }

      nearest_gene[
        query_id
      ] <- GenomicRanges::mcols(
        gene_ranges
      )$gene_symbol[
        selected_subject
      ]
    }
  }

  log_message(
    "Regions overlapping an expressed gene body: ",
    sum(!is.na(nearest_gene))
  )

  ###########################################
  # Assign nearest TSS to remaining regions
  ###########################################

  remaining_regions <- which(
    is.na(nearest_gene)
  )

  if (length(remaining_regions) > 0) {
    nearest_hits <- GenomicRanges::distanceToNearest(
      regions[remaining_regions],
      gene_TSS,
      ignore.strand = TRUE
    )

    query_index <- S4Vectors::queryHits(
      nearest_hits
    )

    subject_index <- S4Vectors::subjectHits(
      nearest_hits
    )

    distances <- S4Vectors::mcols(
      nearest_hits
    )$distance

    valid_hits <- distances <= max_distance

    original_query_index <- remaining_regions[
      query_index[valid_hits]
    ]

    nearest_gene[
      original_query_index
    ] <- GenomicRanges::mcols(
      gene_TSS
    )$gene_symbol[
      subject_index[valid_hits]
    ]
  }

  return(nearest_gene)
}


#' Generate ATAC sex-bias category
#'
#' @param overlap_sertoli Logical vector.
#' @param overlap_granulosa Logical vector.
#' @return Character vector with Sertoli, Granulosa or None.
make_ATAC_bias_category <- function(
  overlap_sertoli,
  overlap_granulosa
) {
  category <- rep(
    "None",
    length(overlap_sertoli)
  )

  category[
    overlap_sertoli
  ] <- "Sertoli"

  category[
    !overlap_sertoli &
      overlap_granulosa
  ] <- "Granulosa"

  return(category)
}


#' Match gene annotations
#'
#' @param genes Vector of gene symbols.
#' @param annotation_table Annotation dataframe.
#' @param value_column Column containing values to retrieve.
#' @param default Value used when no match is found.
#' @return Character vector with the same length as genes.
match_gene_annotation <- function(
  genes,
  annotation_table,
  value_column,
  default = "None"
) {
  values <- rep(
    default,
    length(genes)
  )

  valid_genes <- !is.na(genes) &
    genes != ""

  matched_index <- match(
    genes[valid_genes],
    annotation_table$gene
  )

  matched <- !is.na(matched_index)

  values[
    which(valid_genes)[matched]
  ] <- as.character(
    annotation_table[
      matched_index[matched],
      value_column
    ]
  )

  values[
    is.na(values) |
      values == ""
  ] <- default

  return(values)
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
  log_error(
    "ATAC DAR table contains no regions"
  )
}

if (!"cluster" %in% colnames(DAR)) {
  log_error(
    "Column 'cluster' not found in ATAC DAR table"
  )
}

if (!"region" %in% colnames(DAR)) {
  log_error(
    "Column 'region' not found in ATAC DAR table"
  )
}

DAR_regions <- GenomicRanges::GRanges(
  DAR$region
)

names(DAR_regions) <- DAR$region

log_message(
  "Number of ATAC DAR regions: ",
  length(DAR_regions)
)


###########################################
# Load overlap datasets
###########################################

log_message("Importing TRIM28 peaks")

TRIM28_peaks <- rtracklayer::import(
  TRIM28_file
)

log_message("Importing FOXL2 peaks")

FOXL2_peaks <- rtracklayer::import(
  FOXL2_file
)

log_message("Importing SOX9 peaks")

SOX9_peaks <- rtracklayer::import(
  SOX9_file
)

log_message("Importing DMRT1 peaks")

DMRT1_peaks <- rtracklayer::import(
  DMRT1_file
)

log_message(
  "Importing Sertoli-biased ATAC peaks"
)

ATAC_Sertoli_peaks <- rtracklayer::import(
  ATAC_Sertoli_file
)

log_message(
  "Importing Granulosa-biased ATAC peaks"
)

ATAC_Granulosa_peaks <- rtracklayer::import(
  ATAC_Granulosa_file
)

log_message(
  "Importing H3K9me3 ctrl peaks"
)

H3K9me3_ctrl_peaks <- rtracklayer::import(
  H3K9me3_ctrl_file
)

log_message(
  "Importing H3K9me3 cKO peaks"
)

H3K9me3_cKO_peaks <- rtracklayer::import(
  H3K9me3_cKO_file
)

log_message(
  "TRIM28 peaks: ",
  length(TRIM28_peaks)
)

log_message(
  "FOXL2 peaks: ",
  length(FOXL2_peaks)
)

log_message(
  "SOX9 peaks: ",
  length(SOX9_peaks)
)

log_message(
  "DMRT1 peaks: ",
  length(DMRT1_peaks)
)

log_message(
  "Sertoli-biased ATAC peaks: ",
  length(ATAC_Sertoli_peaks)
)

log_message(
  "Granulosa-biased ATAC peaks: ",
  length(ATAC_Granulosa_peaks)
)

log_message(
  "H3K9me3 ctrl peaks: ",
  length(H3K9me3_ctrl_peaks)
)

log_message(
  "H3K9me3 cKO peaks: ",
  length(H3K9me3_cKO_peaks)
)


###########################################
# Load gene-associated datasets
###########################################

log_message("Loading expressed gene list")

expressed_genes <- read.table(
  expressed_genes_file,
  header = FALSE,
  stringsAsFactors = FALSE
)[, 1]

expressed_genes <- unique(
  expressed_genes[
    !is.na(expressed_genes) &
      expressed_genes != ""
  ]
)

log_message(
  "Expressed genes loaded: ",
  length(expressed_genes)
)

log_message(
  "Loading 8-week DEG annotation"
)

DEG_8weeks <- load_DEG_direction(
  DEG_8weeks_file
)

log_message(
  "Loading 7-month DEG annotation"
)

DEG_7months <- load_DEG_direction(
  DEG_7months_file
)

log_message(
  "Loading sex-biased gene annotation"
)

sex_bias <- read.delim(
  sex_bias_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (
  !"gene" %in% colnames(sex_bias) ||
    !"sex.bias" %in% colnames(sex_bias)
) {
  log_error(
    "Sex-bias table must contain columns 'gene' and 'sex.bias'"
  )
}

sex_bias <- sex_bias[
  !is.na(sex_bias$gene) &
    sex_bias$gene != "",
  c(
    "gene",
    "sex.bias"
  ),
  drop = FALSE
]

sex_bias <- sex_bias[
  !duplicated(sex_bias$gene),
  ,
  drop = FALSE
]

sex_bias$sex.bias <- ifelse(
  tolower(sex_bias$sex.bias) == "sertoli",
  "Sertoli",
  ifelse(
    tolower(sex_bias$sex.bias) == "granulosa",
    "Granulosa",
    "None"
  )
)


###########################################
# Prepare genome annotation
###########################################

log_message(
  "Loading genome annotation"
)

genome_gtf <- rtracklayer::import(
  genome
)

log_message(
  "Preparing gene ID to symbol mapping"
)

gene2symbol <- GenomicRanges::mcols(
  genome_gtf
)[, c(
  "gene_id",
  "gene_name"
)]

gene2symbol <- unique(
  as.data.frame(
    gene2symbol
  )
)

gene2symbol <- gene2symbol[
  !is.na(gene2symbol$gene_id) &
    !is.na(gene2symbol$gene_name),
  ,
  drop = FALSE
]

gene2symbol <- gene2symbol[
  !duplicated(gene2symbol$gene_id),
  ,
  drop = FALSE
]

rownames(gene2symbol) <- gene2symbol$gene_id

log_message("Building TxDb object")

TxDb <- txdbmaker::makeTxDbFromGFF(
  genome
)


###########################################
# Build expressed gene annotation
###########################################

log_message(
  "Building expressed gene annotation"
)

gene_ranges <- GenomicFeatures::genes(
  TxDb
)

gene_ids <- names(
  gene_ranges
)

gene_symbols <- gene2symbol[
  gene_ids,
  "gene_name"
]

GenomicRanges::mcols(
  gene_ranges
)$gene_symbol <- gene_symbols

gene_ranges <- gene_ranges[
  !is.na(
    GenomicRanges::mcols(
      gene_ranges
    )$gene_symbol
  )
]

gene_ranges <- gene_ranges[
  GenomicRanges::mcols(
    gene_ranges
  )$gene_symbol %in%
    expressed_genes
]

gene_TSS <- GenomicRanges::promoters(
  gene_ranges,
  upstream = 0,
  downstream = 1
)

log_message(
  "Expressed genes matched to genome annotation: ",
  length(gene_ranges)
)


###########################################
# Compute overlaps
###########################################

log_message("Computing overlaps")

overlap_TRIM28 <- get_overlap_status(
  DAR_regions,
  TRIM28_peaks
)

overlap_FOXL2 <- get_overlap_status(
  DAR_regions,
  FOXL2_peaks
)

overlap_SOX9 <- get_overlap_status(
  DAR_regions,
  SOX9_peaks
)

overlap_DMRT1 <- get_overlap_status(
  DAR_regions,
  DMRT1_peaks
)

overlap_ATAC_Sertoli <- get_overlap_status(
  DAR_regions,
  ATAC_Sertoli_peaks
)

overlap_ATAC_Granulosa <- get_overlap_status(
  DAR_regions,
  ATAC_Granulosa_peaks
)

overlap_H3K9me3_ctrl <- get_overlap_status(
  DAR_regions,
  H3K9me3_ctrl_peaks
)

overlap_H3K9me3_cKO <- get_overlap_status(
  DAR_regions,
  H3K9me3_cKO_peaks
)

ATAC_bias <- make_ATAC_bias_category(
  overlap_sertoli = overlap_ATAC_Sertoli,
  overlap_granulosa = overlap_ATAC_Granulosa
)

log_message(
  "Regions overlapping TRIM28: ",
  sum(overlap_TRIM28)
)

log_message(
  "Regions overlapping FOXL2: ",
  sum(overlap_FOXL2)
)

log_message(
  "Regions overlapping SOX9: ",
  sum(overlap_SOX9)
)

log_message(
  "Regions overlapping DMRT1: ",
  sum(overlap_DMRT1)
)

log_message(
  "Regions overlapping H3K9me3 ctrl: ",
  sum(overlap_H3K9me3_ctrl)
)

log_message(
  "Regions overlapping H3K9me3 cKO: ",
  sum(overlap_H3K9me3_cKO)
)

log_message(
  "Regions with Sertoli-biased ATAC: ",
  sum(ATAC_bias == "Sertoli")
)

log_message(
  "Regions with Granulosa-biased ATAC: ",
  sum(ATAC_bias == "Granulosa")
)

log_message(
  "Regions with no sex-biased ATAC overlap: ",
  sum(ATAC_bias == "None")
)


###########################################
# Annotate genomic regions
###########################################

annotation_table <- annotate_regions(
  regions = DAR_regions,
  TxDb = TxDb,
  gene2symbol = gene2symbol,
  promoter = promoter
)

log_message(
  "Annotation rows: ",
  nrow(annotation_table)
)


###########################################
# Find nearest expressed gene
###########################################

nearest_expressed_gene <- get_nearest_expressed_gene(
  regions = DAR_regions,
  gene_ranges = gene_ranges,
  gene_TSS = gene_TSS,
  max_distance = nearest_gene_max_distance
)

log_message(
  "Nearest expressed gene vector length: ",
  length(nearest_expressed_gene)
)

log_message(
  "Regions with an expressed gene assigned: ",
  sum(!is.na(nearest_expressed_gene))
)


###########################################
# Annotate nearest expressed genes
###########################################

log_message(
  "Adding DEG and sex-bias annotations"
)

DEG_8weeks_status <- match_gene_annotation(
  genes = nearest_expressed_gene,
  annotation_table = DEG_8weeks,
  value_column = "direction",
  default = "None"
)

DEG_7months_status <- match_gene_annotation(
  genes = nearest_expressed_gene,
  annotation_table = DEG_7months,
  value_column = "direction",
  default = "None"
)

sex_bias_status <- match_gene_annotation(
  genes = nearest_expressed_gene,
  annotation_table = sex_bias,
  value_column = "sex.bias",
  default = "None"
)


###########################################
# Check output vector lengths
###########################################

log_message(
  "Checking output vector lengths"
)

region_id <- make_region_id(
  DAR_regions
)

expected_length <- length(
  DAR_regions
)

output_lengths <- c(
  Chromosome = length(
    GenomicRanges::seqnames(
      DAR_regions
    )
  ),
  Start = length(
    GenomicRanges::start(
      DAR_regions
    )
  ),
  End = length(
    GenomicRanges::end(
      DAR_regions
    )
  ),
  Region = length(
    region_id
  ),
  DAR_change = length(
    DAR$cluster
  ),
  TRIM28_overlap = length(
    overlap_TRIM28
  ),
  FOXL2_overlap = length(
    overlap_FOXL2
  ),
  SOX9_overlap = length(
    overlap_SOX9
  ),
  DMRT1_overlap = length(
    overlap_DMRT1
  ),
  H3K9me3_ctrl_overlap = length(
    overlap_H3K9me3_ctrl
  ),
  H3K9me3_cKO_overlap = length(
    overlap_H3K9me3_cKO
  ),
  ATAC_sex_bias = length(
    ATAC_bias
  ),
  Nearest_expressed_gene_50kb = length(
    nearest_expressed_gene
  ),
  Expression_change_8w = length(
    DEG_8weeks_status
  ),
  Expression_change_7m = length(
    DEG_7months_status
  ),
  Gene_sex_bias = length(
    sex_bias_status
  )
)

log_message(
  "Expected number of rows: ",
  expected_length
)

for (column_name in names(output_lengths)) {
  log_message(
    column_name,
    ": ",
    output_lengths[[column_name]]
  )
}

invalid_lengths <- output_lengths !=
  expected_length

if (any(invalid_lengths)) {
  log_error(
    "Invalid output vector length for: ",
    paste(
      names(output_lengths)[invalid_lengths],
      collapse = ", "
    )
  )
}


###########################################
# Build output table
###########################################

log_message(
  "Building final summary table"
)

output_table_df <- data.frame(
  Chromosome = as.character(
    GenomicRanges::seqnames(
      DAR_regions
    )
  ),
  Start = GenomicRanges::start(
    DAR_regions
  ),
  End = GenomicRanges::end(
    DAR_regions
  ),
  Region = region_id,
  `ATAC change` = DAR$cluster,
  `TRIM28 overlap` = ifelse(
    overlap_TRIM28,
    "Yes",
    "No"
  ),
  `FOXL2 overlap` = ifelse(
    overlap_FOXL2,
    "Yes",
    "No"
  ),
  `SOX9 overlap` = ifelse(
    overlap_SOX9,
    "Yes",
    "No"
  ),
  `DMRT1 overlap` = ifelse(
    overlap_DMRT1,
    "Yes",
    "No"
  ),
  `H3K9me3 Ctrl overlap` = ifelse(
    overlap_H3K9me3_ctrl,
    "Yes",
    "No"
  ),
  `H3K9me3 cKO overlap` = ifelse(
    overlap_H3K9me3_cKO,
    "Yes",
    "No"
  ),
  `ATAC sex bias` = ATAC_bias,
  `Nearest expressed gene 50 kb` = nearest_expressed_gene,
  `Expression change 8w` = DEG_8weeks_status,
  `Expression change 7m` = DEG_7months_status,
  `Gene sex bias` = sex_bias_status,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


###########################################
# Add ChIPseeker annotation
###########################################

log_message(
  "Merging ChIPseeker annotation"
)

colnames(annotation_table) <- c(
  "Region",
  "Genomic annotation",
  "Nearest gene",
  "TSS distance"
)

output_table_df <- merge(
  output_table_df,
  annotation_table,
  by = "Region",
  all.x = TRUE,
  sort = FALSE
)


###########################################
# Restore column order
###########################################

output_table_df <- output_table_df[, c(
  "Chromosome",
  "Start",
  "End",
  "Region",
  "ATAC change",
  "TRIM28 overlap",
  "FOXL2 overlap",
  "SOX9 overlap",
  "DMRT1 overlap",
  "H3K9me3 Ctrl overlap",
  "H3K9me3 cKO overlap",
  "ATAC sex bias",
  "Genomic annotation",
  "Nearest gene",
  "TSS distance",
  "Nearest expressed gene 50 kb",
  "Expression change 8w",
  "Expression change 7m",
  "Gene sex bias"
)]


###########################################
# Replace missing values
###########################################

output_table_df[[
  "Nearest expressed gene 50 kb"
]][
  is.na(
    output_table_df[[
      "Nearest expressed gene 50 kb"
    ]]
  )
] <- "None"

output_table_df[[
  "Expression change 8w"
]][
  is.na(
    output_table_df[[
      "Expression change 8w"
    ]]
  )
] <- "None"

output_table_df[[
  "Expression change 7m"
]][
  is.na(
    output_table_df[[
      "Expression change 7m"
    ]]
  )
] <- "None"

output_table_df[[
  "Gene sex bias"
]][
  is.na(
    output_table_df[[
      "Gene sex bias"
    ]]
  ) |
    output_table_df[[
      "Gene sex bias"
    ]] == ""
] <- "None"


###########################################
# Final checks
###########################################

log_message(
  "Final number of output rows: ",
  nrow(output_table_df)
)

log_message(
  "Final number of output columns: ",
  ncol(output_table_df)
)

log_message(
  "Regions without assigned expressed gene: ",
  sum(
    output_table_df[[
      "Nearest expressed gene 50 kb"
    ]] == "None"
  )
)

log_message(
  "8-week up-regulated nearest genes: ",
  sum(
    output_table_df[[
      "Expression change 8w"
    ]] == "Up"
  )
)

log_message(
  "8-week down-regulated nearest genes: ",
  sum(
    output_table_df[[
      "Expression change 8w"
    ]] == "Down"
  )
)

log_message(
  "7-month up-regulated nearest genes: ",
  sum(
    output_table_df[[
      "Expression change 7m"
    ]] == "Up"
  )
)

log_message(
  "7-month down-regulated nearest genes: ",
  sum(
    output_table_df[[
      "Expression change 7m"
    ]] == "Down"
  )
)


###########################################
# Save output table
###########################################

log_message(
  "Writing output table: ",
  output_table
)

write.table(
  output_table_df,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

log_message(
  "Analysis completed successfully"
)