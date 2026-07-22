source(".Rprofile")

suppressPackageStartupMessages({
	library("GenomicRanges")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# Load inputs
raw_peak_folder <- snakemake@input[["raw_peaks"]]
domains_folder <- snakemake@input[["domains"]]
genome <- snakemake@input[["genome"]]

# Load parameters
AB <- snakemake@params[["AB"]]
promoter <- snakemake@params[["promoter"]]
distance <- snakemake@params[["distance"]]
minimum_rep <- snakemake@params[["minimum_rep"]]
tmp_folder <- snakemake@params[["tmp"]]

# get outputs
consensus_peak_table <- snakemake@output[["consensus_peak_table"]]
consensus_domain_table <- snakemake@output[["consensus_domain_table"]]

consensus_peak_bed <- snakemake@output[["consensus_peak_bed"]]
consensus_domain_bed <- snakemake@output[["consensus_domain_bed"]]

gtf_files <- snakemake@output[["gtf_files"]]

# gtf_files <- snakemake@output[["gtf_files"]]

# Get the minimal distance between peaks that defines a domain
distance_between_peaks <- distance[AB]

# Prepare the genome for peak annotation
genome_gtf <- rtracklayer::import(genome)
gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id
TxDb <- txdbmaker::makeTxDbFromGFF(genome)

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

export_table <- function(gr, file) {
    region <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))
  cols <- c("annotation", "geneId", "distanceToTSS")
  df <- data.frame(
    chromosome = seqnames(gr),
    start = start(gr),
    end = end(gr),
    region = region,
    as.data.frame(mcols(gr[,cols])),
    stringsAsFactors = FALSE
  )
  write.table(
    df, file = file, sep = "\t",
    quote = FALSE, row.names = FALSE, col.names = TRUE
  )

}


get_consensus <- function(gr_list, min_overlap, distance){
  for (i in seq_along(gr_list)) {
    mcols(gr_list[[i]])$source_file <- paste0("file_", i)
  }

  all_gr <- do.call(c, gr_list)
  blocks <- disjoin(all_gr)

  presence_matrix <- sapply(gr_list, function(gr) {
    !is.na(findOverlaps(blocks, gr, select = "first"))
  })

  support_count <- rowSums(presence_matrix)
  kept_blocks <- blocks[support_count >= min_overlap]
  fused_regions <- reduce(kept_blocks, min.gapwidth = as.numeric(distance))

  peak_anno <- ChIPseeker::as.GRanges(
    ChIPseeker::annotatePeak(
      fused_regions,
      genomicAnnotationPriority = c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Intergenic"),
      tssRegion = c(-promoter, 0),
      TxDb = TxDb,
      level = "gene",
      overlap = "all"
    )
  )

  peak_anno$geneId <- gene2symbol[peak_anno$geneId, "gene_name"]

  return(peak_anno)
}

count_peaks <- function(file_path, peak_type){
	new_bed_files <- list.files(path = file_path, pattern = paste0(peak_type, "..bed"))
	peak_nb <- data.frame(sapply(paste0(file_path, "/", new_bed_files), R.utils::countLines))
	rownames(peak_nb) <- new_bed_files
	colnames(peak_nb) <- NULL
	write.csv(peak_nb, file=paste0(file_path, "/", peak_type, "_counts.csv"), sep = "\t", quote=FALSE, col.names=FALSE)
}

make_SAF <- function(gr){
	df <- as.data.frame(gr)
	genes <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))

	saf <- data.frame(
		genes,
		df
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

###################################################################################

# Peaks
bed_files <- list.files(path = raw_peak_folder, pattern = AB, full.names = TRUE)

gr_list <- lapply(bed_files, function(file){
  bed <- read.csv(file, header = FALSE, sep="\t")
	GR <- GRanges(
	  seqnames = bed[[1]],
	  ranges = IRanges(start = bed[[2]], end = bed[[3]]),
	  strand = "*"
	)

   return(GR)
})

print("Get consensus peaks")
consensus <- get_consensus(gr_list, minimum_rep, 1)
rtracklayer::export(
	consensus,
	consensus_peak_bed,
	format="BED"
)
export_table(
	consensus, 
	consensus_peak_table
)



# Domains
bed_files <- list.files(path = domains_folder, pattern = AB, full.names = TRUE)
print(bed_files)
print(domains_folder)


gr_list <- lapply(bed_files, rtracklayer::import)
print("get consensus domains cond")
consensus <- get_consensus(gr_list, minimum_rep, distance_between_peaks)
rtracklayer::export(
	consensus,
	consensus_domain_bed,
	format="BED"
)
export_table(
	consensus, 
	consensus_domain_table
)

make_SAF(consensus)


# count_peaks(consensus_cond, "peak")
# count_peaks(consensus_cond, "domain")


