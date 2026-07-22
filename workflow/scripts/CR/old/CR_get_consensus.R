source(".Rprofile")

suppressPackageStartupMessages({
	library("GenomicRanges")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_peak_folder <- snakemake@input[["raw_peaks"]]
domains_folder <- snakemake@input[["domains"]]

promoter <- snakemake@params[["promoter"]]
distance <- snakemake@params[["distance"]]
minimum_rep <- snakemake@params[["minimum_rep"]]
tmp_folder <- snakemake@params[["tmp"]]

consensus_cond <- snakemake@output[["consensus_cond"]]
consensus_AB <- snakemake@output[["consensus_AB"]]

dir.create(consensus_cond)
dir.create(consensus_AB)

# Get condition from file names
raw_bed_files <- list.files(path = raw_peak_folder, pattern = ".bed")
conditions <- unique(gsub("_R..seacr.peaks.stringent.bed", "", raw_bed_files))

# Get AB from file names
ABs <- unique(sapply(strsplit(conditions, "_"), `[`, 2))

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
	print(as.numeric(distance))
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

  print(peak_anno)

  return(peak_anno)
}

count_peaks <- function(file_path, peak_type){
	new_bed_files <- list.files(path = file_path, pattern = paste0(peak_type, "..bed"))
	peak_nb <- data.frame(sapply(paste0(file_path, "/", new_bed_files), R.utils::countLines))
	rownames(peak_nb) <- new_bed_files
	colnames(peak_nb) <- NULL
	write.csv(peak_nb, file=paste0(file_path, "/", peak_type, "_counts.csv"), sep = "\t", quote=FALSE, col.names=FALSE)
}

make_SAF <- function(gr, AB){
	df <- as.data.frame(gr)
	genes <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))

	saf <- data.frame(
		genes,
		df
	)

	write.table(
		saf, 
		file = paste0(tmp_folder, "/", AB, "_consensus_peaks.gtf"),
		quote = FALSE,
		row.names = FALSE,
		col.names = FALSE,
		sep = "\t"
	)
}

###################################################################################

gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
alt_gtf_url <- "http://ftp.cbi.pku.edu.cn/pub/mirror/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

genome_gtf <- tryCatch({
  rtracklayer::import(gtf_url)
 }, error = function(e) {
  rtracklayer::import(alt_gtf_url)
 })

gene2symbol <- GenomicRanges::mcols(genome_gtf)[, c("gene_id", "gene_name")]
gene2symbol <- unique(gene2symbol)
rownames(gene2symbol) <- gene2symbol$gene_id

TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)

# Consensus peaks and domains per condition
for (condition in conditions){

	# Peaks
	bed_files <- list.files(path = raw_peak_folder, pattern = condition, full.names = TRUE)
  gr_list <- lapply(bed_files, function(file){
    bed <- read.csv(file, header = FALSE, sep="\t")
		GR <- GRanges(
		  seqnames = bed[[1]],
		  ranges = IRanges(start = bed[[2]], end = bed[[3]]),
		  strand = "*"
		)

    return(GR)
    }
  )

	print("get consensus peaks cond")
	consensus <- get_consensus(gr_list, minimum_rep, 1)
	rtracklayer::export(
		consensus,
		paste0(consensus_cond, "/", condition, "_consensus_peaks.bed"),
		format="BED"
	)
	export_table(
		consensus, 
		paste0(consensus_cond, "/", condition, "_consensus_peaks.txt")
	)

	AB <- sapply(strsplit(condition, "_"), `[`, 2)
	dist <- distance[AB]
	# Domains
	bed_files <- list.files(path = domains_folder, pattern = condition, full.names = TRUE)
  gr_list <- lapply(bed_files, rtracklayer::import)
	print("get consensus domains cond")
	consensus <- get_consensus(gr_list, minimum_rep, dist)
	rtracklayer::export(
		consensus,
		paste0(consensus_cond, "/", condition, "_consensus_domains.bed"),
		format="BED"
	)
	export_table(
		consensus, 
		paste0(consensus_cond, "/", condition, "_consensus_domains.txt")
	)
}



# Consensus peaks and domains per AB
for (AB in ABs){
	dist <- distance[AB]
	# Peaks
	bed_files <- list.files(path = consensus_cond, pattern = paste0(AB, "_consensus_peaks.bed"), full.names = TRUE)
  gr_list <- lapply(bed_files, rtracklayer::import)
	print("get consensus peaks AB")
	consensus <- get_consensus(gr_list, minimum_rep, 1)
	rtracklayer::export(
		consensus,
		paste0(consensus_AB, "/", AB, "_consensus_peaks.bed"),
		format="BED"
	)
	export_table(
		consensus, 
		paste0(consensus_AB, "/", AB, "_consensus_peaks.txt")
	)

	make_SAF(consensus, AB)

	# Domains
	bed_files <- list.files(path = consensus_cond, pattern = paste0(AB, "_consensus_domains.bed"), full.names = TRUE)
  # gr_list <- lapply(bed_files, function(file){
  #   bed <- read.csv(file, header = FALSE, sep="\t")
  #   GR <- GRanges(bed[,6])
  #   return(GR)
  #   }
  # )
  
  gr_list <- lapply(bed_files, rtracklayer::import)
  print("get consensus domains AB")
	consensus <- get_consensus(gr_list, minimum_rep, dist)
	rtracklayer::export(
		consensus,
		paste0(consensus_AB, "/", AB, "_consensus_domains.bed"),
		format="BED"
	)
	export_table(
		consensus, 
		paste0(consensus_AB, "/", AB, "_consensus_domains.txt")
	)
}

count_peaks(consensus_cond, "peak")
count_peaks(consensus_cond, "domain")
count_peaks(consensus_AB, "peak")
count_peaks(consensus_AB, "domain")

