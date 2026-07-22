source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

raw_counts <- read.csv(file = snakemake@input[["counts"]], row.names = 1)
samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)

# Promoter region, i.e. distance to TSS
promoter <- snakemake@params[["promoter"]]

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

adj.pval <- snakemake@params[["adjpval"]]
log2FC <- snakemake@params[["log2FC"]]

save_folder <- snakemake@params[["save_folder"]]

###########################################
#                                         #
#    DESeq2 analysis between conditions   #
#                                         #
###########################################

DESeqObj <- DESeq2::DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = samplesheet,
  design = ~conditions
)

# Get conditions
conditions <- unique(samplesheet$conditions)

if (length(conditions)==2) {
  DERs <- DESeq2::DESeq(DESeqObj)
} else if (length(conditions)>2) {
  DERs <- DESeq2::DESeq(DESeqObj, test = "LRT", reduced = ~1)
} else {
  stop("Not enough conditions to perform differential analysis.")
}

resDERs <- DESeq2::results(DERs)
resDERs$symbol <- GenomicRanges::mcols(resDERs)$symbol

sig.DE <- subset(resDERs, padj < adj.pval)
sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

DER_GR <- GenomicRanges::GRanges(rownames(sig.DE))

TxDb <- txdbmaker::makeTxDbFromGFF(gtf_url)

DE_anno <- as.data.frame(
  ChIPseeker::annotatePeak(
    DER_GR,
    genomicAnnotationPriority = c("Promoter", "5UTR", "Exon", "Intron", "3UTR", "Downstream", "Intergenic"),
    tssRegion = c(-promoter, 0),
    TxDb = TxDb,
    level = "gene",
    overlap = "all"
  )
)

DE_anno$geneId <- gene2symbol[DE_anno$geneId, "gene_name"]

filtered_DERs <- data.frame(
  sig.DE[, -c(3:4)],
  annotation = DE_anno$annotation,
  nearest.gene = DE_anno$geneId,
  distanceToTSS = DE_anno$distanceToTSS
)

unfiltered_DERS <- data.frame(resDERs)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################
write.table(filtered_DERs[order(filtered_DERs$padj), ], file = snakemake@output[["tsv"]], quote = FALSE, sep="\t")

write.table(unfiltered_DERS, file = paste0(snakemake@output[["tsv"]], "_no_filter"), quote = FALSE, sep="\t")

filtered_DERs <- rownames(filtered_DERs)

save(filtered_DERs, file = snakemake@output[["sig_DERs"]])

rtracklayer::export(GenomicRanges::GRanges(filtered_DERs), snakemake@output[["DER_bed"]], "bed")

