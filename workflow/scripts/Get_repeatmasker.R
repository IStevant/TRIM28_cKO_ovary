source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("rtracklayer")
  library("GenomicRanges")
  library("GenomicFeatures")
  library("readr")
  library("R.utils")
  library("data.table")
})


# Download RepeatMasker table for mm10 (UCSC)
rmsk_url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz"

# Load and parse RepeatMasker table
rmsk_table <- fread(cmd = paste("curl -s", rmsk_url, "| gunzip -c"), sep = "\t", header = FALSE, fill = TRUE)
colnames(rmsk_table)[c(6:8,10:11)] <- c("chrom", "start", "end", "strand", "repName")
colnames(rmsk_table)[12:13] <- c("repClass", "repFamily")
colnames(rmsk_table)[3] <- c("Kimura")

# Convert to GRanges object
rmsk_gr <- GRanges(
  seqnames = rmsk_table$chrom,
  ranges = IRanges(start = rmsk_table$start + 1, end = rmsk_table$end),
  strand = rmsk_table$strand,
  repName = rmsk_table$repName,
  repClass = rmsk_table$repClass,
  repFamily = rmsk_table$repFamily,
  Kimura_pct = rmsk_table$Kimura/10
)

df <- data.frame(
  seqnames=seqnames(rmsk_gr),
  start=format(start(rmsk_gr)-1, scientific=F),
  end=format(end(rmsk_gr), scientific=F),
  strand=strand(rmsk_gr),
  repClass = rmsk_table$repClass,
  repFamily = rmsk_table$repFamily,
  repName = rmsk_table$repName,
  Kimura_pct = rmsk_table$Kimura/10)
  

# List repeat class to analyse
keptRepeats <- c(
  "DNA",
  "LINE",
  "LTR",
  "SINE"
)

# Filter repeats
# df <- df[df$repClass %in% keptRepeats,]

write.table(df, file=snakemake@output[["repeats"]], quote=F, sep="\t", row.names=F)
