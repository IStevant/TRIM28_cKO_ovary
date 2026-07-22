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
})


# FANTOM5 enhancers (mm10)
fantom_url <- "https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/mouse_permissive_enhancers_phase_1_and_2.bed.gz"
enh_fantom <- import(fantom_url)

# ENCODE mouse cCREs (enhancer-like regions)
encode_url <- "https://www.encodeproject.org/files/ENCFF167FJQ/@@download/ENCFF167FJQ.bed.gz"

# Read the BED file manually
encode_df <- read_tsv(encode_url,
                      col_names = FALSE,
                      comment = "#",
                      col_types = cols(
                        X1 = col_character(),  # chrom
                        X2 = col_integer(),    # start
                        X3 = col_integer(),    # end
                        X4 = col_character(),  # name
                        X5 = col_integer(),    # score
                        X6 = col_character(),  # strand (often ".")
                        X7 = col_integer(),    # thickStart
                        X8 = col_integer(),    # thickEnd
                        X9 = col_character(),  # itemRgb
                        X10 = col_character()  # label (e.g. "CA", "CA-CTCF")
                      ))

# Optional: filter to enhancers only (e.g. CA elements)
encode_enhancers <- encode_df[grepl("^CA", encode_df$X10), ]

# Convert to GRanges
enh_encode <- GRanges(
  seqnames = encode_enhancers$X1,
  ranges = IRanges(start = encode_enhancers$X2 + 1, end = encode_enhancers$X3),
  name = encode_enhancers$X4,
  label = encode_enhancers$X10
)

# Combine and reduce overlapping enhancers
# Keep only seqnames, start, end — discard metadata
enh_fantom_clean <- GRanges(seqnames = seqnames(enh_fantom),
                            ranges = ranges(enh_fantom))

enh_encode_clean <- GRanges(seqnames = seqnames(enh_encode),
                            ranges = ranges(enh_encode))

# Now safely combine and reduce
enhancers_all <- reduce(c(enh_fantom_clean, enh_encode_clean))

df <- data.frame(
  seqnames=seqnames(enhancers_all),
  start=format(start(enhancers_all)-1, scientific=F),
  end=format(end(enhancers_all), scientific=F),
  name=c(rep(".", length(enhancers_all))),
  score=c(rep(".", length(enhancers_all))),
  strand=strand(enhancers_all))


write.table(df, file=snakemake@output[["enhancers"]], quote=F, sep="\t", row.names=F, col.names=F)
