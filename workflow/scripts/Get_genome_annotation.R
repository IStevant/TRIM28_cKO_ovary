source(".Rprofile")


###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("R.utils")
})


# Alternative URLs for the mouse genome
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
alt_gtf_url <- "http://ftp.cbi.pku.edu.cn/pub/mirror/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"


# Download the genome from internet, if the server is down, use the mirror website instead
tryCatch({
  download.file(gtf_url, destfile = snakemake@output[["genome"]])
}, error = function(e) {
  download.file(alt_gtf_url, destfile = snakemake@output[["genome"]])
})


