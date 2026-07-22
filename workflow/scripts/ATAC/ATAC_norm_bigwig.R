source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
suppressPackageStartupMessages({
  library("dplyr")
  library("doParallel")
  library("foreach")
  library("rtracklayer")
})

doParallel::registerDoParallel(cores = 12)

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

size_factors <- snakemake@input[["size_factors"]]
bigwig_folder <- snakemake@params[["bigwig_folder"]]
new_bigwig_folder <- snakemake@params[["new_bigwig_folder"]]
output_file <- snakemake@output[["output_file"]]

###########################################
#                                         #
#           Normalize bigwig              #
#                                         #
###########################################

pattern <- "*.mLb.clN.bigWig"

bigwig_files <- list.files(path = bigwig_folder, pattern = pattern)
samples <- gsub(".mLb.clN.bigWig", "", bigwig_files)

size_factors <- read.csv(size_factors, row.names = 1)

# norm_bigwig <- foreach(sample = samples) %dopar% {
norm_bigwig <- for (sample in samples) {
  print(sample)
  bw <- paste0(sample, ".mLb.clN.bigWig")
  print(paste(bigwig_folder, bw, sep = "/"))
  bw_data <- rtracklayer::import.bw(paste(bigwig_folder, bw, sep = "/"))
  print("Import OK")
  scale_folder <- paste(bigwig_folder, "scale", sep = "/")
  scale_file <- list.files(path = scale_folder, pattern = paste0(sample, ".*"))
  scale <- read.csv(paste(scale_folder, scale_file, sep = "/"), header = FALSE)[1, 1]
  print(paste("Original scaling factor:",scale))
  denorm_bw <- bw_data
  denorm_bw$score <- bw_data$score / scale
  sizeFactor <- size_factors[grep(tolower(sample), tolower(rownames(size_factors))), "SizeFactor"]
  norn_bw <- denorm_bw
  norn_bw$score <- denorm_bw$score * (1 / sizeFactor)
  new_bw <- paste0(sample, ".bw")
  print(paste("Writting new bw file to:", new_bigwig_folder, "/", new_bw))
  rtracklayer::export.bw(norn_bw, paste(new_bigwig_folder, new_bw, sep = "/"))
}

write.csv(size_factors, file = output_file)
