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
AB <- snakemake@params[["AB"]]
bigwig_folder <- snakemake@params[["bigwig_folder"]]
new_bigwig_folder <- snakemake@params[["new_bigwig_folder"]]
output_file <- snakemake@output[["output_file"]]

###########################################
#                                         #
#           Normalize bigwig              #
#                                         #
###########################################

pattern <- "*.bigWig"

bigwig_files <- list.files(path = bigwig_folder, pattern = pattern)
print(bigwig_files)
print(AB)
bigwig_files <- grep(AB, bigwig_files, value = TRUE)
samples <- gsub(".bigWig", "", bigwig_files)
print(samples)

size_factors <- read.csv(size_factors, row.names = 1)

# norm_bigwig <- foreach(sample = samples) %dopar% {
norm_bigwig <- for (sample in samples) {
  print(sample)
  bw <- paste0(sample, ".bigWig")
  print(paste(bigwig_folder, bw, sep = "/"))
  bw_data <- rtracklayer::import.bw(paste(bigwig_folder, bw, sep = "/"))
  print("Import OK")
  sizeFactor <- size_factors[grep(tolower(sample), tolower(rownames(size_factors))), "SizeFactor"]
  norn_bw <- bw_data
  print(sample)
  print(sizeFactor)
  norn_bw$score <- bw_data$score * (1 / sizeFactor)
  new_bw <- paste0(sample, ".bw")
  print(paste("Writting new bw file to:", new_bigwig_folder, "/", new_bw))
  rtracklayer::export.bw(norn_bw, paste(new_bigwig_folder, new_bw, sep = "/"))
}

write.csv(size_factors, file = output_file)
