source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

peak_folder <- snakemake@params[["raw_peaks"]]
distance <- snakemake@params[["distance"]]
output_folder <- snakemake@output[["output_folder"]]
output_file <- snakemake@output[["peak_counts"]]

###########################################
#                                         #
#              Merge peaks                #
#                                         #
###########################################

# For each defined distance (for each antibody, a minimal distance is set in the "analysis_parameters.yaml" file)
## Get the antibody name
## Get the associated bed files
## Get the associated minimal distance
## Merge the peaks when the distance between the peaks are inferior to the minimal distance and generate the new bed file
bedfiles <- lapply(seq_along(distance), function(AB) {
	AB_name <- names(distance[AB])
	bed_files <- list.files(path = peak_folder, pattern = AB_name)
	dist <- distance[AB]
	for (bed in bed_files){
		bedtoolsr::bt.merge(
			i = paste0(peak_folder, "/", bed),
			d = dist,
			output = paste0(output_folder, "/", bed)
		)
	}
})

# Get the newly created beds
new_bed_files <- list.files(path = output_folder, pattern = "*.bed")

# Create a dataframe and count the number of domains per file
peak_nb <- data.frame(sapply(paste0(output_folder, "/", new_bed_files), R.utils::countLines))
rownames(peak_nb) <- new_bed_files
colnames(peak_nb) <- NULL

# Save the count file
write.csv(peak_nb, file=output_file, sep = "\t", quote=FALSE, col.names=FALSE)
