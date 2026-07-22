source(".Rprofile")

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

summary <- snakemake@input[["summary_counts"]]
merged_peaks <- snakemake@input[["merged_peaks"]]
peak_count <- snakemake@input[["peak_count"]]
domain_count <- snakemake@input[["domain_count"]]
output_pdf <- snakemake@output[["pdf_report"]]
output_png <- snakemake@output[["png_report"]]

###########################################
#                                         #
#             Load libraries              #
#                                         #
###########################################

suppressPackageStartupMessages({
	library("ggplot2")
	library("randomcoloR")
	library("scales")
	library("cowplot")
})

###########################################
#                                         #
#                Load data                #
#                                         #
###########################################

mapping_summary <- read.csv(summary)
samples <- mapping_summary$Sample
conditions <- gsub("_R.*", "\\1", samples)
group <- unique(conditions)

mapping_summary <- data.frame(
	group=conditions,
	mapping_summary
)

merged_peaks <- read.csv(merged_peaks, header = FALSE, row.names = 1)
rownames(merged_peaks) <- gsub(".seacr.peaks.stringent.bed", "", rownames(merged_peaks))
colnames(merged_peaks) <- "Domains"
merged_peaks$Sample <- rownames(merged_peaks)
merged_peaks$group <- gsub("_R.*", "\\1", rownames(merged_peaks))
merged_peaks <- merged_peaks[match(mapping_summary[-1, "Sample"], merged_peaks$Sample),]

consensus_peaks <- read.csv(peak_count, header = FALSE, row.names = 1)
rownames(consensus_peaks) <- gsub("_consensus_peaks.bed", "", rownames(consensus_peaks))
colnames(consensus_peaks) <- "Peaks"
consensus_peaks$group <- rownames(consensus_peaks)
consensus_peaks <- consensus_peaks[match(unique(mapping_summary[-1, "group"]), consensus_peaks$group),]

consensus_domains <- read.csv(domain_count, header = FALSE, row.names = 1)
rownames(consensus_domains) <- gsub("_consensus_domains.bed", "", rownames(consensus_domains))
colnames(consensus_domains) <- "Domains"
consensus_domains$group <- rownames(consensus_domains)
consensus_domains <- consensus_domains[match(unique(mapping_summary[-1, "group"]), consensus_domains$group),]
print(consensus_domains)

# Pick colors
set.seed(1234)
colours <- randomcoloR::distinctColorPalette(length(group))
names(colours) <- group

###########################################
#                                         #
#                Functions                #
#                                         #
###########################################

ratio <- 1
fontsize <- 12
text_size <- 4.2
legend_position <- "none"

plot_reads <- function(mapping_summary){
	print("Plot total reads")

	data <- mapping_summary
	data$label <- paste0(round(data$Total/1000000,1), "M")
	data$Total <- data$Total/1000000

	reads <- ggplot(data, aes(Total, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Million fragments counts") +
		theme_light() +
		ggtitle("Read pairs per samples\n(fragments)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(reads)
}

plot_trimmed <- function(mapping_summary){
	print("Plot trimmed reads")
	data <- mapping_summary
	data$label <- paste0(round((data$Trimmed*100)/data$Total,2), "%")
	data$Trimmed <- (data$Trimmed*100)/data$Total

	mapping <- ggplot(data, aes(Trimmed, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of kept fragments after trimming\n(rel. to total frag.)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(mapping)
}


plot_mapped <- function(mapping_summary){
	print("Plot mapped reads")
	data <- mapping_summary
	data$label <- paste0(round((data$Mapped*100)/data$Total,2), "%")
	data$Mapped <- (data$Mapped*100)/data$Total

	mapping <- ggplot(data, aes(Mapped, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of mapped fragment\n(rel. to trimmed frag.)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(mapping)
}

plot_multimapped <- function(mapping_summary){
	print("Plot multimapped reads")
	data <- mapping_summary
	data$label <- paste0(round((data$Multimapped*100)/data$Total,2), "%")
	data$Multimapped <- (data$Multimapped*100)/data$Total

	mapping <- ggplot(data, aes(Multimapped, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#444444", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of Multimapped fragments\n(rel. to mapped frag.)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(mapping)
}

plot_dups <- function(mapping_summary){
	print("Plot duplicated reads")
	data <- mapping_summary
	data$label <- paste0(round((data$Duplicated*100)/data$Total,2), "%")
	data$Duplicated <- (data$Duplicated*100)/data$Total

	mapping <- ggplot(data, aes(Duplicated, factor(Sample, level = mapping_summary$Sample), fill=group)) +
		geom_bar(stat="identity") +
		geom_text(aes(label=label), hjust=-0.15, color="#444444", size=text_size) +
		scale_x_continuous(labels = comma, limits=c(0,100)) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Percentages") +
		theme_light() +
		ggtitle("% of Duplicated fragments\n(rel. to mapped frag.)") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)

	return(mapping)
}

plot_peak_calling <- function(mapping_summary){
	print("Plot peaks")
	data <- mapping_summary
	data$label <- scales::comma(data$Peaks)

	data <- data[-1,]

	colours <- colours[-1]
	group <- group[-1]

	peak_calling <- ggplot(data, aes(Peaks, factor(Sample, level = mapping_summary$Sample))) +
		geom_bar(stat="identity", aes(fill=group)) +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Peak counts") +
		theme_light() +
		ggtitle("Number of called peaks") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)
	return(peak_calling)
}

plot_consensus_peaks <- function(consensus_peaks){
	print("Plot consensus peaks")

	data <- consensus_peaks
	data$label <- scales::comma(data$Peaks)

	colours <- colours[-1]
	group <- group[-1]

	peak_calling <- ggplot(data, aes(Peaks, factor(group, level = consensus_peaks$group))) +
		geom_bar(stat="identity", aes(fill=group)) +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Peak counts") +
		theme_light() +
		ggtitle("Number of consensus peaks") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)
	return(peak_calling)
}

plot_merged_peaks <- function(merged_peaks){
	print("Plot domains")
	data <- merged_peaks
	data$label <- scales::comma(data$Domains)
	colours <- colours[-1]
	group <- group[-1]

	peak_calling <- ggplot(data, aes(Domains, factor(Sample, level = merged_peaks$Sample))) +
		geom_bar(stat="identity", aes(fill=group)) +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Domain counts") +
		theme_light() +
		ggtitle("Number of enriched domains") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)
	return(peak_calling)
}


plot_consensus_domains <- function(consensus_domains){
	print("Plot consensus domains")

	data <- consensus_domains

	data$label <- scales::comma(data$Domains)

	colours <- colours[-1]
	group <- group[-1]

	peak_calling <- ggplot(data, aes(Domains, factor(group, level = consensus_domains$group))) +
		geom_bar(stat="identity", aes(fill=group)) +
		geom_text(aes(label=label), hjust=1.1, color="white", size=text_size) +
		scale_x_continuous(labels = comma) +
		scale_fill_manual(values=colours) +
		scale_y_discrete(limits=rev) +
		ylab("Samples") + 
		xlab("Domain counts") +
		theme_light() +
		ggtitle("Number of consensus domains") +
		theme(
			plot.title = element_text(size=fontsize, hjust = 0.5, face="bold"),
			axis.text=element_text(size=fontsize),
			axis.title=element_text(size=fontsize),
			legend.title=element_blank(),
			legend.text=element_text(size=fontsize),
			legend.position=legend_position,
			aspect.ratio=ratio
		)
	return(peak_calling)
}


reads <- plot_reads(mapping_summary)
trimmed <- plot_trimmed(mapping_summary)
mapping <- plot_mapped(mapping_summary)
multimapped <- plot_multimapped(mapping_summary)
dups <- plot_dups(mapping_summary)
peak <- plot_peak_calling(mapping_summary)
consensus_peaks <- plot_consensus_peaks(consensus_peaks)
domains <- plot_merged_peaks(merged_peaks)
consensus_domains <- plot_consensus_domains(consensus_domains)

figure <- plot_grid(
	reads,
	trimmed,
	mapping,
	multimapped,
	dups,
	peak,
	consensus_peaks,
	domains,
	consensus_domains,
	ncol=3, 
	align="hv", 
	axis="bt", 
	labels="AUTO",
	label_size = 16
)

save_plot(
	output_pdf, 
	figure, 
	dpi=300, 
	base_width=40, 
	base_height=32,
	units = c("cm")
)

save_plot(
	output_png, 
	figure, 
	dpi=300, 
	bg = "white", 
	base_width=40, 
	base_height=32,
	units = c("cm")
)
