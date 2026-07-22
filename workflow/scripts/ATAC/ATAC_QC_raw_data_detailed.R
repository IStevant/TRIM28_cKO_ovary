source(".Rprofile")

###########################################
# Load libraries
###########################################

suppressPackageStartupMessages({
  library("ChIPpeakAnno")
  library("cowplot")
  library("data.table")
  library("eulerr")
  library("futile.logger")
  library("GenomicRanges")
  library("ggplot2")
  library("gridExtra")
  library("IRanges")
  library("randomcoloR")
  library("rtracklayer")
  library("scales")
  library("shades")
  library("stringr")
  library("UpSetR")
  library("gVenn")
  library("ComplexHeatmap")
})


###########################################
# Logging functions
###########################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}


###########################################
# Snakemake inputs, outputs and parameters
###########################################

log_message("Reading Snakemake inputs, outputs and parameters.")

summary_file <- snakemake@input[["summary_counts"]]

mapping_dir <- snakemake@params[["mapping_folder"]]
peak_type <- snakemake@params[["peak_type"]]
min_peak_rep <- snakemake@params[["min_peak_rep"]]
raw_peak_folder <- snakemake@params[["raw_peak_folder"]]

output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]


###########################################
# Global plotting parameters
###########################################

ratio <- 0.4
fontsize <- 12
text_size <- 3.8
legend_position <- "none"


###########################################
# Helper functions
###########################################

find_output_file <- function(files, group, extension) {
  expected_name <- paste0("ATAC_raw_data_detailed_QC_", group, ".", extension)
  selected_file <- files[basename(files) == expected_name]

  if (length(selected_file) != 1) {
    log_error(
      "Could not uniquely identify output file for group: ",
      group,
      " and extension: ",
      extension,
      "\nExpected file name: ",
      expected_name,
      "\nAvailable files:\n",
      paste(files, collapse = "\n")
    )
  }

  return(selected_file)
}


get_sample_label <- function(sample) {
  stringr::str_to_sentence(sapply(strsplit(sample, "_"), `[`, 2))
}


get_group_from_sample <- function(sample) {
  gsub("_REP.*", "\\1", sample)
}


###########################################
# Load data
###########################################

log_message("Loading mapping summary table.")

mapping_summary <- read.csv(summary_file)

samples <- mapping_summary$Sample
groups <- get_group_from_sample(samples)
group_list <- unique(groups)

mapping_summary <- data.frame(
  group = groups,
  mapping_summary
)

log_message("Detected groups: ", paste(group_list, collapse = ", "))

set.seed(1234)
colours <- randomcoloR::distinctColorPalette(length(group_list))
names(colours) <- group_list


###########################################
# Mapping summary plots
###########################################

plot_reads <- function(mapping_summary, group) {
  log_message("Plotting total reads for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)
  data$label <- scales::comma(data$Total)
  data$Total <- data$Total / 1000000

  ggplot(data, aes(Total, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Million read counts") +
    ggtitle("Read counts per samples (R1+R2)") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


plot_mapped <- function(mapping_summary, group) {
  log_message("Plotting mapped reads for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)
  data$label <- paste0(round((data$Mapped * 100) / data$Total, 2), "%")
  data$Mapped <- (data$Mapped * 100) / data$Total

  ggplot(data, aes(Mapped, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma, limits = c(0, 100)) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle("% of mapped reads per samples\n(rel. to total reads)") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


plot_duplicated <- function(mapping_summary, group) {
  log_message("Plotting duplicated reads for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)
  data$dup <- (data$Duplicated * 100) / data$Mapped
  data$label <- paste0(round(data$dup, 2), "%")

  ggplot(data, aes(dup, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = -0.15, color = "#444444", size = text_size) +
    scale_x_continuous(labels = comma, limits = c(0, 100)) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle("% of duplicated reads\n(rel. to mapped reads)") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


plot_MT <- function(mapping_summary, group) {
  log_message("Plotting mitochondrial reads for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)
  data$MT <- (data$Mitochondrial * 100) / data$Mapped
  data$label <- paste0(round(data$MT, 2), "%")

  ggplot(data, aes(MT, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = -0.15, color = "#444444", size = text_size) +
    scale_x_continuous(labels = comma, limits = c(0, 100)) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle("% of mitochondrial reads\n(rel. to mapped reads)") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


plot_peak_calling <- function(mapping_summary, group) {
  log_message("Plotting peak counts for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)
  data$label <- scales::comma(data$Peaks)

  ggplot(data, aes(Peaks, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 4)) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Peak counts") +
    ggtitle("Number of peaks per sample") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


plot_frip <- function(mapping_summary, group) {
  log_message("Plotting FRiP scores for ", group)

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- get_sample_label(data$Sample)

  ggplot(data, aes(FRiP, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = round(FRiP, 2)), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("FRiP score") +
    ggtitle("Fraction of reads in peaks per sample") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position,
      aspect.ratio = ratio
    )
}


###########################################
# Insert size plots
###########################################

plot_insert_size_single <- function(path, sample) {
  table <- read.table(
    file = file.path(path, sample),
    skip = 10,
    header = TRUE
  )

  title <- gsub(".mLb.*", "\\1", sample)
  title <- stringr::str_to_sentence(paste("Insert size of", sapply(strsplit(title, "_"), `[`, 2)))

  group <- get_group_from_sample(sample)
  colour <- colours[group]

  ggplot(table, aes(x = insert_size, y = All_Reads.fr_count)) +
    geom_area(fill = colour) +
    scale_y_continuous(labels = comma) +
    geom_vline(
      xintercept = c(147, 294, 441),
      linetype = "dashed",
      colour = "#666666"
    ) +
    annotate("text", x = 0, y = max(table$All_Reads.fr_count), hjust = 1.1, label = "Nucleosome free", color = "#666666", size = 3.5, angle = 90) +
    annotate("text", x = 167, y = max(table$All_Reads.fr_count), hjust = 1.1, label = "Mononucleosome", color = "#666666", size = 3.5, angle = 90) +
    annotate("text", x = 314, y = max(table$All_Reads.fr_count), hjust = 1.1, label = "Dinucleosome", color = "#666666", size = 3.5, angle = 90) +
    ylab("Coverage") +
    xlab("Insert size (nt)") +
    xlim(0, 800) +
    ggtitle(title) +
    theme_light() +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 11),
      legend.title = element_blank()
    )
}


plot_insert_size <- function(group) {
  log_message("Plotting insert size distributions for ", group)

  path <- file.path(mapping_dir, "bwa", "merged_library", "picard_metrics")
  file_names <- list.files(path = path, pattern = paste0(group, ".*insert_size_metrics"))

  if (length(file_names) == 0) {
    log_error("No insert size metric files found for group: ", group)
  }

  lapply(file_names, function(file) {
    plot_insert_size_single(path, file)
  })
}


###########################################
# Peak annotation plot
###########################################

plot_anno <- function(group) {
  log_message("Plotting peak annotation for ", group)

  annotation_file <- file.path(
    mapping_dir,
    "multiqc",
    paste0(peak_type, "_peak"),
    "multiqc_data",
    "multiqc_mlib_peak_annotation-plot.txt"
  )

  anno <- read.csv(annotation_file, sep = "\t")
  anno_group <- anno[grep(group, anno$Sample), ]

  sample <- stringr::str_to_sentence(sub(".*_(REP[0-9]+)", "\\1", anno_group$Sample))

  group_colour <- colours[group]
  anno_col <- group_colour %>%
    brightness(0.85) %>%
    saturation(seq(from = 0.2, to = 1, length.out = 4))

  data <- data.frame(
    sample = as.factor(sample),
    count = c(
      anno_group$Intergenic,
      anno_group$promoter.TSS,
      anno_group$exon,
      anno_group$intron
    ),
    annotation = c(
      rep("Intergenic", nrow(anno_group)),
      rep("Promoter", nrow(anno_group)),
      rep("Exon", nrow(anno_group)),
      rep("Intron", nrow(anno_group))
    )
  )

  data$proportion <- ave(data$count, data$sample, FUN = function(x) x / sum(x))
  data$label <- scales::percent(data$proportion, accuracy = 1)
  data$label[data$proportion < 0.03] <- ""

  ggplot(data, aes(x = sample, y = count, fill = annotation)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(
      aes(label = label),
      position = position_fill(0.5),
      size = 3.5
    ) +
    scale_fill_manual(values = alpha(anno_col, 0.6)) +
    xlab("Samples") +
    ylab("Proportions") +
    ggtitle("Peak annotation") +
    theme_light() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      aspect.ratio = 0.8
    )
}


###########################################
# Replicate overlap plots
###########################################

get_overlap <- function(group) {
  log_message("Computing peak overlap for ", group)

  files <- list.files(
    path = raw_peak_folder,
    pattern = paste0(group, ".*", peak_type, "Peak"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    log_error("No peak files found for group: ", group)
  }

  peaks <- lapply(files, rtracklayer::import)
  names(peaks) <- stringr::str_to_sentence(stringr::str_extract(basename(files), "REP[0-9]+"))

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  makeVennDiagram(
    Peaks = peaks,
    NameOfPeaks = names(peaks),
    plot = FALSE
  )
}


plot_venn <- function(group, overlap) {

  log_message("Plotting replicate overlap UpSet plot for ", group)

  files <- list.files(
    path = raw_peak_folder,
    pattern = paste0(group, ".*", peak_type, "Peak"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    log_error("No peak files found for group: ", group)
  }

  peaks <- lapply(
    files,
    rtracklayer::import
  )

  names(peaks) <- stringr::str_to_sentence(
    stringr::str_extract(
      basename(files),
      "REP[0-9]+"
    )
  )

  overlaps <- gVenn::computeOverlaps(
    peaks
  )

  group_colour <- colours[group]

  upset_grob <- grid::grid.grabExpr({
    ComplexHeatmap::draw(
      gVenn::plotUpSet(
        overlaps,
        comb_col = group_colour
      )
    )
  })

  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Peak overlap between replicates",
      fontface = "bold",
      size = 12
    )

  cowplot::plot_grid(
    title,
    upset_grob,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
}

get_summary_table <- function(overlap) {
  log_message("Preparing peak overlap summary table.")

  venn_counts <- as.data.frame(overlap$vennCounts)
  venn_counts <- venn_counts[complete.cases(venn_counts), ]

  presence_matrix <- venn_counts[, -ncol(venn_counts)]
  row_sums <- rowSums(presence_matrix)

  peak_counts <- t(data.frame(
    total = sum(venn_counts[, "Counts"]),
    common_all = sum(venn_counts[row_sums >= ncol(presence_matrix), "Counts"]),
    common_min_rep = sum(venn_counts[row_sums >= min_peak_rep, "Counts"])
  ))

  peak_counts <- data.frame(
    Counts = scales::comma(peak_counts),
    Percent = paste0(round((100 * peak_counts) / peak_counts[1], digits = 2), "%")
  )

  rownames(peak_counts) <- c(
    "Total",
    "Peaks common to all rep.",
    paste("Peaks found in at least", min_peak_rep, "rep.")
  )

  theme_table <- gridExtra::ttheme_default(
    core = list(fg_params = list(cex = 0.9)),
    colhead = list(fg_params = list(cex = 0.9)),
    rowhead = list(fg_params = list(cex = 0.9))
  )

  gridExtra::tableGrob(peak_counts, theme = theme_table)
}


get_peak_size_dist <- function(group) {
  log_message("Plotting peak size distribution for ", group)

  files <- list.files(
    path = raw_peak_folder,
    pattern = paste0(group, ".*", peak_type, "Peak"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    log_error("No peak files found for group: ", group)
  }

  peaks <- lapply(files, rtracklayer::import)
  all_peaks <- unlist(GRangesList(peaks))
  consensus_peaks <- reduce(all_peaks)

  peak_sizes <- data.frame(
    size = width(consensus_peaks)
  )

  group_colour <- colours[group]

  ggplot(peak_sizes, aes(x = size)) +
    geom_density(fill = group_colour, color = group_colour, alpha = 0.6) +
    ylab("Density") +
    xlab("Peak size (nt)") +
    scale_x_continuous(labels = comma, limits = c(0, 2000)) +
    ggtitle("Width of the consensus open chromatin regions") +
    theme_light() +
    theme(
      plot.title = element_text(size = fontsize, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = fontsize),
      axis.title = element_text(size = fontsize),
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.position = legend_position
    )
}


###########################################
# Full report per group
###########################################

make_report_per_group <- function(group) {
  log_message("Preparing detailed ATAC-seq QC report for ", group)

  reads <- plot_reads(mapping_summary, group)
  mapping <- plot_mapped(mapping_summary, group)
  duplicate <- plot_duplicated(mapping_summary, group)
  mitochondrial <- plot_MT(mapping_summary, group)
  peak <- plot_peak_calling(mapping_summary, group)
  frip <- plot_frip(mapping_summary, group)
  insert <- plot_insert_size(group)
  annotation <- plot_anno(group)

  overlap <- get_overlap(group)
  venn_upset <- plot_venn(group, overlap)
  summary_table <- get_summary_table(overlap)
  size_plot <- get_peak_size_dist(group)

  log_message("Assembling report figure for ", group)

  title <- ggdraw() +
    draw_label(group, fontface = "bold", size = 16)

  row_1 <- plot_grid(
    reads,
    mapping,
    duplicate,
    labels = c("A", "B", "C"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 3
  )

  row_2 <- plot_grid(
    mitochondrial,
    peak,
    frip,
    labels = c("D", "E", "F"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 3
  )

  row_3 <- plot_grid(
    plotlist = insert,
    labels = c("G", rep("", max(length(insert) - 1, 0))),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = length(insert)
  )

  row_4 <- plot_grid(
    annotation,
    venn_upset,
    labels = c("H", "I"),
    label_size = 12,
    ncol = 2,
    rel_widths = c(1, 2.5)
  )

  row_5 <- plot_grid(
    summary_table,
    size_plot,
    labels = c("J", "K"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 2
  )

  figure <- plot_grid(
    title,
    row_1,
    row_2,
    row_3,
    row_4,
    row_5,
    ncol = 1,
    rel_heights = c(0.3, 1, 1, 1, 1.5, 1)
  )

  pdf_file <- find_output_file(output_pdf, group, "pdf")
  png_file <- find_output_file(output_png, group, "png")

  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

  log_message("Saving PDF report: ", pdf_file)

  save_plot(
    pdf_file,
    figure,
    base_width = 35,
    base_height = 31,
    units = "cm",
    dpi = 300,
    device = cairo_pdf
  )

  log_message("Saving PNG report: ", png_file)

  save_plot(
    png_file,
    figure,
    base_width = 35,
    base_height = 31,
    units = "cm",
    dpi = 300,
    bg = "white"
  )

  log_message("Completed detailed ATAC-seq QC report for ", group)

  invisible(TRUE)
}


###########################################
# Run analysis
###########################################

log_message("Starting detailed ATAC-seq QC report generation.")

lapply(group_list, make_report_per_group)

log_message("Detailed ATAC-seq QC report generation completed.")