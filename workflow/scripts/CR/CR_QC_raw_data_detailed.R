source(".Rprofile")

###########################################
# Load libraries
###########################################

suppressPackageStartupMessages({
  library("ChIPpeakAnno")
  library("cowplot")
  library("data.table")
  library("eulerr")
  library("GenomicRanges")
  library("ggplot2")
  library("gridExtra")
  library("IRanges")
  library("randomcoloR")
  library("rtracklayer")
  library("scales")
  library("shades")
  library("stringr")
  library("txdbmaker")
  library("UpSetR")
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
raw_peaks_dir <- snakemake@input[["raw_peaks"]]
consensus_count_file <- snakemake@input[["consensus_peaks"]]
consensus_peaks_dir <- dirname(consensus_count_file)
repeat_masker_file <- snakemake@input[["repeatMasker"]]
genome_file <- snakemake@input[["genome"]]

mapping_dir <- snakemake@params[["mapping_folder"]]
ctrl_sample <- snakemake@params[["ctrl_AB"]]
min_peak_rep <- snakemake@params[["min_peak_rep"]]
promoter <- snakemake@params[["promoter"]]

output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]


###########################################
# Global plotting parameters
###########################################

ratio <- 0.4
fontsize <- 12
text_size <- 3.8
legend_position <- "none"

options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)


###########################################
# Helper functions
###########################################

read_bed_as_granges <- function(file) {
  bed <- read.table(
    file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )

  GRanges(
    seqnames = bed[[1]],
    ranges = IRanges(
      start = bed[[2]],
      end = bed[[3]]
    ),
    strand = "*"
  )
}


find_output_file <- function(files, group, extension) {
  expected_name <- paste0("CR_QC_detailed_", group, ".", extension)
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


get_condition_files <- function(group, type = c("raw", "consensus_peaks", "consensus_domains")) {
  type <- match.arg(type)

  if (type == "raw") {
    return(list.files(
      path = raw_peaks_dir,
      pattern = paste0(group, ".*bed"),
      full.names = TRUE
    ))
  }

  if (type == "consensus_peaks") {
    return(list.files(
      path = consensus_peaks_dir,
      pattern = paste0(group, ".*peaks.bed"),
      full.names = TRUE
    ))
  }

  if (type == "consensus_domains") {
    return(list.files(
      path = consensus_peaks_dir,
      pattern = paste0(group, ".*domains.bed"),
      full.names = TRUE
    ))
  }
}


###########################################
# Load data
###########################################

log_message("Loading mapping summary table.")

mapping_summary <- read.csv(summary_file)

samples <- mapping_summary$Sample
conditions <- gsub("_R.*", "\\1", samples)
groups <- unique(conditions)
groups <- groups[groups != ctrl_sample]

mapping_summary <- data.frame(
  group = conditions,
  mapping_summary
)

log_message(paste("Detected groups:", paste(groups, collapse = ", ")))

set.seed(1234)
colours <- randomcoloR::distinctColorPalette(length(groups))
names(colours) <- groups

log_message("Creating TxDb object from genome annotation.")
txdb <- txdbmaker::makeTxDbFromGFF(genome_file)


###########################################
# Mapping summary plots
###########################################

plot_reads <- function(mapping_summary, group) {
  log_message(paste("Plotting total reads for", group))

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- stringr::str_to_sentence(sapply(strsplit(data$Sample, "_"), `[`, 3))
  data$label <- scales::comma(data$Total)
  data$Total <- data$Total / 1000000

  ggplot(data, aes(Total, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Million fragment counts") +
    ggtitle("Read pairs per samples\n(fragments)") +
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
  log_message(paste("Plotting mapped reads for", group))

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- stringr::str_to_sentence(sapply(strsplit(data$Sample, "_"), `[`, 3))
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
    ggtitle("% of mapped fragment\n(rel. to trimmed frag.)") +
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
  log_message(paste("Plotting peak counts for", group))

  data <- mapping_summary[mapping_summary$group %in% group, ]
  data$Sample <- stringr::str_to_sentence(sapply(strsplit(data$Sample, "_"), `[`, 3))
  data$label <- scales::comma(data$Peaks)

  ggplot(data, aes(Peaks, Sample, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = label), hjust = 1.1, color = "white", size = text_size) +
    scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 4)) +
    scale_fill_manual(values = colours) +
    scale_y_discrete(limits = rev) +
    ylab("Samples") +
    xlab("Peak counts") +
    ggtitle("Number of called peaks") +
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

plot_insert_size_single <- function(table, sample, group) {
  data <- data.frame(
    insert_size = as.numeric(colnames(table)),
    count = as.numeric(as.vector(t(table[sample, ])))
  )

  title <- stringr::str_to_sentence(paste("Insert size of", sample))
  colour <- colours[group]

  ggplot(data, aes(x = insert_size, y = count)) +
    geom_area(fill = colour) +
    scale_y_continuous(labels = comma) +
    geom_vline(
      xintercept = c(147, 294, 441),
      linetype = "dashed",
      colour = "#666666"
    ) +
    annotate("text", x = 0, y = max(data$count), hjust = 1.1, label = "Nucleosome free", color = "#666666", size = 3.5, angle = 90) +
    annotate("text", x = 167, y = max(data$count), hjust = 1.1, label = "Mononucleosome", color = "#666666", size = 3.5, angle = 90) +
    annotate("text", x = 314, y = max(data$count), hjust = 1.1, label = "Dinucleosome", color = "#666666", size = 3.5, angle = 90) +
    ylab("Coverage") +
    xlab("Insert size (nt)") +
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
  log_message(paste("Plotting insert size distributions for", group))

  file <- file.path(
    mapping_dir,
    "04_reporting",
    "multiqc",
    "multiqc_data",
    "multiqc_fragment_lengths-plot.txt"
  )

  table <- read.table(
    file = file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )

  sample_name <- gsub("_", "", group)
  replicates <- rownames(table[grep(sample_name, rownames(table)), ])

  lapply(replicates, function(rep) {
    plot_insert_size_single(table, rep, group)
  })
}


###########################################
# Genomic annotation plots
###########################################

plot_anno <- function(group) {
  log_message(paste("Plotting genomic annotations for", group))

  raw_peak_files <- get_condition_files(group, "raw")
  consensus_peak_files <- get_condition_files(group, "consensus_peaks")
  consensus_domain_files <- get_condition_files(group, "consensus_domains")

  raw_regions <- lapply(raw_peak_files, read_bed_as_granges)
  consensus_regions <- lapply(c(consensus_peak_files, consensus_domain_files), read_bed_as_granges)

  anno_peaks <- lapply(raw_regions, function(gr) {
    ChIPseeker::annotatePeak(
      gr,
      tssRegion = c(-promoter, 0),
      TxDb = txdb,
      overlap = "all"
    )
  })

  names(anno_peaks) <- sub(".*(R[0-9]).*", "\\1", raw_peak_files)

  anno_consensus <- lapply(consensus_regions, function(gr) {
    ChIPseeker::annotatePeak(
      gr,
      tssRegion = c(-promoter, 0),
      TxDb = txdb,
      overlap = "all"
    )
  })

  names(anno_consensus) <- sub(
    ".*_(domains|peaks)\\.bed$",
    "\\1",
    c(consensus_peak_files, consensus_domain_files)
  )

  anno_peaks <- lapply(seq_along(anno_peaks), function(i) {
    stat <- anno_peaks[[i]]@annoStat
    stat$Conditions <- names(anno_peaks)[i]
    stat
  })

  anno_consensus <- lapply(seq_along(anno_consensus), function(i) {
    stat <- anno_consensus[[i]]@annoStat
    stat$Conditions <- names(anno_consensus)[i]
    stat
  })

  data <- rbindlist(c(anno_peaks, anno_consensus))

  labels <- round(data$Frequency, digits = 2)
  labels[labels < 4] <- 0
  labels <- paste0(labels, "%")
  labels[labels == "0%"] <- " "

  colors <- c(
    "#FB8500",
    "#FFB703",
    "#021547",
    "#2155a0",
    "#219ebc",
    "#8ecae6"
  )

  feature_levels <- c("Intergenic", "Promoter", "3' UTR", "Exon", "Intron", "5' UTR")

  ggplot(
    data,
    aes(
      fill = factor(Feature, levels = feature_levels),
      x = factor(Conditions, levels = unique(Conditions)),
      y = Frequency
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(
        label = labels,
        color = factor(Feature, levels = feature_levels)
      ),
      size = 5,
      position = position_stack(vjust = 0.5)
    ) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = alpha(colors, 0.8)) +
    scale_color_manual(values = c("black", "black", "white", "white", "black", "black"), guide = "none") +
    xlab("") +
    ggtitle("Enriched region\ngenomic location") +
    theme_light() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, margin = margin(r = 10, unit = "pt")),
      strip.text.x = element_text(size = 12, face = "bold")
    )
}


###########################################
# Replicate overlap plots
###########################################

get_overlap <- function(group) {
  log_message(paste("Computing peak overlaps for", group))

  files <- list.files(
    path = raw_peaks_dir,
    pattern = paste0(group, ".*.bed"),
    full.names = TRUE
  )

  peaks <- lapply(files, function(file) {
    read_bed_as_granges(file)
  })

  names(peaks) <- stringr::str_to_sentence(stringr::str_extract(basename(files), "R[0-9]+"))

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  makeVennDiagram(
    Peaks = peaks,
    NameOfPeaks = names(peaks),
    plot = FALSE
  )
}


plot_venn <- function(group, overlap) {
  log_message(paste("Plotting replicate overlap for", group))

  counts <- overlap$vennCounts
  rep_names <- colnames(counts[, -ncol(counts)])

  venn_data <- setNames(
    counts[, "Counts"],
    apply(counts[, -ncol(counts)], 1, function(x) {
      paste(rep_names[x == 1], collapse = "&")
    })
  )

  venn_data <- venn_data[names(venn_data) != ""]

  venn <- euler(venn_data)

  group_colour <- colours[group]
  fill_col <- group_colour %>%
    brightness(0.85) %>%
    saturation(seq(from = 0.2, to = 1, length.out = length(rep_names)))

  border_col <- group_colour %>%
    brightness(0.7) %>%
    saturation(seq(from = 0.2, to = 1, length.out = length(rep_names)))

  venn_plot <- plot(
    venn,
    quantities = list(fontsize = 10),
    edges = list(
      col = as.vector(border_col),
      lex = 2
    ),
    fills = list(
      fill = as.vector(fill_col),
      alpha = 0.45
    )
  )

  title <- ggdraw() +
    draw_label("Peak overlap between replicates", fontface = "bold", size = 12)

  plot_grid(
    title,
    plot_grid(venn_plot, rel_widths = c(0.8)),
    ncol = 1,
    nrow = 2,
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


###########################################
# Peak and domain size distributions
###########################################

plot_region_size_distribution <- function(group, region_type = c("peaks", "domains")) {
  region_type <- match.arg(region_type)

  log_message(paste("Plotting", region_type, "size distribution for", group))

  file_type <- ifelse(region_type == "peaks", "consensus_peaks", "consensus_domains")
  files <- get_condition_files(group, file_type)

  regions <- lapply(files, rtracklayer::import)
  names(regions) <- paste0("Rep", seq_along(files))

  all_regions <- unlist(GRangesList(regions))
  reduced_regions <- reduce(all_regions)

  region_sizes <- data.frame(
    size = width(reduced_regions)
  )

  mean_size <- mean(width(reduced_regions))
  group_colour <- colours[group]

  ggplot(region_sizes, aes(x = size)) +
    geom_density(fill = group_colour, color = group_colour, alpha = 0.6) +
    geom_vline(xintercept = mean_size, colour = "#666666") +
    annotate(
      "text",
      x = mean_size + 500,
      y = max(density(width(reduced_regions))$y) * 0.85,
      hjust = 0,
      label = paste("Mean size:", scales::comma(mean_size), "bp"),
      color = "#444444",
      size = 4.5
    ) +
    ylab("Density") +
    xlab("Region size (bp)") +
    scale_x_continuous(labels = comma, limits = c(0, 15000)) +
    ggtitle(paste("Length of the", region_type)) +
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
# Repeat overlap plots
###########################################

load_repeats <- function() {
  repeats <- read.csv(
    repeat_masker_file,
    sep = "\t",
    header = TRUE
  )

  makeGRangesFromDataFrame(
    repeats,
    keep.extra.columns = TRUE
  )
}


get_regions_for_repeat_analysis <- function(group) {
  raw_peak_files <- get_condition_files(group, "raw")
  consensus_peak_files <- get_condition_files(group, "consensus_peaks")
  consensus_domain_files <- get_condition_files(group, "consensus_domains")

  raw_regions <- lapply(raw_peak_files, read_bed_as_granges)
  names(raw_regions) <- sub(".*(R[0-9]).*", "\\1", raw_peak_files)

  consensus_regions <- lapply(
    c(consensus_peak_files, consensus_domain_files),
    read_bed_as_granges
  )

  names(consensus_regions) <- sub(
    ".*_(domains|peaks)\\.bed$",
    "\\1",
    c(consensus_peak_files, consensus_domain_files)
  )

  c(raw_regions, consensus_regions)
}


plot_repeat_pct <- function(group) {
  log_message(paste("Plotting repeat overlap percentage for", group))

  repeats <- load_repeats()
  regions <- get_regions_for_repeat_analysis(group)
  conditions <- names(regions)

  anno <- lapply(regions, function(gr) {
    hits <- findOverlaps(gr, repeats)
    overlapping_indices <- unique(queryHits(hits))

    percent_overlap <- length(overlapping_indices) / length(gr) * 100
    percent_no_overlap <- 100 - percent_overlap

    data.frame(
      Percentage = c(percent_overlap, percent_no_overlap),
      Feature = c("Repeats", "Other")
    )
  })

  names(anno) <- conditions

  anno <- lapply(seq_along(anno), function(i) {
    stat <- anno[[i]]
    stat$Conditions <- conditions[i]
    stat$labels <- paste0(round(stat$Percentage, digits = 2), "%")
    stat
  })

  data <- rbindlist(anno)

  colors <- c(
    "#3D405B",
    "#E07A5F"
  )

  ggplot(
    data,
    aes(
      fill = Feature,
      x = factor(Conditions, levels = unique(conditions)),
      y = Percentage
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = labels, color = Feature), size = 5, position = position_stack(vjust = 0.5)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = alpha(colors, 0.8)) +
    scale_color_manual(values = c("white", "black"), guide = "none") +
    xlab("Conditions") +
    ggtitle("Overlap with repeats") +
    theme_light() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, margin = margin(r = 10, unit = "pt")),
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.box.spacing = unit(0, "mm")
    )
}


plot_repeat_families <- function(group) {
  log_message(paste("Plotting repeat family overlap for", group))

  repeats <- load_repeats()
  regions <- get_regions_for_repeat_analysis(group)
  conditions <- names(regions)

  anno <- lapply(regions, function(gr) {
    hits <- findOverlaps(gr, repeats)
    overlap_classes <- mcols(repeats)$repClass[subjectHits(hits)]

    peaks_with_class <- data.frame(
      peak = queryHits(hits),
      repClass = overlap_classes
    )

    peaks_with_class_unique <- unique(peaks_with_class)
    count_by_class <- table(peaks_with_class_unique$repClass)
    percent_by_class <- 100 * count_by_class / length(gr)

    data.frame(
      Percentage = as.numeric(percent_by_class),
      Feature = names(percent_by_class)
    )
  })

  names(anno) <- conditions

  anno <- lapply(seq_along(anno), function(i) {
    stat <- anno[[i]]
    stat$Conditions <- conditions[i]
    stat$labels <- paste0(round(stat$Percentage, digits = 2), "%")
    stat
  })

  data <- rbindlist(anno)

  ggplot(
    data,
    aes(
      fill = Feature,
      x = factor(Conditions, levels = unique(conditions)),
      y = Percentage
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = labels), size = 5, position = position_stack(vjust = 0.5), color = "white") +
    scale_y_continuous(labels = scales::comma) +
    xlab("Conditions") +
    ggtitle("Overlap with repeats") +
    theme_light() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12, margin = margin(r = 10, unit = "pt")),
      strip.text.x = element_text(size = 12, face = "bold"),
      legend.box.spacing = unit(0, "mm")
    )
}


###########################################
# Full report per group
###########################################

make_report_per_group <- function(group) {
  log_message(paste("Preparing detailed QC report for", group))

  reads <- plot_reads(mapping_summary, group)
  mapping <- plot_mapped(mapping_summary, group)
  peak <- plot_peak_calling(mapping_summary, group)
  insert <- plot_insert_size(group)

  annotation <- plot_anno(group)

  overlap <- get_overlap(group)
  venn_upset <- plot_venn(group, overlap)
  summary_table <- get_summary_table(overlap)

  peak_size_plot <- plot_region_size_distribution(group, "peaks")
  domain_size_plot <- plot_region_size_distribution(group, "domains")

  repeats <- plot_repeat_pct(group)
  repeats_class <- plot_repeat_families(group)

  log_message(paste("Assembling report figure for", group))

  title <- ggdraw() +
    draw_label(group, fontface = "bold", size = 16)

  row_1 <- plot_grid(
    reads,
    mapping,
    peak,
    labels = c("A", "B", "C"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 3
  )

  row_2 <- plot_grid(
    plotlist = insert,
    labels = c("D", rep("", max(length(insert) - 1, 0))),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = length(insert)
  )

  row_3 <- plot_grid(
    venn_upset,
    summary_table,
    labels = c("E", "F"),
    label_size = 12,
    ncol = 2
  )

  row_4 <- plot_grid(
    peak_size_plot,
    domain_size_plot,
    labels = c("G", "H"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 2
  )

  row_5 <- plot_grid(
    annotation,
    repeats,
    labels = c("I", "J"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 2
  )

  row_6 <- plot_grid(
    repeats_class,
    labels = c("K"),
    align = "h",
    axis = "bt",
    label_size = 12,
    ncol = 1
  )

  figure <- plot_grid(
    title,
    row_1,
    row_2,
    row_3,
    row_4,
    row_5,
    row_6,
    ncol = 1,
    rel_heights = c(0.3, 1, 1, 1, 1.5, 1.5, 1.5)
  )

  pdf_file <- find_output_file(output_pdf, group, "pdf")
  png_file <- find_output_file(output_png, group, "png")

  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

  log_message(paste("Saving PDF report:", pdf_file))

  save_plot(
    pdf_file,
    figure,
    base_width = 35,
    base_height = 50,
    units = "cm",
    dpi = 300,
    device = cairo_pdf
  )

  log_message(paste("Saving PNG report:", png_file))

  save_plot(
    png_file,
    figure,
    base_width = 35,
    base_height = 50,
    units = "cm",
    dpi = 300,
    bg = "white"
  )

  log_message(paste("Completed detailed QC report for", group))

  invisible(TRUE)
}


###########################################
# Run analysis
###########################################

log_message("Starting detailed CUT&RUN QC report generation.")

lapply(groups, make_report_per_group)

log_message("Detailed CUT&RUN QC report generation completed.")