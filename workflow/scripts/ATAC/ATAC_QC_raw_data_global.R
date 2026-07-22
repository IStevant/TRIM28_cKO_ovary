source(".Rprofile")

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
  library("randomcoloR")
  library("scales")
})


###########################################################################
# Logging functions
###########################################################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}


###########################################################################
# Snakemake inputs and outputs
###########################################################################

log_message("Reading Snakemake inputs and outputs.")

summary_file <- snakemake@input[["summary_counts"]]

output_pdf <- snakemake@output[["pdf_report"]]
output_png <- snakemake@output[["png_report"]]


###########################################################################
# Plot parameters
###########################################################################

ratio <- 1.8
fontsize <- 12
text_size <- 4.2
legend_position <- "none"


###########################################################################
# Load data
###########################################################################

log_message("Loading mapping summary table: ", summary_file)

if (!file.exists(summary_file)) {
  log_error("Mapping summary file does not exist: ", summary_file)
}

mapping_summary <- read.csv(summary_file)

required_columns <- c(
  "Sample",
  "Total",
  "Mapped",
  "Duplicated",
  "Mitochondrial",
  "Peaks",
  "FRiP"
)

missing_columns <- setdiff(
  required_columns,
  colnames(mapping_summary)
)

if (length(missing_columns) > 0) {
  log_error(
    "The mapping summary table is missing the following columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

samples <- mapping_summary$Sample
conditions <- gsub("_REP.*", "\\1", samples)
groups <- unique(conditions)

mapping_summary <- data.frame(
  group = conditions,
  mapping_summary
)

log_message(
  "Detected ",
  length(groups),
  " experimental group(s): ",
  paste(groups, collapse = ", ")
)

set.seed(1234)

colours <- randomcoloR::distinctColorPalette(
  length(groups)
)

names(colours) <- groups


###########################################################################
# Plot functions
###########################################################################

# Plot the total number of reads per sample.
plot_reads <- function(mapping_summary) {

  log_message("Plotting total read counts.")

  data <- mapping_summary

  data$label <- paste0(
    round(data$Total / 1000000, 1),
    "M"
  )

  data$Total <- data$Total / 1000000

  reads <- ggplot(
    data,
    aes(
      Total,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label),
      hjust = 1.1,
      color = "white",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Million read counts") +
    ggtitle("Read counts per samples (R1+R2)") +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(reads)
}


# Plot the percentage of mapped reads relative to total reads.
plot_mapped <- function(mapping_summary) {

  log_message("Plotting mapped read percentages.")

  data <- mapping_summary

  data$label <- paste0(
    round(
      (data$Mapped * 100) / data$Total,
      2
    ),
    "%"
  )

  data$Mapped <- (
    data$Mapped * 100
  ) / data$Total

  mapping <- ggplot(
    data,
    aes(
      Mapped,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label),
      hjust = 1.1,
      color = "white",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(0, 100)
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle(
      "% of mapped reads per samples\n",
      "(rel. to total reads)"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(mapping)
}


# Plot the percentage of reads retained after filtering.
plot_filtered <- function(mapping_summary) {

  log_message("Plotting retained read percentages after filtering.")

  required_filter_columns <- c(
    "Mapped",
    "Kept.after.filtering"
  )

  missing_filter_columns <- setdiff(
    required_filter_columns,
    colnames(mapping_summary)
  )

  if (length(missing_filter_columns) > 0) {
    log_error(
      "Cannot plot filtered reads. Missing columns: ",
      paste(missing_filter_columns, collapse = ", ")
    )
  }

  data <- mapping_summary

  data$label1 <- paste0(
    round(
      (
        data$Kept.after.filtering * 100
      ) / data$Mapped,
      2
    ),
    "%"
  )

  data$label2 <- paste0(
    "(",
    round(
      data$Kept.after.filtering / 1000000,
      1
    ),
    "M)"
  )

  data$Kept.after.filtering <- (
    data$Kept.after.filtering * 100
  ) / data$Mapped

  filtered <- ggplot(
    data,
    aes(
      Kept.after.filtering,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label1),
      hjust = 1.1,
      color = "white",
      size = text_size
    ) +
    geom_text(
      aes(label = label2),
      hjust = -0.15,
      color = "#444444",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(0, 100)
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle(
      "% of kept reads after filters\n",
      "(rel. to mapped reads)"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(filtered)
}


# Plot the percentage of duplicated reads relative to mapped reads.
plot_duplicated <- function(mapping_summary) {

  log_message("Plotting duplicated read percentages.")

  data <- mapping_summary

  data$dup <- (
    data$Duplicated * 100
  ) / data$Mapped

  data$label <- paste0(
    round(data$dup, 2),
    "%"
  )

  duplicated <- ggplot(
    data,
    aes(
      dup,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label),
      hjust = -0.15,
      color = "#444444",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(0, 100)
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle(
      "% of duplicated reads\n",
      "(rel. to mapped reads)"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(duplicated)
}


# Plot the percentage of mitochondrial reads relative to mapped reads.
plot_MT <- function(mapping_summary) {

  log_message("Plotting mitochondrial read percentages.")

  data <- mapping_summary

  data$MT <- (
    data$Mitochondrial * 100
  ) / data$Mapped

  data$label <- paste0(
    round(data$MT, 2),
    "%"
  )

  MT <- ggplot(
    data,
    aes(
      MT,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label),
      hjust = -0.15,
      color = "#444444",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(0, 100)
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Percentages") +
    ggtitle(
      "% of mitochondrial reads\n",
      "(rel. to mapped reads)"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(MT)
}


# Plot the number of called peaks per sample.
plot_peak_calling <- function(mapping_summary) {

  log_message("Plotting peak counts.")

  data <- mapping_summary

  data$label <- scales::comma(
    data$Peaks
  )

  peak_calling <- ggplot(
    data,
    aes(
      Peaks,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = label),
      hjust = 1.1,
      color = "white",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Peak counts") +
    ggtitle("Number of peaks per sample") +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(peak_calling)
}


# Plot the fraction of reads in peaks per sample.
plot_frip <- function(mapping_summary) {

  log_message("Plotting FRiP scores.")

  data <- mapping_summary

  frip <- ggplot(
    data,
    aes(
      FRiP,
      factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = round(FRiP, 2)),
      hjust = 1.1,
      color = "white",
      size = text_size
    ) +
    scale_x_continuous(
      labels = comma
    ) +
    scale_fill_manual(
      values = colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("FRiP score") +
    ggtitle("Fraction of reads in peaks per sample") +
    theme_light() +
    theme(
      plot.title = element_text(
        size = fontsize,
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        size = fontsize
      ),
      axis.title = element_text(
        size = fontsize
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = fontsize
      ),
      legend.position = legend_position,
      aspect.ratio = ratio
    )

  return(frip)
}


###########################################################################
# Generate plots
###########################################################################

log_message("Generating mapping summary plots.")

reads <- plot_reads(mapping_summary)
mapping <- plot_mapped(mapping_summary)

# filtered <- plot_filtered(mapping_summary)

duplicate <- plot_duplicated(mapping_summary)
mitochondrial <- plot_MT(mapping_summary)
peak <- plot_peak_calling(mapping_summary)
frip <- plot_frip(mapping_summary)


###########################################################################
# Assemble figure
###########################################################################

log_message("Assembling the mapping summary figure.")

figure <- plot_grid(
  reads,
  mapping,
  # filtered,
  duplicate,
  mitochondrial,
  peak,
  frip,
  ncol = 3,
  align = "hv",
  axis = "bt",
  labels = "AUTO",
  label_size = 16
)


###########################################################################
# Export figure
###########################################################################

dir.create(
  dirname(output_pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(output_png),
  recursive = TRUE,
  showWarnings = FALSE
)

log_message("Saving PDF report: ", output_pdf)

save_plot(
  output_pdf,
  figure,
  dpi = 300,
  base_width = 32,
  base_height = 27,
  units = "cm"
)

log_message("Saving PNG report: ", output_png)

save_plot(
  output_png,
  figure,
  dpi = 300,
  bg = "white",
  base_width = 32,
  base_height = 27,
  units = "cm"
)

log_message("Mapping summary report generation completed.")