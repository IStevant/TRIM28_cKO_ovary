source(".Rprofile")

suppressPackageStartupMessages({
  library("cowplot")
  library("ggplot2")
  library("randomcoloR")
  library("scales")
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
# Snakemake inputs and outputs
###########################################

log_message("Reading Snakemake inputs and outputs.")

summary_file <- snakemake@input[["summary_counts"]]
merged_peak_file <- snakemake@input[["merged_peaks"]]
peak_count_file <- snakemake@input[["peak_count"]]
domain_count_file <- snakemake@input[["domain_count"]]

output_pdf <- snakemake@output[["pdf_report"]]
output_png <- snakemake@output[["png_report"]]


###########################################
# Validate inputs
###########################################

input_files <- c(
  summary_file,
  merged_peak_file,
  peak_count_file,
  domain_count_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  log_error(
    "The following input files do not exist: ",
    paste(missing_files, collapse = ", ")
  )
}


###########################################
# Load mapping summary
###########################################

log_message("Loading mapping summary: ", summary_file)

mapping_summary <- read.csv(
  summary_file,
  header = TRUE
)

required_columns <- c(
  "Sample",
  "Total",
  "Trimmed",
  "Mapped",
  "Multimapped",
  "Duplicated",
  "Peaks"
)

missing_columns <- setdiff(
  required_columns,
  colnames(mapping_summary)
)

if (length(missing_columns) > 0) {
  log_error(
    "The mapping summary is missing the following columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

samples <- mapping_summary$Sample

conditions <- gsub(
  "_R.*",
  "\\1",
  samples
)

groups <- unique(
  conditions
)

mapping_summary <- data.frame(
  group = conditions,
  mapping_summary
)


###########################################
# Load peak and domain counts
###########################################

log_message("Loading merged peak counts: ", merged_peak_file)

merged_peaks <- read.csv(
  merged_peak_file,
  header = FALSE,
  row.names = 1,
  sep="\t"
)

rownames(merged_peaks) <- gsub(
  "\\.seacr\\.peaks\\.stringent\\.bed$",
  "",
  rownames(merged_peaks)
)

colnames(merged_peaks) <- "Domains"

merged_peaks$Sample <- rownames(
  merged_peaks
)

merged_peaks$group <- gsub(
  "_R.*",
  "\\1",
  rownames(merged_peaks)
)

merged_peaks <- merged_peaks[
  match(
    mapping_summary[-1, "Sample"],
    merged_peaks$Sample
  ),
  ,
  drop = FALSE
]


log_message("Loading consensus peak counts: ", peak_count_file)

consensus_peak_counts <- read.csv(
  peak_count_file,
  header = FALSE,
  row.names = 1
)

rownames(consensus_peak_counts) <- gsub(
  "_consensus_peaks\\.bed$",
  "",
  rownames(consensus_peak_counts)
)

colnames(consensus_peak_counts) <- "Peaks"

consensus_peak_counts$group <- rownames(
  consensus_peak_counts
)

consensus_peak_counts <- consensus_peak_counts[
  match(
    unique(mapping_summary[-1, "group"]),
    consensus_peak_counts$group
  ),
  ,
  drop = FALSE
]


log_message("Loading consensus domain counts: ", domain_count_file)

consensus_domain_counts <- read.csv(
  domain_count_file,
  header = FALSE,
  row.names = 1
)

rownames(consensus_domain_counts) <- gsub(
  "_consensus_domains\\.bed$",
  "",
  rownames(consensus_domain_counts)
)

colnames(consensus_domain_counts) <- "Domains"

consensus_domain_counts$group <- rownames(
  consensus_domain_counts
)

consensus_domain_counts <- consensus_domain_counts[
  match(
    unique(mapping_summary[-1, "group"]),
    consensus_domain_counts$group
  ),
  ,
  drop = FALSE
]


###########################################
# Plot parameters
###########################################

ratio <- 1
fontsize <- 12
text_size <- 4.2
legend_position <- "none"

set.seed(1234)

colours <- randomcoloR::distinctColorPalette(
  length(groups)
)

names(colours) <- groups


###########################################
# Functions
###########################################

#' Plot total fragment counts
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_reads <- function(mapping_summary) {

  log_message("Plotting total fragment counts.")

  data <- mapping_summary

  data$label <- paste0(
    round(
      data$Total / 1000000,
      1
    ),
    "M"
  )

  data$Total <- data$Total / 1000000

  reads <- ggplot(
    data,
    aes(
      x = Total,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
    xlab("Million fragment counts") +
    ggtitle(
      "Read pairs per sample\n(fragments)"
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

  return(reads)
}


#' Plot retained fragments after trimming
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_trimmed <- function(mapping_summary) {

  log_message("Plotting retained fragments after trimming.")

  data <- mapping_summary

  data$label <- paste0(
    round(
      (data$Trimmed * 100) / data$Total,
      2
    ),
    "%"
  )

  data$Trimmed <- (
    data$Trimmed * 100
  ) / data$Total

  trimmed <- ggplot(
    data,
    aes(
      x = Trimmed,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
    xlab("Percentage") +
    ggtitle(
      "% of retained fragments after trimming\n",
      "(relative to total fragments)"
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

  return(trimmed)
}


#' Plot mapped fragment percentages
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_mapped <- function(mapping_summary) {

  log_message("Plotting mapped fragment percentages.")

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

  mapped <- ggplot(
    data,
    aes(
      x = Mapped,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
    xlab("Percentage") +
    ggtitle(
      "% of mapped fragments\n",
      "(relative to trimmed fragments)"
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

  return(mapped)
}


#' Plot multimapped fragment percentages
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_multimapped <- function(mapping_summary) {

  log_message("Plotting multimapped fragment percentages.")

  data <- mapping_summary

  data$label <- paste0(
    round(
      (data$Multimapped * 100) / data$Total,
      2
    ),
    "%"
  )

  data$Multimapped <- (
    data$Multimapped * 100
  ) / data$Total

  multimapped <- ggplot(
    data,
    aes(
      x = Multimapped,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
    xlab("Percentage") +
    ggtitle(
      "% of multimapped fragments\n",
      "(relative to mapped fragments)"
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

  return(multimapped)
}


#' Plot duplicated fragment percentages
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_duplicates <- function(mapping_summary) {

  log_message("Plotting duplicated fragment percentages.")

  data <- mapping_summary

  data$label <- paste0(
    round(
      (data$Duplicated * 100) / data$Total,
      2
    ),
    "%"
  )

  data$Duplicated <- (
    data$Duplicated * 100
  ) / data$Total

  duplicates <- ggplot(
    data,
    aes(
      x = Duplicated,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
    xlab("Percentage") +
    ggtitle(
      "% of duplicated fragments\n",
      "(relative to mapped fragments)"
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

  return(duplicates)
}


#' Plot called peak counts
#'
#' @param mapping_summary Mapping summary table.
#'
#' @return A ggplot object.
plot_peak_calling <- function(mapping_summary) {

  log_message("Plotting called peak counts.")

  data <- mapping_summary

  data$label <- scales::comma(
    data$Peaks
  )

  data <- data[
    -1,
    ,
    drop = FALSE
  ]

  plot_colours <- colours[-1]

  peak_calling <- ggplot(
    data,
    aes(
      x = Peaks,
      y = factor(
        Sample,
        levels = mapping_summary$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
      values = plot_colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Peak counts") +
    ggtitle("Number of called peaks") +
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


#' Plot consensus peak counts
#'
#' @param consensus_peaks Consensus peak count table.
#'
#' @return A ggplot object.
plot_consensus_peaks <- function(consensus_peaks) {

  log_message("Plotting consensus peak counts.")

  data <- consensus_peaks

  data$label <- scales::comma(
    data$Peaks
  )

  plot_colours <- colours[-1]

  consensus_peak_plot <- ggplot(
    data,
    aes(
      x = Peaks,
      y = factor(
        group,
        levels = consensus_peaks$group
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
      values = plot_colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Conditions") +
    xlab("Peak counts") +
    ggtitle("Number of consensus peaks") +
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

  return(consensus_peak_plot)
}


#' Plot enriched domain counts
#'
#' @param merged_peaks Merged peak count table.
#'
#' @return A ggplot object.
plot_merged_peaks <- function(merged_peaks) {

  log_message("Plotting enriched domain counts.")

  data <- merged_peaks

  data$label <- scales::comma(
    data$Domains
  )

  plot_colours <- colours[-1]

  domain_plot <- ggplot(
    data,
    aes(
      x = Domains,
      y = factor(
        Sample,
        levels = merged_peaks$Sample
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
      values = plot_colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Samples") +
    xlab("Domain counts") +
    ggtitle("Number of enriched domains") +
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

  return(domain_plot)
}


#' Plot consensus domain counts
#'
#' @param consensus_domains Consensus domain count table.
#'
#' @return A ggplot object.
plot_consensus_domains <- function(consensus_domains) {

  log_message("Plotting consensus domain counts.")

  data <- consensus_domains

  data$label <- scales::comma(
    data$Domains
  )

  plot_colours <- colours[-1]

  consensus_domain_plot <- ggplot(
    data,
    aes(
      x = Domains,
      y = factor(
        group,
        levels = consensus_domains$group
      ),
      fill = group
    )
  ) +
    geom_bar(
      stat = "identity"
    ) +
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
      values = plot_colours
    ) +
    scale_y_discrete(
      limits = rev
    ) +
    ylab("Conditions") +
    xlab("Domain counts") +
    ggtitle("Number of consensus domains") +
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

  return(consensus_domain_plot)
}


###########################################
# Generate plots
###########################################

reads_plot <- plot_reads(
  mapping_summary
)

trimmed_plot <- plot_trimmed(
  mapping_summary
)

mapped_plot <- plot_mapped(
  mapping_summary
)

multimapped_plot <- plot_multimapped(
  mapping_summary
)

duplicate_plot <- plot_duplicates(
  mapping_summary
)

peak_plot <- plot_peak_calling(
  mapping_summary
)

consensus_peak_plot <- plot_consensus_peaks(
  consensus_peak_counts
)

domain_plot <- plot_merged_peaks(
  merged_peaks
)

consensus_domain_plot <- plot_consensus_domains(
  consensus_domain_counts
)


###########################################
# Assemble figure
###########################################

log_message("Assembling mapping and peak-calling QC figure.")

figure <- cowplot::plot_grid(
  reads_plot,
  trimmed_plot,
  mapped_plot,
  multimapped_plot,
  duplicate_plot,
  peak_plot,
  consensus_peak_plot,
  domain_plot,
  consensus_domain_plot,
  ncol = 3,
  align = "hv",
  axis = "bt",
  labels = "AUTO",
  label_size = 16
)


###########################################
# Export figure
###########################################

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

cowplot::save_plot(
  output_pdf,
  figure,
  dpi = 300,
  base_width = 40,
  base_height = 32,
  units = "cm"
)

log_message("Saving PNG report: ", output_png)

cowplot::save_plot(
  output_png,
  figure,
  dpi = 300,
  bg = "white",
  base_width = 40,
  base_height = 32,
  units = "cm"
)

log_message("Mapping and peak-calling QC report completed.")