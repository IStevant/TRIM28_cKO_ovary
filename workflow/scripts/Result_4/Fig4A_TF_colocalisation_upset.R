source(".Rprofile")

###########################################
# Logging functions
###########################################

log_message <- function(...) {
  message("[INFO] ", ...)
}

log_error <- function(...) {
  stop("[ERROR] ", ...)
}

log_message("Starting ovarian TF UpSet plot script")

###########################################
# Libraries
###########################################

suppressPackageStartupMessages({
  library("GenomicRanges")
  library("S4Vectors")
  library("rtracklayer")
  library("ggplot2")
  library("cowplot")
  library("ComplexHeatmap")
  library("grid")
})

###########################################
# Load input data
###########################################

log_message("Importing TF peak files")

TF_list <- list(
  FOXL2 = rtracklayer::import(snakemake@input[["FOXL2"]]),
  NR5A2 = rtracklayer::import(snakemake@input[["NR5A2"]]),
  ESR2  = rtracklayer::import(snakemake@input[["ESR2"]]),
  RUNX = rtracklayer::import(snakemake@input[["RUNX"]])
)

TF_names <- names(TF_list)

for (TF_name in TF_names) {
  log_message(TF_name, " peaks: ", length(TF_list[[TF_name]]))
}

###########################################
# Functions
###########################################

make_region_id <- function(regions) {
  paste0(
    as.character(GenomicRanges::seqnames(regions)),
    ":",
    GenomicRanges::start(regions),
    "-",
    GenomicRanges::end(regions)
  )
}

build_TF_intersection_matrix <- function(TF_list) {
  log_message("Building merged TF regulatory regions")

  all_TF_regions <- unlist(
    GenomicRanges::GRangesList(TF_list),
    use.names = FALSE
  )

  merged_regions <- GenomicRanges::reduce(
    all_TF_regions
  )

  log_message("Merged TF regions: ", length(merged_regions))

  overlap_matrix <- sapply(
    names(TF_list),
    function(TF_name) {
      GenomicRanges::countOverlaps(
        merged_regions,
        TF_list[[TF_name]],
        ignore.strand = TRUE
      ) > 0
    }
  )

  n_TFs <- rowSums(overlap_matrix)

  combination <- apply(
    overlap_matrix,
    1,
    function(x) {
      TFs <- names(x)[x]

      if (length(TFs) == 0) {
        return("No TF")
      }

      paste(TFs, collapse = " + ")
    }
  )

  region_table <- data.frame(
    chromosome = as.character(GenomicRanges::seqnames(merged_regions)),
    start = GenomicRanges::start(merged_regions),
    end = GenomicRanges::end(merged_regions),
    region = make_region_id(merged_regions),
    overlap_matrix,
    n_TFs = n_TFs,
    TF_combination = combination,
    stringsAsFactors = FALSE
  )

  list(
    regions = merged_regions,
    overlap_matrix = overlap_matrix,
    region_table = region_table
  )
}

make_filtered_upset_matrix <- function(
  overlap_matrix,
  min_TFs = 3
) {
  keep_regions <- rowSums(overlap_matrix) >= min_TFs

  filtered_matrix <- overlap_matrix[
    keep_regions,
    ,
    drop = FALSE
  ]

  if (nrow(filtered_matrix) == 0) {
    log_error("No TF intersections with at least ", min_TFs, " TFs")
  }

  filtered_matrix
}

make_combination_summary <- function(
  region_table,
  min_TFs = 3
) {
  summary_table <- region_table[
    region_table$n_TFs >= min_TFs,
    ,
    drop = FALSE
  ]

  combination_summary <- as.data.frame(
    table(summary_table$TF_combination),
    stringsAsFactors = FALSE
  )

  colnames(combination_summary) <- c(
    "TF_combination",
    "n_regions"
  )

  combination_summary <- combination_summary[
    combination_summary$n_regions > 0,
    ,
    drop = FALSE
  ]

  combination_summary$n_TFs <- vapply(
    strsplit(combination_summary$TF_combination, " \\+ "),
    length,
    numeric(1)
  )

  combination_summary <- combination_summary[
    order(
      combination_summary$n_TFs,
      combination_summary$n_regions,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  combination_summary
}

draw_upset_plot <- function(
  filtered_matrix
) {
  log_message("Drawing UpSet plot for intersections with >=3 TFs")

  combination_matrix <- ComplexHeatmap::make_comb_mat(
    filtered_matrix
  )

  upset <- ComplexHeatmap::UpSet(
    combination_matrix,
    comb_order = order(
      ComplexHeatmap::comb_degree(combination_matrix),
      ComplexHeatmap::comb_size(combination_matrix),
      decreasing = TRUE
    ),
    top_annotation = ComplexHeatmap::upset_top_annotation(
      combination_matrix,
      add_numbers = TRUE
    ),
    right_annotation = ComplexHeatmap::upset_right_annotation(
      combination_matrix,
      add_numbers = TRUE,
    ),
    row_names_side = "left",
    column_title = "Ovarian TF intersections with at least 3 TFs",
    heatmap_legend_param = list(
      title = "Overlap"
    )
  )

  grid::grid.grabExpr(
    ComplexHeatmap::draw(upset)
  )
}

###########################################
# Run analysis
###########################################

TF_intersections <- build_TF_intersection_matrix(
  TF_list = TF_list
)

filtered_matrix <- make_filtered_upset_matrix(
  overlap_matrix = TF_intersections$overlap_matrix,
  min_TFs = 3
)

filtered_region_table <- TF_intersections$region_table[
  TF_intersections$region_table$n_TFs >= 3,
  ,
  drop = FALSE
]

combination_summary <- make_combination_summary(
  region_table = TF_intersections$region_table,
  min_TFs = 3
)

log_message("Regions with >=3 TFs: ", nrow(filtered_region_table))

###########################################
# Save tables
###########################################

write.table(
  filtered_region_table,
  file = snakemake@output[["table"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# write.table(
#   combination_summary,
#   file = snakemake@output[["summary"]],
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )

###########################################
# Save figure
###########################################

upset_plot <- draw_upset_plot(
  filtered_matrix = filtered_matrix
)

cowplot::save_plot(
  filename = snakemake@output[["pdf"]],
  plot = upset_plot,
  base_width = 18,
  base_height = 12,
  units = "cm",
  dpi = 300
)

cowplot::save_plot(
  filename = snakemake@output[["png"]],
  plot = upset_plot,
  base_width = 18,
  base_height = 12,
  units = "cm",
  dpi = 300,
  bg = "white"
)

log_message("Analysis completed successfully")