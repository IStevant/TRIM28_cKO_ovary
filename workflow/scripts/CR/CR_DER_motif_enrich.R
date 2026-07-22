source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("monaLisa")
  library("BSgenome.Mmusculus.UCSC.mm10")
  library("doParallel")
  library("foreach")
  library("dplyr")
  library("cowplot")
  library("grid")
  library("ggplot2")
  library("gtable")
  library("SummarizedExperiment")
  library("ComplexHeatmap")
  library("circlize")
  library("motifStack")
  library("randomcoloR")
  library("seriation")
})

doParallel::registerDoParallel(cores = 24)

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

# filtered_DERs
filtered_DERs <- read.csv(file = snakemake@input[["sig_DERs"]], header = TRUE, row.names = 1)
# Load the mouse TFs list
TFs <- as.vector(read.csv(snakemake@input[["TF_genes"]], header = FALSE)[, 1])
# Run analysis using the genome bakground or calculating the enrichment compared to the conditions
background <- snakemake@params[["background"]]
res_table <- snakemake@output[["res_table"]]
# Print the logos of the TFs on the heatmap
logos <- snakemake@params[["logos"]]
genome_version <- "mm10"

samplesheet <- read.csv(file = snakemake@input[["samplesheet"]], row.names = 1)

conditions <- unique(samplesheet$conditions)

conditions_color <- distinctColorPalette(length(unique(samplesheet$conditions)))
names(conditions_color) <- sort(unique(samplesheet$conditions))

# Load Jaspar 2024 DEtabase
JASPAR <- JASPAR2020::JASPAR2020
JASPAR@db <- JASPAR2024::JASPAR2024() %>% .@db

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################

get_midpoint_gr <- function(peaks) {
  mids <- floor((start(peaks) + end(peaks))/2)
  GRanges(seqnames=seqnames(peaks),
          ranges=IRanges(start=mids, end=mids),
          strand="*")
}

#' Get TFBS motifs enrichment using the monaLisa package.
#' @param DERs Table of the differentially accessible regions, with the peak coordinates as rownames.
#' @param res_table Minimum value. Default is 5.
#' @return Return a monaLisa enrichment object.
get_enriched_TFs <- function(DERs, res_table) {

DERs <- filtered_DERs

  # genes <- TFs

  pwms <- TFBSTools::getMatrixSet(
    JASPAR,
    opts = list(
      matrixtype = "PWM",
      tax_group = "vertebrates"
    )
  )


  peaks <- DERs
  clusters <- peaks$x
  # seqlevels_keep <- paste0("chr", c(1:19, "X"))
  # peak_gr <- keepSeqlevels(peak_gr, seqlevels_keep, pruning.mode = "coarse")
  # summits <- get_midpoint_gr(peak_gr)
  # peak_gr <- resize(summits, width = 3000, fix = "center")

  # generate GRanges objects
  peak_gr <- GRanges(rownames(peaks))
  # summits <- get_midpoint_gr(peak_gr)
  # peak_gr <- resize(summits, width = 1000, fix = "center")

  # Get the peak sequences
  if (genome_version == "mm10") {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10, peak_gr)
  } else {
    sequences <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm39, peak_gr)
  }

  # Define which sequences are male or female specific
  bins <- clusters
  bins <- factor(bins)

  # Calculate motif enrichments
  # If background is "genome", run the analysis against radom genomic regions
  # Else, compare the conditions
  if (background == "genome") {
    if (genome_version == "mm10") {
      se2 <- calcBinnedMotifEnrR(
        seqs = sequences,
        bins = bins,
        pwmL = pwms,
        background = "genome",
        genome = BSgenome.Mmusculus.UCSC.mm10,
        genome.regions = NULL, # sample from full genome
        genome.oversample = 2,
        BPPARAM = BiocParallel::MulticoreParam(24)
      )
    } else {
      se2 <- calcBinnedMotifEnrR(
        seqs = sequences,
        bins = bins,
        pwmL = pwms,
        background = "genome",
        genome = BSgenome.Mmusculus.UCSC.mm39,
        genome.regions = NULL, # sample from full genome
        genome.oversample = 2,
        BPPARAM = BiocParallel::MulticoreParam(24)
      )
    }
  } else {
    se2 <- calcBinnedMotifEnrR(
      seqs = sequences,
      bins = bins,
      pwmL = pwms,
      BPPARAM = BiocParallel::MulticoreParam(24)
    )
  }

  if (length(as.vector(rowData(se2)$motif.name)) > 1){
    # Get gene expression of the TFs
    TF_names <- tolower(as.vector(rowData(se2)$motif.name))
    names(TF_names) <- rowData(se2)$motif.id

    # for (TF in TF_names) {
    #   if (TF %in% rownames(genes)) {
    #     TF_exp[TF, ] <- unlist(genes[TF, ])
    #   }
    # }

      # Select the motifs enriched with a -Log10Padj > 10
      sel2 <- apply(
        SummarizedExperiment::assay(se2, "negLog10Padj"), 1,
        function(x) max(abs(x), 0, na.rm = TRUE)
      ) >= 10
      seSel <- se2[sel2, ]

    TF_summary <- data.frame(
      TF.name = seSel@elementMetadata$motif.name,
      TF.matrix = rownames(SummarizedExperiment::assay(seSel, "log2enr")),
      SummarizedExperiment::assay(seSel, "log2enr"),
      10^(-SummarizedExperiment::assay(seSel, "negLog10Padj"))
    )

    write.table(TF_summary, file = res_table, row.names=FALSE, quote=FALSE, sep="\t")

  } else {
    write.table("No TFBS enrichment found", file = res_table, row.names=FALSE, quote=FALSE)
  }

  return(seSel)
}

#' Merge the enrichment result by TFBS motif similarity and plot the results as heatmap.
#' @param seSel monaLisa enrichment object.
#' @return Return a grid object.
merge_TF_motifs <- function(seSel) {
  # Cluster motifs by enrichment
  TF_enrichment <- SummarizedExperiment::assay(seSel, "log2enr")
  TF_pVal <- 10^(-SummarizedExperiment::assay(seSel, "negLog10Padj"))
  colnames(TF_pVal) <- paste0("p-val_", colnames(TF_pVal))
  TF_enrichment <- cbind(TF_enrichment, TF_pVal)

  print(TF_enrichment)

  nbCluster <- 1

  hcl <- hclust(dist(TF_enrichment), method = "ward.D2")
  clustering <- cutree(hcl, k = nbCluster)

  motif_sig <- lapply(1:nbCluster, function(cl) {
    TFs <- names(clustering[clustering == cl])
    if (length(TFs) > 1) {
      matrices <- seSel@elementMetadata$motif.pfm[TFs]
      pfms <- universalmotif::convert_motifs(matrices, class = "motifStack-pfm")
      TF_names <- make.unique(unlist(lapply(pfms, function(x) x$name)), sep = "XXX")
      TF_names <- gsub("::", "YYY", TF_names)
      TF_names <- gsub("-", "ZZZ", TF_names)
      pfms <- lapply(seq_along(pfms), function(i) {
        pfms[[i]]$name <- TF_names[[i]]
        return(pfms[[i]])
      })
      names(pfms) <- TF_names
      # pfms <- pfms[unique(names(pfms))]
      print("Cluster motifs")
      hc <- motifStack::clusterMotifs(pfms)
      print("Build the tree")
      phylog <- ade4::hclust2phylog(hc)
      print("Order the tree")
      leaves <- names(phylog$leaves)
      pfms <- pfms[leaves]
      # extract the motif signatures
      print("Build the signatures")
      motifSig <- motifSignature(pfms, phylog, cutoffPval = 0.0001, min.freq = 1)
      ## get the signatures from object of motifSignature
      print("Get the signatures")
      sig <- signatures(motifSig)
    } else {
      sig <- universalmotif::convert_motifs(seSel@elementMetadata$motif.pfm[TFs], class = "motifStack-pcm")
    }
    return(sig)
  })

  motif_sig_enr <- lapply(motif_sig, function(cl) {
    lapply(cl, function(motif) {
      TF_names <- motif@name
      # TF_names <- gsub("XXX", ".", TF_names)
      # TF_names <- gsub("YYY", "::", TF_names)
      # TF_names <- gsub("ZZZ", "-", TF_names)
      lapply(TF_names, function(TF) {
        tf_vector <- toupper(unlist(strsplit(TF, ";")))
        enrichment <- TF_enrichment
        enrichment_names <- make.unique(toupper(seSel@elementMetadata$motif.name))
        enrichment_names <- gsub("[.]", "XXX", enrichment_names)
        enrichment_names <- gsub("::", "YYY", enrichment_names)
        enrichment_names <- gsub("-", "ZZZ", enrichment_names)
        rownames(enrichment) <- enrichment_names
        tf_enrichment <- enrichment[tf_vector, , drop = FALSE]

        if (nrow(tf_enrichment) > 1) {
          extract_prefix <- function(gene) {
            prefix <- sub("([a-zA-Z]+).*", "\\1", gene)
            if (prefix == "NR") {
              return(gene)
            } else if (grepl("^HOX", prefix)) {
              prefix <- stringr::str_sub(prefix, end = -2)
              return(prefix)
            } else {
              return(prefix)
            }
          }

          grouped_genes <- split(tf_vector, sapply(tf_vector, extract_prefix))

          new_tf_vector <- sapply(grouped_genes, function(gene_group) {
            common_prefix <- extract_prefix(gene_group[1])
            if (common_prefix != "HOX") {
              suffixes <- gsub(paste0("^", common_prefix), "", gene_group)
              numeric_suffixes <- sort(as.numeric(suffixes[grepl("^\\d+$", suffixes)]), na.last = TRUE)
              non_numeric_suffixes <- sort(suffixes[!grepl("^\\d+$", suffixes)])
              all_suffixes <- c(non_numeric_suffixes, numeric_suffixes)
              concatenated_suffixes <- paste(all_suffixes, collapse = "/")
              paste0(common_prefix, concatenated_suffixes)
            } else {
              paste0(common_prefix, "s")
            }
          })

          new_tf_vector <- new_tf_vector[order(new_tf_vector)]
          paste(new_tf_vector, collapse = ";")
          tf_names <- paste(new_tf_vector, collapse = ";")

          tf_enrichment <- as.data.frame(colMedians(tf_enrichment))
          colnames(tf_enrichment) <- tf_names
          tf_enrichment <- as.data.frame(t(tf_enrichment))
          tf_enrichment$mat_names <- TF
        } else {
          tf_enrichment <- as.data.frame(tf_enrichment)
          tf_enrichment$mat_names <- TF
        }

        return(tf_enrichment)
      })
    })
  })

  print("Finish merging motifs")
  enrichment <- do.call("rbind", do.call("rbind", do.call("rbind", motif_sig_enr)))
  enrichment <- enrichment[!duplicated(enrichment), ]


  merged_pfms <- do.call("c", do.call("c", do.call("c", motif_sig)))
  names(merged_pfms) <- unlist(lapply(merged_pfms, function(mat) mat@name))

  motifs <- merged_pfms[enrichment$mat_names]
  motifs_pfms <- universalmotif::convert_motifs(motifs, class = "TFBSTools-PFMatrix")

  # # Save merged motifs
  # fileConn <- file(snakemake@output[["res_table"]], "w")

  # for (i in seq_along(motifs_pfms)) {
  #   pfm <- motifs_pfms[[i]]
  #   # tf_name <- pfm@name
  #   tf_name <- rownames(enrichment)[i]
  #   writeLines(paste0(">", tf_name), fileConn)
  #   pfm_matrix <- pfm@profileMatrix
  #   bases <- rownames(pfm_matrix)
  #   for (j in 1:nrow(pfm_matrix)) {
  #     pfm_values <- paste(pfm_matrix[j, ], collapse = "   ")
  #     writeLines(paste(bases[j], "[", pfm_values, "]"), fileConn)
  #   }
  # }

  maxwidth <- max(unlist(lapply(motifs_pfms, function(x) ncol(x@profileMatrix))))

  print("Prepare heatmap")

  grobL <- lapply(motifs_pfms, seqLogoGrob, xmax = maxwidth, xjust = "center")

  enrichment_names <- rownames(enrichment)

  enrichment_names <- gsub("XXX", ".", enrichment_names)
  enrichment_names <- gsub("YYY", "::", enrichment_names)
  enrichment_names <- gsub("ZZZ", "-", enrichment_names)

  rownames(enrichment) <- enrichment_names
  names(grobL) <- enrichment_names

  matrix <- enrichment[, grep("p-val", colnames(enrichment), invert=TRUE)]
  matrix <- matrix[, 1: ncol(matrix)-1]

  matrix[matrix > 1] <- 1
  matrix[matrix < (-1)] <- (-1)


  hmSeqlogo <- HeatmapAnnotation(
    logo = annoSeqlogo(
      grobL = grobL, which = "row",
      space = unit(0.5, "mm"),
      width = unit(1.5, "inch")
    ),
    show_legend = FALSE,
    show_annotation_name = FALSE,
    which = "row"
  )

  bincols <- rep("black", ncol(matrix))
  names(bincols) <- colnames(matrix)
  conditions <- colnames(matrix)

  condition_anno <- HeatmapAnnotation(
    Stages = anno_block(
      gp = gpar(fill = bincols, col = 0),
      labels = names(bincols),
      labels_gp = gpar(col = "white", fontsize = 14),
      height = unit(7, "mm")
    )
  )

  cold <- colorRampPalette(c("#04bbc6", "#52d0cf", "#8addd8", "#c1f0e0", "#fffee8"))
  warm <- colorRampPalette(c("#fffee8", "#ffd9cb", "#ffb1ad", "#f5808d", "#eb2d62"))
  TYP <- c(cold(12), warm(12))

  mypalette <- colorRamp2(
    breaks = seq(-1, 1, length.out = 24),
    colors = TYP
  )

  if (nrow(matrix) < 10) {
    cell_height <- 10
  } else {
    cell_height <- 6
  }

  print("Order heatmap rows")
  # matrix_noNA <- matrix
  # row_dend <- hclust(dist(matrix_noNA), method = "ward.D2")
  # row_dend <- reorder(row_dend, dist(matrix_noNA), method = "OLO")
  # matrix <- matrix[order(-matrix[,1], matrix[,2]), ]


  if (background == "genome") {
    print("Replace non sig. enrichment by NA")
    pVal <- enrichment[,grep("p-val", colnames(enrichment))]
    matrix[pVal>0.01] <- -10
    matrix <- matrix[order(-matrix[,1], matrix[,2]), ]
    matrix[matrix == -10] <- NA
  } else {
    matrix <- matrix[order(-matrix[,1], matrix[,2]), ]
  }

  print("Draw heatmap")

  ht_list <- Heatmap(
    matrix,
    # clustering_method_rows = "ward.D2",
    name = "Log2 enrichment",
    right_annotation = hmSeqlogo,
    top_annotation = condition_anno,
    column_split = conditions,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_title = NULL,
    row_title = NULL,
    col = mypalette,
    width = ncol(matrix) * unit(15, "mm"),
    height = nrow(matrix) * unit(cell_height, "mm"),
    row_names_max_width = unit(25, "cm"),
    heatmap_legend_param = list(direction = "horizontal")
  )

  ht_list_2 <- grid.grabExpr(draw(ht_list, heatmap_legend_side = "bottom"), wrap.grobs = TRUE)

  return(ht_list_2)
}

###########################################
#                                         #
#                Analysis                 #
#                                         #
###########################################

# For each condition, get TFBS motif enrichments
enrichments <- get_enriched_TFs(filtered_DERs, res_table)

if ( length(enrichments@elementMetadata$motif.name) > 5 ) {
  figure <- merge_TF_motifs(enrichments)
} else {
  figure <- ggplot() + theme_void()
}


###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 40,
  base_height = 30,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 40,
  base_height = 30,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
