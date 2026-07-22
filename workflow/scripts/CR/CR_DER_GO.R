source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("ggrepel")
  library("ggplot2")
  library("doParallel")
  library("foreach")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("cowplot")
})



###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

AB <- snakemake@params[["AB"]]

DER <- read.table(snakemake@input[["DER"]], header = TRUE, row.names = 1, sep="\t")

adj.pval <- snakemake@params[["adjpval"]]
log2FC <- snakemake@params[["log2FC"]]
path <- snakemake@params[["path"]]

n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 12))
doParallel::registerDoParallel(cores = n_cores)

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################


#' Prepare the volcano plot.
#' @param dds DESeq2 analysis result object.
#' @param p.adj adjusted p-value threshold used in the DE analysis.
#' @param log2FC Log2(Fold change) threshold used in the DE analysis.
#' @param path Path to the GO term result table file.
#' @return Ggplot object.
plot_DEG <- function(DER, p.adj, log2FC, path) {

  # print(head(DER))

  res <- mutate(
    DER,
    Diff.Exp. = dplyr::case_when(
      log2FoldChange >= log2FC ~ "Up in WT",
      log2FoldChange <= (-log2FC) ~ "Up in KO",
      TRUE ~ "non sig."
    )
  )

  # sig.DE <- subset(res, padj < p.adj)
  # sig.DE <- subset(sig.DE, abs(log2FoldChange) > log2FC)

  # Select genes overlaped with differentially enriched regions
  overlapped_genes <- res[!grepl("Intergenic | Downstream", res$nearest.gene), ]

  de_genes <- data.frame(
    gene = overlapped_genes$nearest.gene,
    cluster = overlapped_genes$Diff.Exp.
  )
  gene_file <- paste0(path, "/CR_", AB, "_genes_DER.csv")
  write.table(de_genes, file = gene_file, quote = FALSE, sep = "\t", row.names = FALSE)

  res_file <- paste0(path, "/CR_", AB, "_GO_genes_DER.csv")
  GO_terms <- GO_term_per_cluster(de_genes, res_file)
  go_term_plot <- go_plot(GO_terms, nb_terms = 10)
  # KEGG_terms <- KEGG_term_per_cluster(de_genes)
  # if (!is.null(KEGG_terms)){
  # 	KEGG_term_plot <- kegg_plot(KEGG_terms, nb_terms=5)
  # } else {
  # 	KEGG_term_plot <- NULL
  # }

  return(go_term_plot)
}

#' Get GO term enrichment.
#' @param de_genes Vector of gene names.
#' @param res_file Path and name of the result file.
#' @return dataframe.
GO_term_per_cluster <- function(de_genes, res_file) {
  print("Calculate GO term over-representation...")
  formula_res <- compareCluster(
    gene ~ cluster,
    data = de_genes,
    fun = "enrichGO",
    keyType = "SYMBOL",
    OrgDb = "org.Mm.eg.db",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  # print("Calculate GO term semantic similarities...")
  # lineage1_ego <- simplify(
  #   formula_res,
  #   cutoff = 0.6,
  #   by = "p.adjust",
  #   select_fun = min
  # )

  # write.table(lineage1_ego, file = res_file, quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(formula_res, file = res_file, quote = FALSE, sep = "\t", row.names = FALSE)

  # return(lineage1_ego)
  return(formula_res)
}

#' Get KEGG pathway enrichment.
#' @param de_genes Vector of gene names.
#' @return dataframe.
KEGG_term_per_cluster <- function(de_genes) {
  print("Calculate KEGG term over-representation...")
  entrez_genes <- bitr(de_genes$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  de_gene_clusters <- de_genes[de_genes$gene %in% entrez_genes$SYMBOL, c("gene", "cluster")]
  de_gene_clusters <- data.frame(
    ENTREZID = entrez_genes$ENTREZID[which(entrez_genes$SYMBOL == de_gene_clusters$gene)],
    cluster = de_gene_clusters$cluster
  )

  formula_res <- compareCluster(
    ENTREZID ~ cluster,
    data = de_gene_clusters,
    fun = "enrichKEGG",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  )
  return(formula_res)
}

#' Plot GO term enrichment.
#' @param go_res Dataframe containing the GO term enrichment result.
#' @param nb_terms Number of top GO term per cluster to show.
#' @return Ggplot object
go_plot <- function(go_res, nb_terms = 5) {
  # Make the first GO term letter as capital letter
  go_res@compareClusterResult[, 4] <- gsub("^([a-z])", "\\U\\1", go_res[, 4], perl = TRUE)
  options(enrichplot.colours = c("#77BFA3", "#98C9A3", "#BFD8BD", "#DDE7C7", "#EDEEC9"))
  plot <- clusterProfiler::dotplot(go_res, showCategory = nb_terms)
  x_labels <- levels(as.data.frame(go_res)$Cluster)
  # Important: Reload ggplot2 to apply new theme
  library(ggplot2)
  plot <- plot +
    geom_point(shape = 21) +
    ggtitle(paste0("Enriched biological process GO terms (Top ", nb_terms, ")")) +
    scale_x_discrete(labels = x_labels) +
    labs(color = "Adj. p-value", size = "Gene ratio") +
    guides(color = guide_colorbar(reverse = TRUE), size = guide_legend(reverse = TRUE)) +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      axis.title.x = element_blank()
    )
  return(plot)
}

#' Plot GO term enrichment.
#' @param kegg_res Dataframe containing the KEGG enrichment result.
#' @return Ggplot object
kegg_plot <- function(kegg_res, nb_terms = 5) {
  plot <- dotplot(kegg_res)
  return(plot)
}
#################################################################################################################################

###########################################
#                                         #
#         Volcano + GO sex DEGs           #
#                                         #
###########################################

DER_plot <- plot_DEG(DER, adj.pval, log2FC, path)

figure <- plot_grid(
  plotlist = DER_plot,
  labels = "AUTO",
  ncol = 1
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  DER_plot,
  base_width = 20,
  base_height = 25,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  DER_plot,
  base_width = 20,
  base_height = 25,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)
