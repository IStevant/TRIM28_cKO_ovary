source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################

suppressPackageStartupMessages({
  library("ggrepel")
  library("ggplot2")
  library("cowplot")
  library("GenomicRanges")
  library("dplyr")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

AB <- snakemake@params[["AB"]]

WT_peaks <- rtracklayer::import(snakemake@input[["WT_peaks"]])
DER <- read.csv(snakemake@input[["DER"]], header = TRUE, row.names = 1)

repeat_masker <- snakemake@input[["repeatMasker"]]

###########################################
#                                         #
#               Functions                 #
#                                         #
###########################################


peaks_overlap_family <- function(peaks, rmsk) {
  hits <- findOverlaps(peaks, rmsk)
  fam <- mcols(rmsk)$repFamily[subjectHits(hits)]
  tab <- table(fam)
  res <- data.frame(
    repFamily = names(tab),
    n_peaks_overlap = as.integer(tab),
    stringsAsFactors = FALSE
  )
  res$percent_of_peaks <- 100 * res$n_peaks_overlap / length(peaks)
  res <- res[order(res$percent_of_peaks, decreasing=TRUE), ]
  rownames(res) <- NULL
  res
}


family_coverage_by_elements <- function(peaks, rmsk) {
  hits <- findOverlaps(rmsk, peaks)
  covered <- rep(FALSE, length(rmsk))
  covered[queryHits(hits)] <- TRUE
  fam <- mcols(rmsk)$repFamily

  all_fam <- unique(fam)
  n_total <- tapply(fam, fam, length)
  n_covered <- tapply(covered, fam, sum)
  n_covered <- n_covered[all_fam]
  n_total <- n_total[all_fam]
  n_covered[is.na(n_covered)] <- 0

  pct <- 100 * n_covered / n_total
  res <- data.frame(repFamily = all_fam,
                    n_total = n_total,
                    n_covered = n_covered,
                    pct_covered = pct,
                    stringsAsFactors = FALSE)
  res <- res[order(res$pct_covered, decreasing = TRUE), ]
  rownames(res) <- NULL
  res
}


fisher_by_family <- function(target_idx, background_idx, repFamily) {
  fams <- sort(unique(as.character(repFamily)))
  res <- data.frame(repFamily = fams,
                    a = integer(length(fams)),
                    b = integer(length(fams)),
                    c = integer(length(fams)),
                    d = integer(length(fams)),
                    odds = NA_real_, pval = NA_real_, stringsAsFactors = FALSE)
  n_target <- length(target_idx)
  n_background <- length(background_idx)
  
  for(i in seq_along(fams)) {
    fam <- fams[i]
    fam_idx <- which(as.character(repFamily) == fam)
    a <- sum(fam_idx %in% target_idx)
    b <- n_target - a
    bg_idx <- setdiff(background_idx, target_idx)
    c <- sum(fam_idx %in% bg_idx)
    d <- length(bg_idx) - c
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft <- tryCatch(fisher.test(mat, alternative = "two.sided"), error = function(e) NULL) # greater
    if(!is.null(ft)) {
      res$odds[i] <- unname(ft$estimate)
      res$pval[i] <- ft$p.value
    }
    res$a[i] <- a; res$b[i] <- b; res$c[i] <- c; res$d[i] <- d
  }
  res$padj <- p.adjust(res$pval, method = "BH")
  res
}


plot_volcano_TE <- function(res, te_class, title, pval_cutoff = 0.05) {
  res$log2OR <- log2(res$odds)
  res$negLogP <- -log10(res$pval)
  res$class <- te_class[match(res$repFamily, names(te_class))]
  res$sig <- res$padj <= pval_cutoff
  res <- res[is.finite(res$log2OR) & is.finite(res$negLogP), ]
  # res <- res[res$log2OR>0,]
  
  ggplot(res, aes(x = log2OR, y = negLogP, color = class)) +
    geom_point(aes(alpha = sig), size = 2) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = subset(res, sig & negLogP >= -log10(pval_cutoff)),
      aes(label = repFamily),
      size = 3, max.overlaps = 15
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = FALSE) +
    labs(
      title = title,
      x = "log2(Odds Ratio)",
      y = "-log10(p-value)",
      color = "TE class"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}


###########################################
#                                         #
#               Draw plots                #
#                                         #
###########################################


# Load repeats as GRanges
repeats <- read.csv(repeat_masker, sep="\t", header=TRUE)
rmsk <- GenomicRanges::makeGRangesFromDataFrame(
  repeats,
  keep.extra.columns = TRUE
)

print(head(DER))

peaks_lost <- GRanges(rownames(DER[DER$x == "a",, drop = FALSE]))
print(head(peaks_lost))
peaks_gained <- GRanges(rownames(DER[DER$x == "b",, drop = FALSE]))


hits_WT <- findOverlaps(rmsk, WT_peaks)
covered_WT <- rep(FALSE, length(rmsk)); covered_WT[queryHits(hits_WT)] <- TRUE

hits_lost <- findOverlaps(rmsk, peaks_lost)
covered_lost <- rep(FALSE, length(rmsk)); covered_lost[queryHits(hits_lost)] <- TRUE

hits_gained <- findOverlaps(rmsk, peaks_gained)
covered_gained <- rep(FALSE, length(rmsk)); covered_gained[queryHits(hits_gained)] <- TRUE

repFamily <- mcols(rmsk)$repFamily
repClass  <- mcols(rmsk)$repClass
names(repClass) <- as.character(repFamily)
all_idx <- seq_along(repFamily)

res_WT     <- fisher_by_family(which(covered_WT), all_idx, repFamily)
res_lost   <- fisher_by_family(which(covered_lost), all_idx, repFamily)
res_gained <- fisher_by_family(which(covered_gained), all_idx, repFamily)


p1 <- plot_volcano_TE(res_WT, repClass, paste("TE families marked by", AB, "(WT vs genome)"))
p2 <- plot_volcano_TE(res_lost, repClass, paste("TE families losing", AB, "(KO vs genome)"))
p3 <- plot_volcano_TE(res_gained, repClass, paste("TE families gaining", AB, "(KO vs genome)"))



figure <- plot_grid(
  plotlist = list(p1, p2, p3),
  labels = "AUTO",
  ncol = 1, 
  align = "v"
)

###########################################
#                                         #
#               Save files                #
#                                         #
###########################################

save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 30,
  base_height = 30,
  units = c("cm"),
  dpi = 300
)

save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 30,
  base_height = 30,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

