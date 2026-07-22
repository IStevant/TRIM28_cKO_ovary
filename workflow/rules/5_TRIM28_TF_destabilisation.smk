"""
Author: Isabelle Stévant
Affiliation: CNRS UMR 9002
Date: 15/05/2026
Licence: MIT
"""

configfile: "analysis_parameters.yaml"


# =============================================================================
# Configuration
# =============================================================================

ATAC_PROCESSED_DATA = "results/ATAC/processed_data"
ATAC_TABLES         = "results/ATAC/tables"
ATAC_DAR_TABLES     = "results/Result_3/tables"

CHIP_PROCESSED_DATA = "results/ChIP/processed_data"

OUTPUT_TABLES  = "results/Result_5/tables"
PROCESSED_DATA = "results/Result_5/processed_data"
OUTPUT_PNG     = "results/Result_5/graphs/PNG"
OUTPUT_PDF     = "results/Result_5/graphs/PDF"
TMP            = "results/Result_5/tmp"

SUMO = ["SUMO1", "SUMO2"]

ATAC_NORM_BIGWIG_FOLDER   = f"{ATAC_PROCESSED_DATA}/ATAC_bigwig"
ATAC_MERGED_BIGWIG_FOLDER = f"{ATAC_NORM_BIGWIG_FOLDER}/merged"

CR_NORM_BIGWIG_FOLDER   = f"{CR_PROCESSED_DATA}/CR_bigwig"
CR_MERGED_BIGWIG_FOLDER = f"{CR_NORM_BIGWIG_FOLDER}/merged"

CHIP_NORM_BIGWIG_FOLDER   = f"{CHIP_PROCESSED_DATA}/ChIP_bigwig"
CHIP_MERGED_BIGWIG_FOLDER = f"{CHIP_NORM_BIGWIG_FOLDER}/merged"

GENOME = "results/data/gencode.vM25.annotation.gtf.gz"


# =============================================================================
# Shared inputs
# =============================================================================

TF_INPUTS = {
    "FOXL2": config["FOXL2_Ctrl"],
    "NR5A2": config["NR5A2"],
    "ESR2": config["ESR2"],
    "RUNX": config["RUNX"],
}


# =============================================================================
# Main output list
# =============================================================================

rule_Result5_input_list = [
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_overlap_regions.csv",
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_summary.csv",
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_fisher_summary.csv",
    f"{OUTPUT_PDF}/Fig5A_TF_TRIM28_overlap_DAR.pdf",
    f"{OUTPUT_PNG}/Fig5A_TF_TRIM28_overlap_DAR.png",

    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_overlap_regions_heatmap.csv",
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_heatmap_table.csv",
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_combination_summary_heatmap.csv",
    f"{OUTPUT_PDF}/SupFig_TF_TRIM28_overlap_DAR_heatmap.pdf",
    f"{OUTPUT_PNG}/SupFig_TF_TRIM28_overlap_DAR_heatmap.png",

    f"{OUTPUT_TABLES}/SupFig_DAR_TRIM28_vs_no_TRIM28.csv",
    f"{OUTPUT_PDF}/SupFig_DAR_TRIM28_vs_no_TRIM28.pdf",
    f"{OUTPUT_PNG}/SupFig_DAR_TRIM28_vs_no_TRIM28.png",

    f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_regions.csv",
    f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_genes.csv",
    f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_gene_region_links.csv",
    f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_gene_summary.csv",
    f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_fisher_summary.csv",
    f"{OUTPUT_PDF}/Fig5D_TRIM28_hub_DEG.pdf",
    f"{OUTPUT_PNG}/Fig5D_TRIM28_hub_DEG.png",

    f"{OUTPUT_PDF}/Fig5E_TRIM28_hub_SUMO_enriched_heatmap.pdf",
    f"{OUTPUT_PNG}/Fig5E_TRIM28_hub_SUMO_enriched_heatmap.png",

    f"{OUTPUT_PNG}/Fig5E_TRIM28_hub_SUMO_DAR_enriched_heatmap.png",
    f"{OUTPUT_PDF}/Fig5E_TRIM28_hub_SUMO_DAR_enriched_heatmap.pdf"
]


# =============================================================================
# Figure 5A - Ovarian TF / TRIM28 overlap in DARs
# =============================================================================

rule Fig5A_Ovarian_TFs_overlap_DAR_bg:
    input:
        ATAC_all        = f"{ATAC_TABLES}/ATAC_all_OCRs.bed",
        ATAC_clustering = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28          = config["TRIM28_peaks"],
        **TF_INPUTS
    output:
        regions        = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_overlap_regions.csv",
        class_summary  = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_summary.csv",
        fisher_summary = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_fisher_summary.csv",
        pdf            = f"{OUTPUT_PDF}/Fig5A_TF_TRIM28_overlap_DAR.pdf",
        png            = f"{OUTPUT_PNG}/Fig5A_TF_TRIM28_overlap_DAR.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_5/Fig5A_TF_TRIM28_colocalisation_DAR.R"


# =============================================================================
# Supplementary figure - TF / TRIM28 overlap heatmap in DARs
# =============================================================================

rule SupFig_Ovarian_TFs_overlap_DAR_heatmap:
    input:
        ATAC_clustering = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28          = config["TRIM28_peaks"],
        **TF_INPUTS
    output:
        regions             = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_overlap_regions_heatmap.csv",
        heatmap_table       = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_heatmap_table.csv",
        combination_summary = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/DAR_combination_summary_heatmap.csv",
        pdf                 = f"{OUTPUT_PDF}/SupFig_TF_TRIM28_overlap_DAR_heatmap.pdf",
        png                 = f"{OUTPUT_PNG}/SupFig_TF_TRIM28_overlap_DAR_heatmap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_5/SupFig_TF_TRIM28_colocalisation_DAR_heatmap.R"


# =============================================================================
# Supplementary figure - DARs with or without TRIM28
# =============================================================================

rule SupFig_DAR_TRIM28_vs_no_TRIM28:
    input:
        ATAC_all        = f"{ATAC_TABLES}/ATAC_all_OCRs.bed",
        ATAC_clustering = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28          = config["TRIM28_peaks"],
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/SupFig_DAR_TRIM28_vs_no_TRIM28.csv",
        pdf   = f"{OUTPUT_PDF}/SupFig_DAR_TRIM28_vs_no_TRIM28.pdf",
        png   = f"{OUTPUT_PNG}/SupFig_DAR_TRIM28_vs_no_TRIM28.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_5/SupFig_DAR_TRIM28_vs_no_TRIM28.R"


# =============================================================================
# Figure 5D - TRIM28 hub and differentially expressed genes
# =============================================================================

rule Fig5D_TRIM28_hub_DEG:
    input:
        DEG             = config["scRNAseq"],
        expressed_genes = config["expressed_genes"],
        ATAC_clustering = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28          = config["TRIM28_peaks"],
        genome          = GENOME,
        **TF_INPUTS
    output:
        regions           = f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_regions.csv",
        genes             = f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_genes.csv",
        gene_region_links = f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_gene_region_links.csv",
        gene_summary      = f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_gene_summary.csv",
        fisher_summary    = f"{OUTPUT_TABLES}/Fig5D_TRIM28_hub_DEG_fisher_summary.csv",
        pdf               = f"{OUTPUT_PDF}/Fig5D_TRIM28_hub_DEG.pdf",
        png               = f"{OUTPUT_PNG}/Fig5D_TRIM28_hub_DEG.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_5/Fig5D_TRIM28_hub_DEG.R"


# =============================================================================
# Figure 5E - TRIM28+TF hub / SUMO enriched heatmap
# =============================================================================

rule Fig5E_hub_SUMO_enriched_heatmap:
    input:
        DEG             = config["scRNAseq"],
        TRIM28          = config["TRIM28_peaks"],
        ATAC_clustering = f"{ATAC_DAR_TABLES}/DAR_enriched_heatmap_clustering.csv",
        ATAC_bw_res     = f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt",
        ChIP_bw_res     = expand(f"{CHIP_MERGED_BIGWIG_FOLDER}/{{SUMO}}_merge_bw.txt", SUMO=SUMO),
        samplesheet     = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        sig_DARs        = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        DAR_anno        = f"{ATAC_DAR_TABLES}/ATAC_DAR.tsv",
        genome          = GENOME,
        **TF_INPUTS
    params:
        ATAC_bw = ATAC_MERGED_BIGWIG_FOLDER,
        ChIP_bw = CHIP_MERGED_BIGWIG_FOLDER,
        H3K27ac_bw = CR_MERGED_BIGWIG_FOLDER,
        SUMO    = SUMO
    output:
        pdf = f"{OUTPUT_PDF}/Fig5E_TRIM28_hub_SUMO_enriched_heatmap.pdf",
        png = f"{OUTPUT_PNG}/Fig5E_TRIM28_hub_SUMO_enriched_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_5/Fig5E_TRIM28_hub_DEG_SUMO_heatmap.R"


rule Fig5E_hub_SUMO_DAR_enriched_heatmap:
    input:
        ATAC_clustering = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28          = config["TRIM28_peaks"],
        samplesheet     = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        genome          = GENOME,
        **TF_INPUTS
    params:
        ATAC_bw = ATAC_MERGED_BIGWIG_FOLDER,
        ChIP_bw = CHIP_MERGED_BIGWIG_FOLDER,
        H3K27ac_bw = CR_MERGED_BIGWIG_FOLDER,
        SUMO    = SUMO
    output:
        pdf = f"{OUTPUT_PDF}/Fig5E_TRIM28_hub_SUMO_DAR_enriched_heatmap.pdf",
        png = f"{OUTPUT_PNG}/Fig5E_TRIM28_hub_SUMO_DAR_enriched_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_5/Fig5E_TRIM28_hub_DAR_SUMO_heatmap.R"