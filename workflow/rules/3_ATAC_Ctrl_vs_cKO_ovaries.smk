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

CR_PROCESSED_DATA = "results/CR/processed_data"
CR_TABLES         = "results/CR/tables"

CHIP_PROCESSED_DATA = "results/ChIP/processed_data"

OUTPUT_TABLES  = "results/Result_3/tables"
PROCESSED_DATA = "results/Result_3/processed_data"
OUTPUT_PNG     = "results/Result_3/graphs/PNG"
OUTPUT_PDF     = "results/Result_3/graphs/PDF"
TMP            = "results/Result_3/tmp"

AB   = ["H3K9me3", "H3K27ac"]
SUMO = ["SUMO1", "SUMO2"]

ATAC_NORM_BIGWIG_FOLDER   = f"{ATAC_PROCESSED_DATA}/ATAC_bigwig"
ATAC_MERGED_BIGWIG_FOLDER = f"{ATAC_NORM_BIGWIG_FOLDER}/merged"

CR_NORM_BIGWIG_FOLDER   = f"{CR_PROCESSED_DATA}/CR_bigwig"
CR_MERGED_BIGWIG_FOLDER = f"{CR_NORM_BIGWIG_FOLDER}/merged"

CHIP_NORM_BIGWIG_FOLDER   = f"{CHIP_PROCESSED_DATA}/ChIP_bigwig"
CHIP_MERGED_BIGWIG_FOLDER = f"{CHIP_NORM_BIGWIG_FOLDER}/merged"

GENOME       = "results/data/gencode.vM25.annotation.gtf.gz"
REPEATMASKER = "results/data/mm10_repeats.bed"
ENHANCERS    = "results/data/mm10_enhancers_Encode_Fantom5.bed"

TFBS_BACKGROUND = ["genome"]


# =============================================================================
# Main output list
# =============================================================================

rule_Result3_input_list = [
    f"{OUTPUT_TABLES}/ATAC_DAR.tsv",
    f"{OUTPUT_TABLES}/ATAC_DAR.bed",
    f"{PROCESSED_DATA}/ATAC_sig_DARs.Robj",

    f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
    f"{OUTPUT_PDF}/ATAC_DAR_heatmap.pdf",
    f"{OUTPUT_PNG}/ATAC_DAR_heatmap.png",

    f"{OUTPUT_PDF}/ATAC_DAR_enriched_heatmap.pdf",
    f"{OUTPUT_PNG}/ATAC_DAR_enriched_heatmap.png",

    f"{OUTPUT_PDF}/ATAC_DAR_feature_overlap.pdf",
    f"{OUTPUT_PNG}/ATAC_DAR_feature_overlap.png",

    f"{OUTPUT_PDF}/DAR_GO.pdf",
    f"{OUTPUT_PNG}/DAR_GO.png",

    f"{OUTPUT_TABLES}/ATAC_DAR_TFBS_genome.csv",
    f"{OUTPUT_PDF}/ATAC_DAR_TF_motifs_genome_bg.pdf",
    f"{OUTPUT_PNG}/ATAC_DAR_TF_motifs_genome_bg.png",

    f"{OUTPUT_TABLES}/DAR_TRIM28",
    f"{OUTPUT_PDF}/DAR_TRIM28_overlap.pdf",
    f"{OUTPUT_PNG}/DAR_TRIM28_overlap.png",

    f"{OUTPUT_TABLES}/DAR_TRIM28_TFBS.csv",
    f"{OUTPUT_PDF}/DAR_TRIM28_TFBS.pdf",
    f"{OUTPUT_PNG}/DAR_TRIM28_TFBS.png",

    f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
    f"{OUTPUT_PDF}/DAR_enriched_heatmap.pdf",
    f"{OUTPUT_PNG}/DAR_enriched_heatmap.png",

    f"{OUTPUT_PDF}/DAR_histone_feature_overlap.pdf",
    f"{OUTPUT_PNG}/DAR_histone_feature_overlap.png",

    f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC.csv",
    f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_overlap.pdf",
    f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_overlap.png",

    f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC/DAR_Sertoli_ATAC_TFBS.csv",
    f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_TFBS.pdf",
    f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_TFBS.png",

    f"{OUTPUT_PDF}/DAR_SUMO_enriched_heatmap.pdf",
    f"{OUTPUT_PNG}/DAR_SUMO_enriched_heatmap.png",

    f"{OUTPUT_PDF}/DAR_TRIM28_overlap_fine_custering.pdf",
    f"{OUTPUT_PNG}/DAR_TRIM28_overlap_fine_custering.png",

    f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_overlap_fine_custering.pdf",
    f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_overlap_fine_custering.png",

    f"{OUTPUT_TABLES}/ATAC_DAR_TFBS_genome_fine_clustering.csv",
    f"{OUTPUT_PDF}/ATAC_DAR_TF_motifs_genome_bg_fine_clustering.pdf",
    f"{OUTPUT_PNG}/ATAC_DAR_TF_motifs_genome_bg_fine_clustering.png",
    f"{OUTPUT_TABLES}/ATAC_ctrl_annotation.tsv",
    f"{OUTPUT_TABLES}/ATAC_cKO_annotation.tsv",

    f"{OUTPUT_PDF}/DAR_sex_bias_gene_TF.pdf",
    f"{OUTPUT_PNG}/DAR_sex_bias_gene_TF.png",
    f"{OUTPUT_TABLES}/ATAC_sex_bias-TF_statistics.csv",
    f"{OUTPUT_TABLES}/ATAC_sex_bias-gene_statistics.csv",

    f"{OUTPUT_TABLES}/ATAC_DAR_DE_TE_clusters.tsv",
    f"{OUTPUT_TABLES}/ATAC_DAR_DE_TE_clusters_summary.tsv"
]

# =============================================================================
# Supplementary data - Differential accessible regions
# =============================================================================

rule Annotate_ATAC_ctrl_peaks:
    input:
        ATAC      = f"{ATAC_TABLES}/peaks/consensus_per_condition/WT_consensus_peaks.txt",
        TRIM28    = config["TRIM28_peaks"],
        FOXL2     = config["FOXL2_Ctrl"],
        H3K9me3   = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        genome    = GENOME
    output:
        table = f"{OUTPUT_TABLES}/ATAC_ctrl_annotation.tsv"
    params:
        promoter = config["CR_promoter_distance"]
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/Result_3/SupData_ATAC_peak_annotation.R"

rule Annotate_ATAC_cKO_peaks:
    input:
        ATAC      = f"{ATAC_TABLES}/peaks/consensus_per_condition/KO_consensus_peaks.txt",
        TRIM28    = config["TRIM28_peaks"],
        FOXL2     = config["FOXL2_Ctrl"],
        H3K9me3   = f"{CR_TABLES}/peaks/consensus_per_condition/KO_H3K9me3_consensus_peaks.bed",
        genome    = GENOME
    output:
        table = f"{OUTPUT_TABLES}/ATAC_cKO_annotation.tsv"
    params:
        promoter = config["CR_promoter_distance"]
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/Result_3/SupData_ATAC_peak_annotation.R"


# =============================================================================
# Supplementary data - Differential accessible regions
# =============================================================================

rule SupData_Get_DARs:
    input:
        counts      = f"{ATAC_TABLES}/ATAC_raw_counts.csv",
        samplesheet = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        genome      = GENOME
    params:
        adjpval  = config["ATAC_adjpval"],
        log2FC   = config["ATAC_log2FC"],
        promoter = config["ATAC_promoter_distance"]
    output:
        tsv          = f"{OUTPUT_TABLES}/ATAC_DAR.tsv",
        sig_DARs_bed = f"{OUTPUT_TABLES}/ATAC_DAR.bed",
        sig_DARs     = f"{PROCESSED_DATA}/ATAC_sig_DARs.Robj"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/SupData_ATAC_DAR.R"


rule SupData_Plot_heatmap_DARs:
    input:
        sig_DARs    = f"{PROCESSED_DATA}/ATAC_sig_DARs.Robj",
        norm_counts = f"{ATAC_TABLES}/ATAC_norm_counts.csv"
    output:
        clusters = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        pdf      = f"{OUTPUT_PDF}/ATAC_DAR_heatmap.pdf",
        png      = f"{OUTPUT_PNG}/ATAC_DAR_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/ATAC_DAR_heatmap.R"


# =============================================================================
# Figure 3A - Enriched heatmap
# =============================================================================

rule Fig3A_Plot_enriched_heatmap_DARs:
    input:
        merged_bigwig_check = f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt",
        samplesheet         = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        clusters            = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv"
    params:
        merged_bigwig_folder = ATAC_MERGED_BIGWIG_FOLDER
    output:
        pdf = f"{OUTPUT_PDF}/ATAC_DAR_enriched_heatmap.pdf",
        png = f"{OUTPUT_PNG}/ATAC_DAR_enriched_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/Fig3A_Plot_DAR_enriched_heatmap.R"


# =============================================================================
# Figure 3B - DAR feature overlap
# =============================================================================

rule Fig3B_DAR_feature_overlap:
    input:
        sig_DARs     = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        repeatMasker = REPEATMASKER,
        enhancers    = ENHANCERS,
        genome       = GENOME
    output:
        pdf = f"{OUTPUT_PDF}/ATAC_DAR_feature_overlap.pdf",
        png = f"{OUTPUT_PNG}/ATAC_DAR_feature_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3B_ATAC_DAR_feature_overlap.R"


# =============================================================================
# Supplementary data - GO terms associated with DARs
# =============================================================================

rule SupData_GO_term_DAR_genes:
    input:
        DAR    = f"{OUTPUT_TABLES}/ATAC_DAR.tsv",
        genome = GENOME
    params:
        adjpval = config["ATAC_adjpval"],
        log2FC  = config["ATAC_log2FC"],
        path    = OUTPUT_TABLES
    output:
        pdf = f"{OUTPUT_PDF}/DAR_GO.pdf",
        png = f"{OUTPUT_PNG}/DAR_GO.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/SupData_ATAC_DAR_GO.R"


# =============================================================================
# Figure 3 - TF motif enrichment in DARs
# =============================================================================

rule Fig3_TFBS_motifs_DAR:
    input:
        samplesheet = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        TF_genes    = config["TF_genes"],
        sig_DARs    = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv"
    params:
        TFBS_log10FDR = config["TFBS_log10FDR"],
        background    = "genome",
        logos         = True
    output:
        res_table    = f"{OUTPUT_TABLES}/ATAC_DAR_TFBS_genome.csv",
        pdf          = f"{OUTPUT_PDF}/ATAC_DAR_TF_motifs_genome_bg.pdf",
        png          = f"{OUTPUT_PNG}/ATAC_DAR_TF_motifs_genome_bg.png"
    threads: 48
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/Fig3_ATAC_DAR_motif_enrich.R"


# =============================================================================
# SupFig - TF motif enrichment in DARs
# =============================================================================

rule SupFig_TFBS_motifs_DAR_fine_clustering:
    input:
        # samplesheet = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        TF_genes    = config["TF_genes"],
        sig_DARs    = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv"
    params:
        TFBS_log10FDR = config["TFBS_log10FDR"],
        background    = "genome",
        logos         = True
    output:
        res_table    = f"{OUTPUT_TABLES}/ATAC_DAR_TFBS_genome_fine_clustering.csv",
        pdf          = f"{OUTPUT_PDF}/ATAC_DAR_TF_motifs_genome_bg_fine_clustering.pdf",
        png          = f"{OUTPUT_PNG}/ATAC_DAR_TF_motifs_genome_bg_fine_clustering.png"
    threads: 48
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/Fig3_ATAC_DAR_motif_enrich_fine_clustering.R"

# =============================================================================
# Figure 3C - DAR / TRIM28 overlap
# =============================================================================

rule Fig3C_DAR_TRIM28_overlap:
    input:
        DAR_peaks    = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28_peaks = config["TRIM28_peaks"],
        genome       = GENOME
    output:
        common_ATAC = directory(f"{OUTPUT_TABLES}/DAR_TRIM28"),
        pdf         = f"{OUTPUT_PDF}/DAR_TRIM28_overlap.pdf",
        png         = f"{OUTPUT_PNG}/DAR_TRIM28_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3C_DAR_TRIM28_overlap.R"

# =============================================================================
# SupFig - DAR / TRIM28 overlap
# =============================================================================

rule SupFig_DAR_TRIM28_overlap:
    input:
        DAR_peaks    = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        TRIM28_peaks = config["TRIM28_peaks"],
        genome       = GENOME
    output:
        common_ATAC = directory(f"{OUTPUT_TABLES}/DAR_TRIM28_fine_clustering"),
        pdf         = f"{OUTPUT_PDF}/DAR_TRIM28_overlap_fine_custering.pdf",
        png         = f"{OUTPUT_PNG}/DAR_TRIM28_overlap_fine_custering.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/SupFig_DAR_TRIM28_overlap_fine_clustering.R"

# =============================================================================
# Figure 3D - TFBS enrichment in DAR / TRIM28 regions
# =============================================================================

rule Fig3D_DAR_TRIM28_TFBS:
    input:
        DAR_peaks    = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28_peaks = config["TRIM28_peaks"]
    output:
        # TFBS_log10FDR = config["TFBS_log10FDR"],
        res_table = f"{OUTPUT_TABLES}/DAR_TRIM28_TFBS.csv",
        pdf       = f"{OUTPUT_PDF}/DAR_TRIM28_TFBS.pdf",
        png       = f"{OUTPUT_PNG}/DAR_TRIM28_TFBS.png"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_3/Fig3D_DAR_TRIM28_TFBS.R"


# =============================================================================
# Figure 3E - Multi-omics enriched heatmap
# =============================================================================

rule Fig3E_Plot_enriched_heatmap_DARs:
    input:
        ATAC_bw_res  = f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt",
        CR_bw_res    = expand(f"{CR_MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt", AB=AB),
        samplesheet  = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        sig_DARs     = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        DAR_anno     = f"{OUTPUT_TABLES}/ATAC_DAR.tsv",
        TRIM28_peaks = config["TRIM28_peaks"]
    params:
        AB      = AB,
        ATAC_bw = ATAC_MERGED_BIGWIG_FOLDER,
        CR_bw   = CR_MERGED_BIGWIG_FOLDER
    output:
        atac_clustering = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        pdf             = f"{OUTPUT_PDF}/DAR_enriched_heatmap.pdf",
        png             = f"{OUTPUT_PNG}/DAR_enriched_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_3/Fig3E_DAR_histones_enriched_heatmap.R"


# =============================================================================
# Figure 3F - Histone feature overlap
# =============================================================================

rule Fig3F_DAR_histone_feature_overlap:
    input:
        atac_clustering = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        repeatMasker    = REPEATMASKER,
        enhancers       = ENHANCERS,
        genome          = GENOME
    output:
        pdf = f"{OUTPUT_PDF}/DAR_histone_feature_overlap.pdf",
        png = f"{OUTPUT_PNG}/DAR_histone_feature_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3F_DAR_histone_feature_overlap.R"


# =============================================================================
# Figure 3G - DAR / Sertoli-specific ATAC overlap
# =============================================================================

rule Fig3G_DAR_Sertoli_ATAC_overlap:
    input:
        ATAC_peaks         = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        ATAC_Sertoli_peaks = config["Sertoli_spe_DAR"],
        ATAC_Granulosa_peaks = config["Granulosa_spe_DAR"],
        genome             = GENOME
    output:
        common_ATAC = f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC.csv",
        pdf         = f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_overlap.pdf",
        png         = f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3G_DAR_Sertoli_ATAC_overlap.R"

# =============================================================================
# SupFig - DAR / Sertoli-specific ATAC overlap
# =============================================================================

rule SupFig_DAR_Sertoli_ATAC_overlap:
    input:
        ATAC_peaks         = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        ATAC_Sertoli_peaks = config["Sertoli_spe_DAR"],
        ATAC_Granulosa_peaks = config["Granulosa_spe_DAR"],
        genome             = GENOME
    output:
        common_ATAC = f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC_fine_custering.csv",
        pdf         = f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_overlap_fine_custering.pdf",
        png         = f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_overlap_fine_custering.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3G_DAR_Sertoli_ATAC_overlap_fine_clustering.R"

# =============================================================================
# Supplementary figure - TFBS enrichment in DAR / Sertoli-specific ATAC regions
# =============================================================================

rule SupFig_DAR_Sertoli_ATAC_TFBS:
    input:
        sig_DARs = f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC.csv"
    params:
        TFBS_log10FDR = config["TFBS_log10FDR"]
    output:
        res_table = f"{OUTPUT_TABLES}/DAR_Sertoli_ATAC/DAR_Sertoli_ATAC_TFBS.csv",
        pdf       = f"{OUTPUT_PDF}/DAR_Sertoli_ATAC_TFBS.pdf",
        png       = f"{OUTPUT_PNG}/DAR_Sertoli_ATAC_TFBS.png"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_3/SupFig_DAR_Sertoli_ATAC_TFBS.R"


# =============================================================================
# Figure 3H - ATAC/SUMO enriched heatmap
# =============================================================================

rule Fig3H_DAR_SUMO_enriched_heatmap:
    input:
        TRIM28_peaks = config["TRIM28_peaks"],
        atac_clustering = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        ATAC_bw_res     = f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt",
        ChIP_bw_res     = expand(f"{CHIP_MERGED_BIGWIG_FOLDER}/{{SUMO}}_merge_bw.txt", SUMO=SUMO),
        samplesheet     = f"{ATAC_TABLES}/ATAC_samplesheet.csv",
        sig_DARs        = f"{OUTPUT_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        DAR_anno        = f"{OUTPUT_TABLES}/ATAC_DAR.tsv",
    params:
        ATAC_bw = ATAC_MERGED_BIGWIG_FOLDER,
        ChIP_bw = CHIP_MERGED_BIGWIG_FOLDER,
        SUMO = SUMO
    output:
        pdf = f"{OUTPUT_PDF}/DAR_SUMO_enriched_heatmap.pdf",
        png = f"{OUTPUT_PNG}/DAR_SUMO_enriched_heatmap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_3/Fig3H_DAR_SUMO_enriched_heatmap.R"


rule ATAC_DAR_overlap_summary:
    input:
        ATAC = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        H3K9me3_ctrl   = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        H3K9me3_cKO   = f"{CR_TABLES}/peaks/consensus_per_condition/KO_H3K9me3_consensus_peaks.bed",
        TRIM28        = config["TRIM28_peaks"],
        FOXL2         = config["FOXL2_Ctrl"],
        SOX9          = config["SOX9"],
        DMRT1          = config["DMRT1"],
        ATAC_Sertoli   = config["Sertoli_spe_DAR"],
        ATAC_Granulosa = config["Granulosa_spe_DAR"],
        genome        = "results/data/gencode.vM25.annotation.gtf.gz",

        expressed_genes = config["expressed_genes"],
        DEG_8weeks = config["scRNAseq"],
        DEG_7months = config["bulk_RNAseq"],
        sex_bias = config["sex_biased_genes"]
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        table = f"{OUTPUT_TABLES}/ATAC_DAR_overlap_summary.tsv"
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/Result_3/SupData_DAR_summary_table.R"




rule ATAC_DAR_sex_bias_gene_TF:
    input:
        ATAC_summary = f"{OUTPUT_TABLES}/ATAC_DAR_overlap_summary.tsv",
    output:
        pdf = f"{OUTPUT_PDF}/DAR_sex_bias_gene_TF.pdf",
        png = f"{OUTPUT_PNG}/DAR_sex_bias_gene_TF.png",
        TF_statistics = f"{OUTPUT_TABLES}/ATAC_sex_bias-TF_statistics.csv",
        gene_statistics = f"{OUTPUT_TABLES}/ATAC_sex_bias-gene_statistics.csv",
        expression_statistics = f"{OUTPUT_TABLES}/ATAC_sex_bias-expression_statistics.csv"
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/Result_3/Fig3E_DAR_sex_bias_gene_TF.R"


rule ATAC_DAR_DE_TE_clusters:
    input:
        TE_up = config["Up_TE_cKO"],
        TE_down = config["Down_TE_cKO"],
        ATAC_clusters = f"{OUTPUT_TABLES}/DAR_enriched_heatmap_clustering.csv",
        repeatMasker = REPEATMASKER
    output:
        table = f"{OUTPUT_TABLES}/ATAC_DAR_DE_TE_clusters.tsv",
        summary = f"{OUTPUT_TABLES}/ATAC_DAR_DE_TE_clusters_summary.tsv"
    threads: 1
    resources:
        mem_mb = 8000
    script:
        "../scripts/Result_3/SupData_DAR_TE_clusters.R"