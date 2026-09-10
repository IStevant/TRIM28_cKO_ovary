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

CR_PROCESSED_DATA   = "results/CR/processed_data"
CR_TABLES           = "results/CR/tables"

CHIP_PROCESSED_DATA = "results/ChIP/processed_data"

OUTPUT_TABLES  = "results/Result_2/tables"
PROCESSED_DATA = "results/Result_2/processed_data"
OUTPUT_PNG     = "results/Result_2/graphs/PNG"
OUTPUT_PDF     = "results/Result_2/graphs/PDF"
TMP            = "results/Result_2/tmp"

AB   = ["H3K9me3"]
SUMO = ["SUMO1", "SUMO2"]

CR_NORM_BIGWIG_FOLDER   = f"{CR_PROCESSED_DATA}/CR_bigwig"
CR_MERGED_BIGWIG_FOLDER = f"{CR_NORM_BIGWIG_FOLDER}/merged"

CHIP_NORM_BIGWIG_FOLDER   = f"{CHIP_PROCESSED_DATA}/ChIP_bigwig"
CHIP_MERGED_BIGWIG_FOLDER = f"{CHIP_NORM_BIGWIG_FOLDER}/merged"

GENOME       = "results/data/gencode.vM25.annotation.gtf.gz"
REPEATMASKER = "results/data/mm10_repeats.bed"
ENHANCERS    = "results/data/mm10_enhancers_Encode_Fantom5.bed"

TFBS_BACKGROUND = ["genome"]

REGIONS = {
    "region1": "chr11:21712388-21978474",
    "region2": "chr4:113211887-113246732",
}

wildcard_constraints:
    AB = "H3K9me3"


# =============================================================================
# Main output list
# =============================================================================

rule_Result2_input_list = [
    expand(f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER.tsv", AB=AB),
    expand(f"{PROCESSED_DATA}/Diff_Analysis/{{AB}}_sig_DERs.Robj", AB=AB),
    expand(f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_sig_DERs.bed", AB=AB),
    expand(f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_heatmap_clusters.csv", AB=AB),
    expand(f"{OUTPUT_PDF}/{{AB}}_DER_heatmap.pdf", AB=AB),
    expand(f"{OUTPUT_PNG}/{{AB}}_DER_heatmap.png", AB=AB),
    expand(f"{OUTPUT_PDF}/{{AB}}_DER_enriched_heatmap.pdf", AB=AB),
    expand(f"{OUTPUT_PNG}/{{AB}}_DER_enriched_heatmap.png", AB=AB),
    expand(f"{OUTPUT_PNG}/Fig2A_genomic_track_example_{{region}}.png", region=REGIONS.keys()),
    expand(f"{OUTPUT_PDF}/Fig2A_genomic_track_example_{{region}}.pdf", region=REGIONS.keys()),
    expand(f"{OUTPUT_PDF}/Fig2B_DER_{{AB}}_TRIM28_overlap.pdf", AB=AB),
    expand(f"{OUTPUT_PNG}/Fig2B_DER_{{AB}}_TRIM28_overlap.png", AB=AB),
    f"{OUTPUT_PDF}/Fig2C_DER_H3K9me3_TRIM28_feature_overlap.pdf",
    f"{OUTPUT_PNG}/Fig2C_DER_H3K9me3_TRIM28_feature_overlap.png",
    expand(f"{OUTPUT_PDF}/{{AB}}_DER_GO.pdf", AB=AB),
    expand(f"{OUTPUT_PNG}/{{AB}}_DER_GO.png", AB=AB),
    f"{OUTPUT_PDF}/Fig2D_H3K9me3_DER_TF_motifs_genome_bg.pdf",
    f"{OUTPUT_PNG}/Fig2D_H3K9me3_DER_TF_motifs_genome_bg.png",
    f"{OUTPUT_TABLES}/TFBS_Enrichment/H3K9me3_DER_TF_motifs_genome_bg.csv",
    f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_WT_TE_enrichment.csv",
    f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_loss_TE_enrichment.csv",
    f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_gain_TE_enrichment.csv",
    f"{OUTPUT_PDF}/H3K9me3_DER_TE.pdf",
    f"{OUTPUT_PNG}/H3K9me3_DER_TE.png",
    f"{OUTPUT_TABLES}/DER_H3K9me3_Sertoli_ATAC_overlap",
    f"{OUTPUT_PDF}/SupFig_DER_H3K9me3_Sertoli_ATAC_overlap.pdf",
    f"{OUTPUT_PNG}/SupFig_DER_H3K9me3_Sertoli_ATAC_overlap.png",
    expand(f"{OUTPUT_PDF}/{{AB}}_DER_SUMO_enriched_heatmap_TRIM28_split.pdf", AB=AB),
    expand(f"{OUTPUT_PNG}/{{AB}}_DER_SUMO_enriched_heatmap_TRIM28_split.png", AB=AB),
    f"{OUTPUT_PDF}/FigS3_DER_H3K9me3_TRIM28_FOXL2_overlap.pdf",
    f"{OUTPUT_PNG}/FigS3_DER_H3K9me3_TRIM28_FOXL2_overlap.png",
    f"{OUTPUT_TABLES}/SupData_DER_H3K9me3_TRIM28_FOXL2_overlap.csv",
    f"{OUTPUT_PDF}/SupFig_H3K9me3_TRIM28_FOXL2_overlap.pdf",
    f"{OUTPUT_PNG}/SupFig_H3K9me3_TRIM28_FOXL2_overlap.png",
    f"{OUTPUT_TABLES}/SupData1_H3K9me3_TRIM28_FOXL2_overlap.csv",
    f"{OUTPUT_TABLES}/H3K9me3_DER_overlap_summary.tsv",
    expand(f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_gene_expression.tsv", AB=AB)
]


# =============================================================================
# Supplementary data  4- Differential enriched regions
# =============================================================================

rule SupData4_Get_DERs:
    input:
        counts      = f"{CR_TABLES}/matrices/{{AB}}_raw_counts.csv",
        samplesheet = f"{CR_TABLES}/matrices/{{AB}}_samplesheet.csv",
        genome      = GENOME
    params:
        adjpval  = config["CR_adjpval"],
        log2FC   = config["CR_log2FC"],
        promoter = config["CR_promoter_distance"]
    output:
        tsv      = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER.tsv",
        sig_DERs = f"{PROCESSED_DATA}/Diff_Analysis/{{AB}}_sig_DERs.Robj",
        DER_bed  = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_sig_DERs.bed"
    threads: 12
    resources:
        mem_mb = 84000
    script:
        "../scripts/Result_2/SupData4_CUTandRUN_DER.R"


rule Plot_heatmap_DERs:
    input:
        sig_DERs    = f"{PROCESSED_DATA}/Diff_Analysis/{{AB}}_sig_DERs.Robj",
        norm_counts = f"{CR_TABLES}/matrices/{{AB}}_norm_counts.csv"
    output:
        clusters = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_heatmap_clusters.csv",
        pdf      = f"{OUTPUT_PDF}/{{AB}}_DER_heatmap.pdf",
        png      = f"{OUTPUT_PNG}/{{AB}}_DER_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/CR_DER_heatmap.R"


# =============================================================================
# Figure 2A - Enriched heatmap
# =============================================================================

rule Fig2A_Plot_enriched_heatmap_DERs:
    input:
        merged_bigwig_check = f"{CR_MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt",
        samplesheet         = f"{CR_TABLES}/matrices/{{AB}}_samplesheet.csv",
        clusters            = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_heatmap_clusters.csv"
    params:
        merged_bigwig_folder = CR_MERGED_BIGWIG_FOLDER
    output:
        pdf = f"{OUTPUT_PDF}/{{AB}}_DER_enriched_heatmap.pdf",
        png = f"{OUTPUT_PNG}/{{AB}}_DER_enriched_heatmap.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/Fig2A_Plot_DER_enriched_heatmap.R"


rule Fig2A_plot_genomic_tracks_example:
    input:
        H3K9me3_WT_bigwig  = f"{CR_MERGED_BIGWIG_FOLDER}/WT_H3K9me3.bw",
        H3K9me3_KO_bigwig  = f"{CR_MERGED_BIGWIG_FOLDER}/KO_H3K9me3.bw",
        H3K9me3_WT_domains = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_domains.bed",
        H3K9me3_KO_domains = f"{CR_TABLES}/peaks/consensus_per_condition/KO_H3K9me3_consensus_domains.bed",
        TRIM28_peaks       = config["TRIM28_peaks"],
        FOXL2_peaks        = config["FOXL2_Ctrl"],
        repeatMasker       = REPEATMASKER,
        CRE                = ENHANCERS,
        genome             = GENOME
    params:
        region = lambda wildcards: REGIONS[wildcards.region]
    output:
        png = f"{OUTPUT_PNG}/Fig2A_genomic_track_example_{{region}}.png",
        pdf = f"{OUTPUT_PDF}/Fig2A_genomic_track_example_{{region}}.pdf"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/Fig2A_Plot_genomic_tracks_example.R"


# =============================================================================
# Figure 2B - DER H3K9me3 / TRIM28 overlap
# =============================================================================

rule Fig2B_DER_H3K9me3_TRIM28_overlap:
    input:
        H3K9me3_DER  = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_heatmap_clusters.csv",
        TRIM28_peaks = config["TRIM28_peaks"],
        genome       = GENOME
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        common_H3K9me3 = directory(f"{OUTPUT_TABLES}/DER_{{AB}}_TRIM28"),
        pdf            = f"{OUTPUT_PDF}/Fig2B_DER_{{AB}}_TRIM28_overlap.pdf",
        png            = f"{OUTPUT_PNG}/Fig2B_DER_{{AB}}_TRIM28_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_2/Fig2B_DER_H3K9me3_TRIM28_overlap.R"


# =============================================================================
# Figure 2C - DER H3K9me3 / TRIM28 genomic feature overlap
# =============================================================================

rule Fig2C_DER_H3K9me3_TRIM28_feature_overlap:
    input:
        H3K9me3_DER  = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv",
        TRIM28_peaks = config["TRIM28_peaks"],
        repeatMasker = REPEATMASKER,
        enhancers    = ENHANCERS,
        genome       = GENOME
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"]
    output:
        pdf = f"{OUTPUT_PDF}/Fig2C_DER_H3K9me3_TRIM28_feature_overlap.pdf",
        png = f"{OUTPUT_PNG}/Fig2C_DER_H3K9me3_TRIM28_feature_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_2/Fig2C_DER_H3K9me3_TRIM28_feature_overlap.R"


# =============================================================================
# Supplementary data - GO terms associated with H3K9me3 DERs
# =============================================================================

rule SupData_GO_term_DER_genes:
    input:
        DER    = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER.tsv",
        genome = GENOME
    params:
        AB      = "{AB}",
        adjpval = config["CR_adjpval"],
        log2FC  = config["CR_log2FC"],
        path    = f"{OUTPUT_TABLES}/Diff_Analysis"
    output:
        pdf = f"{OUTPUT_PDF}/{{AB}}_DER_GO.pdf",
        png = f"{OUTPUT_PNG}/{{AB}}_DER_GO.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/CR_DER_GO.R"


# =============================================================================
# Supplementary data - Gene expression in genes with H3K9me3 DERs
# =============================================================================

rule SupData_Expression_DER_genes:
    input:
        DER    = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER.tsv",
        bulk_DEG = config["bulk_RNAseq"],
        scRNA_DEG = config["scRNAseq"],
        genome = GENOME
    params:
        AB      = "{AB}",
        adjpval = config["CR_adjpval"],
        log2FC  = config["CR_log2FC"],
    output:
        table = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_gene_expression.tsv"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/SupData_DER_genes.R"


# =============================================================================
# Figure 2D - TF motif enrichment in H3K9me3 DERs
# =============================================================================

rule Fig2D_TFBS_motifs_DER:
    input:
        samplesheet = f"{CR_TABLES}/matrices/H3K9me3_samplesheet.csv",
        TF_genes    = config["TF_genes"],
        sig_DERs    = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv"
    params:
        TFBS_log10FDR = config["TFBS_log10FDR"],
        background    = "genome",
        logos         = True
    output:
        pdf       = f"{OUTPUT_PDF}/Fig2D_H3K9me3_DER_TF_motifs_genome_bg.pdf",
        png       = f"{OUTPUT_PNG}/Fig2D_H3K9me3_DER_TF_motifs_genome_bg.png",
        res_table = f"{OUTPUT_TABLES}/TFBS_Enrichment/H3K9me3_DER_TF_motifs_genome_bg.csv"
    threads: 48
    resources:
        mem_mb = 95000
    script:
        "../scripts/Result_2/Fig2D_DER_TFBS.R"


# =============================================================================
# Supplementary figure - TE enrichment in H3K9me3 DERs
# =============================================================================

rule SupFig_TE_DER:
    input:
        WT_peaks     = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        DER          = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv",
        repeatMasker = REPEATMASKER
    output:
        WT_table     = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_WT_TE_enrichment.csv",
        lost_table   = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_loss_TE_enrichment.csv",
        gained_table = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_gain_TE_enrichment.csv",
        pdf          = f"{OUTPUT_PDF}/H3K9me3_DER_TE.pdf",
        png          = f"{OUTPUT_PNG}/H3K9me3_DER_TE.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/H3K9me3_DER_TE_overlap.R"


# =============================================================================
# Supplementary figure - DER H3K9me3 / Sertoli-specific ATAC overlap
# =============================================================================

rule SupFig_DER_H3K9me3_Sertoli_ATAC_overlap:
    input:
        H3K9me3_DER          = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv",
        ATAC_Sertoli_peaks   = config["Sertoli_spe_DAR"],
        ATAC_Granulosa_peaks = config["Granulosa_spe_DAR"],
        genome               = GENOME
    output:
        common_H3K9me3 = directory(f"{OUTPUT_TABLES}/DER_H3K9me3_Sertoli_ATAC_overlap"),
        pdf            = f"{OUTPUT_PDF}/SupFig_DER_H3K9me3_Sertoli_ATAC_overlap.pdf",
        png            = f"{OUTPUT_PNG}/SupFig_DER_H3K9me3_Sertoli_ATAC_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_2/SupFig_DER_H3K9me3_Sertoli_ATAC_overlap.R"


# =============================================================================
# Figure 2H - H3K9me3 DER / SUMO enriched heatmap
# =============================================================================

rule Fig2_DER_SUMO_enriched_heatmap:
    input:
        merged_bigwig_check = f"{CR_MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt",
        samplesheet         = f"{CR_TABLES}/matrices/{{AB}}_samplesheet.csv",
        clusters            = f"{OUTPUT_TABLES}/Diff_Analysis/{{AB}}_DER_heatmap_clusters.csv",
        TRIM28_peaks        = config["TRIM28_peaks"]
    params:
        CR_bw   = CR_MERGED_BIGWIG_FOLDER,
        ChIP_bw = CHIP_MERGED_BIGWIG_FOLDER,
        SUMO    = SUMO
    output:
        pdf = f"{OUTPUT_PDF}/{{AB}}_DER_SUMO_enriched_heatmap_TRIM28_split.pdf",
        png = f"{OUTPUT_PNG}/{{AB}}_DER_SUMO_enriched_heatmap_TRIM28_split.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_2/Fig2_DER_SUMO_enriched_heatmap.R"

# =============================================================================
# Figure S3 / Supplementary Data  - DER H3K9me3, TRIM28 and FOXL2 overlap
# =============================================================================

rule FigS3_SupData_DER_H3K9me3_vs_TRIM28_vs_FOXL2:
    input:
        H3K9me3_DER = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv",
        TRIM28      = config["TRIM28_peaks"],
        FOXL2       = config["FOXL2_Ctrl"],
        genome      = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        pdf   = f"{OUTPUT_PDF}/FigS3_DER_H3K9me3_TRIM28_FOXL2_overlap.pdf",
        png   = f"{OUTPUT_PNG}/FigS3_DER_H3K9me3_TRIM28_FOXL2_overlap.png",
        table = f"{OUTPUT_TABLES}/SupData_DER_H3K9me3_TRIM28_FOXL2_overlap.csv"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_2/FigS3_DER_H3K9me3_TRIM28_FOXL2_overlap.R"

# =============================================================================
# Supplementary Data  - H3K9me3, TRIM28 and FOXL2 overlap
# =============================================================================

rule Supdata_H3K9me3_vs_TRIM28_vs_FOXL2:
    input:
        H3K9me3 = f"{CR_TABLES}/peaks/consensus_per_condition/KO_H3K9me3_consensus_peaks.bed",
        TRIM28  = config["TRIM28_peaks"],
        FOXL2   = config["FOXL2_Ctrl"],
        genome  = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        pdf   = f"{OUTPUT_PDF}/SupFig_H3K9me3_TRIM28_FOXL2_overlap.pdf",
        png   = f"{OUTPUT_PNG}/SupFig_H3K9me3_TRIM28_FOXL2_overlap.png",
        table = f"{OUTPUT_TABLES}/SupData1_H3K9me3_TRIM28_FOXL2_overlap.csv"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_1/FigS1_SupData1_H3K9me3_TRIM28_FOXL2_overlap.R"



rule H3K9me3_DER_overlap_summary:
    input:
        H3K9me3_DER   = f"{OUTPUT_TABLES}/Diff_Analysis/H3K9me3_DER_heatmap_clusters.csv",
        TRIM28        = config["TRIM28_peaks"],
        FOXL2         = config["FOXL2_Ctrl"],
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
        table = f"{OUTPUT_TABLES}/H3K9me3_DER_overlap_summary.tsv"
    threads: 1

    resources:
        mem_mb = 8000

    script:
        "../scripts/Result_2/SupData_DER_summary_table.R"