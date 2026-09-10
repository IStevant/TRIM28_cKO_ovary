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

CR_PROCESSED_DATA = "results/CR/processed_data"
CR_TABLES         = "results/CR/tables"

OUTPUT_TABLES = "results/Result_1/tables"
OUTPUT_PNG    = "results/Result_1/graphs/PNG"
OUTPUT_PDF    = "results/Result_1/graphs/PDF"
TMP           = "results/Result_1/tmp"

CR_MERGED_BIGWIG_FOLDER = f"{CR_PROCESSED_DATA}/CR_bigwig/merged"

TFBS_BACKGROUND = ["genome", "conditions"]


# =============================================================================
# Main output list
# =============================================================================

rule_Result1_input_list = [
    f"{OUTPUT_TABLES}/peaks/H3K9me3_TRIM28_overlap.txt",
    f"{OUTPUT_TABLES}/peaks/TRIM28_H3K9me3_overlap.txt",
    f"{OUTPUT_PDF}/Fig1E_H3K9me3_TRIM28_overlap.pdf",
    f"{OUTPUT_PNG}/Fig1E_H3K9me3_TRIM28_overlap.png",

    f"{OUTPUT_PDF}/Fig1F_genomic_track_example.pdf",
    f"{OUTPUT_PNG}/Fig1F_genomic_track_example.png",

    f"{OUTPUT_PDF}/FigS1_H3K9me3_TRIM28_FOXL2_overlap.pdf",
    f"{OUTPUT_PNG}/FigS1_H3K9me3_TRIM28_FOXL2_overlap.png",
    f"{OUTPUT_TABLES}/SupData1_H3K9me3_TRIM28_FOXL2_overlap.csv",

    f"{OUTPUT_PDF}/Fig1G_H3K9me3_TRIM28_feature_overlap.pdf",
    f"{OUTPUT_PNG}/Fig1G_H3K9me3_TRIM28_feature_overlap.png",

    f"{OUTPUT_TABLES}/TFBS/Multi_TRIM28_only_TFBS.csv",
    f"{OUTPUT_TABLES}/TFBS/Multi_TRIM28_H3K9me3_TFBS.csv",
    f"{OUTPUT_PDF}/Fig1H_H3K9me3_TRIM28_TFBS.pdf",
    f"{OUTPUT_PNG}/Fig1H_H3K9me3_TRIM28_TFBS.png",

    f"{OUTPUT_TABLES}/TRIM28_peak_annotation.csv",
    f"{OUTPUT_TABLES}/FOXL2_peak_annotation.csv",
    f"{OUTPUT_TABLES}/ESR2_peak_annotation.csv",
    f"{OUTPUT_TABLES}/NR5A2_peak_annotation.csv",
    f"{OUTPUT_TABLES}/RUNX1_peak_annotation.csv",
    f"{OUTPUT_TABLES}/SOX9_peak_annotation.csv",
    f"{OUTPUT_TABLES}/DMRT1_peak_annotation.csv"
]


# =============================================================================
# Figure 1E - H3K9me3 / TRIM28 overlap
# =============================================================================

rule Fig1E_H3K9me3_vs_TRIM28_overlap:
    input:
        H3K9me3_peaks = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        TRIM28_peaks  = config["TRIM28_peaks"]
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"]
    output:
        common_h3k9me3 = f"{OUTPUT_TABLES}/peaks/H3K9me3_TRIM28_overlap.txt",
        common_trim28  = f"{OUTPUT_TABLES}/peaks/TRIM28_H3K9me3_overlap.txt",
        pdf            = f"{OUTPUT_PDF}/Fig1E_H3K9me3_TRIM28_overlap.pdf",
        png            = f"{OUTPUT_PNG}/Fig1E_H3K9me3_TRIM28_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_1/FigE_H3K9me3_TRIM28_overlap_venn.R"

# =============================================================================
# Figure 1F - Genomic tracks
# =============================================================================

rule Fig1F_plot_genomic_tracks_example:
    input:
        H3K9me3_bigwig  = f"{CR_MERGED_BIGWIG_FOLDER}/WT_H3K9me3.bw",
        H3K9me3_domains = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_domains.bed",
        TRIM28_peaks    = config["TRIM28_peaks"],
        FOXL2_peaks     = config["FOXL2_Ctrl"],
        repeatMasker    = "results/data/mm10_repeats.bed",
        CRE             = "results/data/mm10_enhancers_Encode_Fantom5.bed",
        genome          = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        region = "chr9:22186341-22299046"
    output:
        png = f"{OUTPUT_PNG}/Fig1F_genomic_track_example.png",
        pdf = f"{OUTPUT_PDF}/Fig1F_genomic_track_example.pdf"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_1/Fig1F_Plot_genomic_tracks_example.R"

# =============================================================================
# Figure S1 / Supplementary Data 2 - H3K9me3, TRIM28 and FOXL2 overlap
# =============================================================================

rule FigS1_SupData1_H3K9me3_vs_TRIM28_vs_FOXL2:
    input:
        H3K9me3 = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        TRIM28  = config["TRIM28_peaks"],
        FOXL2   = config["FOXL2_Ctrl"],
        genome  = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        pdf   = f"{OUTPUT_PDF}/FigS1_H3K9me3_TRIM28_FOXL2_overlap.pdf",
        png   = f"{OUTPUT_PNG}/FigS1_H3K9me3_TRIM28_FOXL2_overlap.png",
        table = f"{OUTPUT_TABLES}/SupData1_H3K9me3_TRIM28_FOXL2_overlap.csv"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_1/FigS1_SupData1_H3K9me3_TRIM28_FOXL2_overlap.R"


# =============================================================================
# Figure 1G - Genomic feature overlap
# =============================================================================

rule Fig1G_H3K9me3_vs_TRIM28_feature_overlap:
    input:
        H3K9me3_peaks = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        TRIM28_peaks  = config["TRIM28_peaks"],
        repeatMasker  = "results/data/mm10_repeats.bed",
        enhancers     = "results/data/mm10_enhancers_Encode_Fantom5.bed",
        genome        = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        distance_to_H3K9me3 = config["distance_to_H3K9me3"],
        promoter            = config["CR_promoter_distance"]
    output:
        pdf = f"{OUTPUT_PDF}/Fig1G_H3K9me3_TRIM28_feature_overlap.pdf",
        png = f"{OUTPUT_PNG}/Fig1G_H3K9me3_TRIM28_feature_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_1/Fig1G_H3K9me3_TRIM28_feature_overlap.R"


# =============================================================================
# Figure 1H - TFBS enrichment
# =============================================================================

rule Fig1H_H3K9me3_vs_TRIM28_TFBS:
    input:
        H3K9me3_peaks = f"{CR_TABLES}/peaks/consensus_per_condition/WT_H3K9me3_consensus_peaks.bed",
        TRIM28_peaks  = config["TRIM28_peaks"]
    params:
        TFBS_log10FDR = config["TFBS_log10FDR"]
    output:
        res_table_TRIM28_only    = f"{OUTPUT_TABLES}/TFBS/Multi_TRIM28_only_TFBS.csv",
        res_table_TRIM28_H3K9me3 = f"{OUTPUT_TABLES}/TFBS/Multi_TRIM28_H3K9me3_TFBS.csv",
        res_table_H3K9me3_only    = f"{OUTPUT_TABLES}/TFBS/Multi_H3K9me3_only_TFBS.csv",
        pdf                      = f"{OUTPUT_PDF}/Fig1H_H3K9me3_TRIM28_TFBS.pdf",
        png                      = f"{OUTPUT_PNG}/Fig1H_H3K9me3_TRIM28_TFBS.png"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/Fig1H_TRIM28_H3K9me3_TFBS.R"


# =============================================================================
# SupData2 - Peak annotation
# =============================================================================

rule SupData2_annotate_TRIM28_peaks:
    input:
        bed    = config["TRIM28_peaks"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/TRIM28_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_FOXL2_peaks:
    input:
        bed    = config["FOXL2_Ctrl"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/FOXL2_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_ESR2_peaks:
    input:
        bed    = config["ESR2"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/ESR2_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_NR5A2_peaks:
    input:
        bed    = config["NR5A2"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/NR5A2_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_RUNX1_peaks:
    input:
        bed    = config["RUNX"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/RUNX1_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_SOX9_peaks:
    input:
        bed    = config["SOX9"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/SOX9_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"

rule SupData2_annotate_DMRT1_peaks:
    input:
        bed    = config["DMRT1"],
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        promoter = config["distance_to_TSS"]
    output:
        table    = f"{OUTPUT_TABLES}/DMRT1_peak_annotation.csv"
    threads: 12
    resources:
        mem_mb = 98000
    script:
        "../scripts/Result_1/SupData2_annotate_peaks.R"