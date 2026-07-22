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

CR_PROCESSED_DATA = "results/CR/processed_data"
CR_TABLES         = "results/CR/tables"

CHIP_PROCESSED_DATA = "results/ChIP/processed_data"

OUTPUT_TABLES  = "results/Result_4/tables"
PROCESSED_DATA = "results/Result_4/processed_data"
OUTPUT_PNG     = "results/Result_4/graphs/PNG"
OUTPUT_PDF     = "results/Result_4/graphs/PDF"
TMP            = "results/Result_4/tmp"

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

GENOMIC_TRACK_REGIONS = {
    "Foxl2": "chr9:98952876-98967592",
    "Esr2": "chr12:76129810-76137448",
    "Nr5a2": "chr1:137026429-137037762",
    "Gata4": "chr14:63242905-63251802",
    "Fog2-Zfpm2": "chr15:40759587-40769891",
    "Bmp2": "chr2:133587988-133601452",
    "Wt1": "chr2:105120311-105139938",
    "Fst": "chr13:114755817-114767576",
    "Emx2": "chr19:59453588-59467927",
    "Wnt4": "chr4:137273516-137293930",
    "Inha": "chr1:75507014-75521222",
}


# =============================================================================
# Main output list
# =============================================================================

rule_Result4_input_list = [
    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/combinaison_TF_overlaps.csv",
    f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap.pdf",
    f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap.png",

    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/combinaison_TF_overlaps_genome_bg.csv",
    f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap_genome_bg.pdf",
    f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap_genome_bg.png",

    directory(f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/pie_chart"),
    f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap_pie_chart.pdf",
    f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap_pie_chart.png",

    f"{OUTPUT_TABLES}/TF_overlaps/TF_overlap_regions_heatmap.csv",
    f"{OUTPUT_PDF}/Fig4A_TF_colocalisation_heatmap.pdf",
    f"{OUTPUT_PNG}/Fig4A_TF_colocalisation_heatmap.png",

    f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/TF_overlap_TRIM28.csv",
    f"{OUTPUT_PDF}/Fig4A_TF_colocalisation_upset.pdf",
    f"{OUTPUT_PNG}/Fig4A_TF_colocalisation_upset.png",

    f"{OUTPUT_TABLES}/TF_overlaps/TF_overlap_regions_upset.csv",
    f"{OUTPUT_PDF}/Fig4B_TF_colocalisation_TRIM28.pdf",
    f"{OUTPUT_PNG}/Fig4B_TF_colocalisation_TRIM28.png",

    expand(f"{OUTPUT_PDF}/ATAC_TRIM28_TF_tracks/{{region_id}}.pdf", region_id=GENOMIC_TRACK_REGIONS.keys()),
    expand(f"{OUTPUT_PNG}/ATAC_TRIM28_TF_tracks/{{region_id}}.png", region_id=GENOMIC_TRACK_REGIONS.keys()),
]


# =============================================================================
# Shared TF inputs
# =============================================================================

TF_INPUTS = {
    "FOXL2": config["FOXL2_Ctrl"],
    "NR5A2": config["NR5A2"],
    "ESR2": config["ESR2"],
    "RUNX": config["RUNX"],
}


# =============================================================================
# Figure 4A - Ovarian TF overlap with TRIM28
# =============================================================================

rule Fig4A_Ovarian_TFs_overlap_TRIM28_ATAC_bg:
    input:
        ATAC_control = f"{ATAC_TABLES}/peaks/consensus_per_condition/WT_consensus_peaks.txt",
        TRIM28       = config["TRIM28_peaks"],
        genome       = GENOME,
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/combinaison_TF_overlaps.csv",
        pdf   = f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap.pdf",
        png   = f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4A_TF_TRIM28_colocalisation.R"


rule Fig4A_Ovarian_TFs_overlap_TRIM28_genome_bg:
    input:
        ATAC_control = f"{ATAC_TABLES}/peaks/consensus_per_condition/WT_consensus_peaks.txt",
        TRIM28       = config["TRIM28_peaks"],
        genome       = GENOME,
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/combinaison_TF_overlaps_genome_bg.csv",
        pdf   = f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap_genome_bg.pdf",
        png   = f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap_genome_bg.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4A_TF_TRIM28_colocalisation_genome_bg.R"


rule Fig4A_Ovarian_TFs_overlap_TRIM28_pie_chart:
    input:
        ATAC_control = f"{ATAC_TABLES}/peaks/consensus_per_condition/WT_consensus_peaks.txt",
        TRIM28       = config["TRIM28_peaks"],
        genome       = GENOME,
        **TF_INPUTS
    output:
        overlaps = directory(f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/pie_chart"),
        pdf      = f"{OUTPUT_PDF}/Fig4A_TF_TRIM28_overlap_pie_chart.pdf",
        png      = f"{OUTPUT_PNG}/Fig4A_TF_TRIM28_overlap_pie_chart.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4A_TF_TRIM28_colocalisation_pie_chart.R"


rule Fig4A_Ovarian_TFs_overlap_heatmap:
    input:
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/TF_overlaps/TF_overlap_regions_heatmap.csv",
        pdf   = f"{OUTPUT_PDF}/Fig4A_TF_colocalisation_heatmap.pdf",
        png   = f"{OUTPUT_PNG}/Fig4A_TF_colocalisation_heatmap.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4A_TF_colocalisation_heatmap.R"


rule Fig4A_Ovarian_TFs_overlap_upset:
    input:
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/TF_TRIM28_overlaps/TF_overlap_TRIM28.csv",
        pdf   = f"{OUTPUT_PDF}/Fig4A_TF_colocalisation_upset.pdf",
        png   = f"{OUTPUT_PNG}/Fig4A_TF_colocalisation_upset.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4A_TF_colocalisation_upset.R"


# =============================================================================
# Figure 4B - TF colocalisation with TRIM28
# =============================================================================

rule Fig4B_TF_colocalisation_TRIM28:
    input:
        TRIM28 = config["TRIM28_peaks"],
        **TF_INPUTS
    output:
        table = f"{OUTPUT_TABLES}/TF_overlaps/TF_overlap_regions_upset.csv",
        pdf   = f"{OUTPUT_PDF}/Fig4B_TF_colocalisation_TRIM28.pdf",
        png   = f"{OUTPUT_PNG}/Fig4B_TF_colocalisation_TRIM28.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/Result_4/Fig4B_TF_colocalisation_TRIM28.R"


# =============================================================================
# Figure 4C - ATAC, TRIM28 and TF genomic track examples
# =============================================================================

rule Fig4C_plot_genomic_tracks_example:
    input:
        ATAC_ctrl_bigwig = f"{ATAC_MERGED_BIGWIG_FOLDER}/WT.bw",
        ATAC_cKO_bigwig  = f"{ATAC_MERGED_BIGWIG_FOLDER}/KO.bw",
        ATAC_peaks       = f"{ATAC_TABLES}/ATAC_all_OCRs.bed",
        ATAC_clustering  = f"{ATAC_DAR_TABLES}/ATAC_DAR_heatmap_clusters.csv",
        TRIM28_peaks     = config["TRIM28_peaks"],
        genome           = GENOME,
        **TF_INPUTS
    params:
        region = lambda wildcards: GENOMIC_TRACK_REGIONS[wildcards.region_id]
    output:
        pdf = f"{OUTPUT_PDF}/ATAC_TRIM28_TF_tracks/{{region_id}}.pdf",
        png = f"{OUTPUT_PNG}/ATAC_TRIM28_TF_tracks/{{region_id}}.png"
    threads: 12
    resources:
        mem_mb = 64000
    script:
        "../scripts/Result_4/Fig4C_TF_TRIM28_track_examples.R"