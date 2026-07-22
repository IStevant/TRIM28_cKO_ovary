"""
Author: Isabelle Stévant
Affiliation: CNRS UMR 9002
Date: 13/05/2026
Licence: MIT
"""

configfile: "analysis_parameters.yaml"


# =============================================================================
# Configuration
# =============================================================================

PROCESSED_DATA = "results/ATAC/processed_data"
OUTPUT_TABLES  = "results/ATAC/tables"
OUTPUT_PNG     = "results/ATAC/graphs/PNG"
OUTPUT_PDF     = "results/ATAC/graphs/PDF"
TMP            = "results/ATAC/tmp"

FROM_NFCORE = config["ATAC_from_nf-core"] == "yes"

ATAC_NORM_BIGWIG_FOLDER = f"{PROCESSED_DATA}/ATAC_bigwig"
MERGED_BIGWIG_FOLDER    = f"{ATAC_NORM_BIGWIG_FOLDER}/merged"

GENOME = config.get("genome_version", "mm10")

TFBS_BACKGROUND = ["genome", "conditions"]

ATAC_GROUPS = config["ATAC_nf-core_cond"]

# =============================================================================
# Input files
# =============================================================================

if FROM_NFCORE:
    MAPPING_FOLDER = config["ATAC_nf-core_path"]
    SAMPLESHEET    = config["ATAC_nf-core_samplesheet"]
    PEAK_TYPE      = config["ATAC_nf-core_peak_type"]

    MULTIQC = f"{MAPPING_FOLDER}/{config['ATAC_nf-core_multiQC'].format(peak_type=PEAK_TYPE)}"
    FEATURECOUNTS = f"{MAPPING_FOLDER}/{config['ATAC_nf-core_featureCounts'].format(peak_type=PEAK_TYPE)}"
    RAW_PEAKS = f"{MAPPING_FOLDER}/{config['ATAC_nf-core_rawPeaks'].format(peak_type=PEAK_TYPE)}"
    BAM = f"{MAPPING_FOLDER}/{config['ATAC_nf-core_bam']}"
    BIGWIG = f"{MAPPING_FOLDER}/{config['ATAC_nf-core_bigwig']}"

else:
    FEATURECOUNTS = config["ATAC_featureCounts"]
    RAW_PEAKS = config["ATAC_rawPeaks"]
    BAM = config["ATAC_bam"]
    BIGWIG = config["ATAC_bigwig"]
    SAMPLESHEET = config["ATAC_samplesheet"]
    PEAK_TYPE = config.get("ATAC_peak_type", "narrowPeak")
    MAPPING_FOLDER = config.get("ATAC_mapping_folder", "")


# =============================================================================
# Main output list
# =============================================================================

rule_ATAC_input_list = [
    f"{OUTPUT_PNG}/ATAC_filtered_region_comparison.png",
    f"{OUTPUT_PDF}/ATAC_filtered_region_comparison.pdf",
    directory(f"{OUTPUT_TABLES}/peaks/consensus_per_condition"),
    f"{OUTPUT_TABLES}/ATAC_raw_counts.csv",
    f"{OUTPUT_TABLES}/ATAC_norm_counts.csv",
    f"{OUTPUT_TABLES}/ATAC_samplesheet.csv",
    f"{OUTPUT_TABLES}/ATAC_size_factors.csv",
    f"{OUTPUT_TABLES}/ATAC_all_OCRs.bed",
    f"{ATAC_NORM_BIGWIG_FOLDER}/ATAC_size_factors.csv",
    f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt",
]

if FROM_NFCORE:
    rule_ATAC_input_list.extend([
        f"{OUTPUT_PDF}/ATAC_raw_data_QC_global.pdf",
        f"{OUTPUT_PNG}/ATAC_raw_data_QC_global.png",
        expand(f"{OUTPUT_PDF}/ATAC_raw_data_detailed_QC_{{group}}.pdf", group=ATAC_GROUPS)
    ])


# =============================================================================
# Peak processing
# =============================================================================

rule ATAC_Get_consensus_regions:
    input:
        raw_peaks = RAW_PEAKS
    params:
        promoter    = config["ATAC_promoter_distance"],
        minimum_rep = config["ATAC_min_peak_rep"],
        tmp         = TMP
    output:
        anno_folder = directory(f"{OUTPUT_TABLES}/peaks/consensus_per_condition")
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/ATAC/ATAC_get_annotated_peaks.R"


# =============================================================================
# Mapping and quality-control reports
# =============================================================================

rule ATAC_Mapping_summary:
    input:
        samplesheet = SAMPLESHEET,
        multiQC     = MULTIQC
    output:
        counts      = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv",
        percentages = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_percent.csv"
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        """
        sh workflow/scripts/ATAC/ATAC_get_mapping_summary.sh \
            {input.samplesheet} \
            {input.multiQC} \
            {output.counts} \
            {output.percentages}
        """


rule ATAC_Global_mapping_report:
    input:
        summary_counts = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv"
    params:
        mapping_folder = MAPPING_FOLDER,
        peak_type      = PEAK_TYPE
    output:
        pdf_report = f"{OUTPUT_PDF}/ATAC_raw_data_QC_global.pdf",
        png_report = f"{OUTPUT_PNG}/ATAC_raw_data_QC_global.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/ATAC/ATAC_QC_raw_data_global.R"



rule ATAC_Detailed_mapping_report:
    input:
        summary_counts = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv"
    params:
        mapping_folder  = MAPPING_FOLDER,
        peak_type       = PEAK_TYPE,
        raw_peak_folder = RAW_PEAKS,
        min_peak_rep    = config["ATAC_min_peak_rep"]
    output:
        pdf = expand(f"{OUTPUT_PDF}/ATAC_raw_data_detailed_QC_{{group}}.pdf", group=ATAC_GROUPS),
        png = expand(f"{OUTPUT_PNG}/ATAC_raw_data_detailed_QC_{{group}}.png", group=ATAC_GROUPS)
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/ATAC/ATAC_QC_raw_data_detailed.R"


# =============================================================================
# Quantification matrices
# =============================================================================

rule ATAC_Filter_matrices:
    input:
        samplesheet = SAMPLESHEET,
        counts      = FEATURECOUNTS
    params:
        minReads = config["ATAC_minReads"]
    output:
        counts       = f"{OUTPUT_TABLES}/ATAC_raw_counts.csv",
        norm_counts  = f"{OUTPUT_TABLES}/ATAC_norm_counts.csv",
        samplesheet  = f"{OUTPUT_TABLES}/ATAC_samplesheet.csv",
        size_factors = f"{OUTPUT_TABLES}/ATAC_size_factors.csv",
        bed          = f"{OUTPUT_TABLES}/ATAC_all_OCRs.bed"
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/ATAC/ATAC_prepare_matrices.R"


# =============================================================================
# BigWig normalisation and merging
# =============================================================================

rule ATAC_Normalize_bigwig:
    input:
        size_factors = f"{OUTPUT_TABLES}/ATAC_size_factors.csv"
    params:
        bigwig_folder     = BIGWIG,
        new_bigwig_folder = ATAC_NORM_BIGWIG_FOLDER
    output:
        output_file = f"{ATAC_NORM_BIGWIG_FOLDER}/ATAC_size_factors.csv"
    threads: 12
    resources:
        mem_mb = 16000
    script:
        "../scripts/ATAC/ATAC_norm_bigwig.R"


rule ATAC_Merge_norm_bigwig:
    input:
        size_factors = f"{ATAC_NORM_BIGWIG_FOLDER}/ATAC_size_factors.csv"
    params:
        bigwig_folder = ATAC_NORM_BIGWIG_FOLDER,
        merged_folder = MERGED_BIGWIG_FOLDER,
        chom_size     = config["chom_size"]
    output:
        output_file = f"{ATAC_NORM_BIGWIG_FOLDER}/merge_bw.txt"
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        sh workflow/scripts/ATAC/ATAC_merge_bigwig.sh \
            {params.bigwig_folder} \
            {params.merged_folder} \
            {params.chom_size} \
            {output.output_file}
        """


# =============================================================================
# ATAC-seq QC and comparison plots
# =============================================================================

rule ATAC_Plot_OCR_comparison:
    input:
        counts    = f"{OUTPUT_TABLES}/ATAC_raw_counts.csv",
        norm_data = f"{OUTPUT_TABLES}/ATAC_norm_counts.csv"
    params:
        raw_peak_folder = RAW_PEAKS,
        peak_type       = PEAK_TYPE,
        min_peak_rep    = config["ATAC_min_peak_rep"],
        minReads        = config["ATAC_minReads"],
        corr_method     = config["ATAC_corr_met"],
        promoter        = config["ATAC_promoter_distance"]
    output:
        pdf = f"{OUTPUT_PDF}/ATAC_filtered_region_comparison.pdf",
        png = f"{OUTPUT_PNG}/ATAC_filtered_region_comparison.png"
    threads: 1
    resources:
        mem_mb = 64000
    script:
        "../scripts/ATAC/ATAC_filtered_region_comparison.R"