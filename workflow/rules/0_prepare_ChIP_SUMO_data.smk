"""
Author: Isabelle Stévant
Affiliation: CNRS UMR 9002
Date: 14/05/2026
Licence: MIT
"""

configfile: "analysis_parameters.yaml"


# =============================================================================
# Configuration
# =============================================================================

PROCESSED_DATA = "results/ChIP/processed_data"
OUTPUT_TABLES  = "results/ChIP/tables"
OUTPUT_PNG     = "results/ChIP/graphs/PNG"
OUTPUT_PDF     = "results/ChIP/graphs/PDF"
TMP            = "results/ChIP/tmp"

AB             = config["ChIP_nf-core_AB"]
CONDITION      = config["ChIP_nf-core_cond"]
MAPPING_FOLDER = config["ChIP_nf-core_path"]
SAMPLESHEET    = config["ChIP_nf-core_samplesheet"]

RAW_PEAKS = f"{MAPPING_FOLDER}/{config['ChIP_nf-core_rawPeaks']}"
BIGWIG    = f"{MAPPING_FOLDER}/{config['ChIP_nf-core_bigwig']}"

PEAK_TYPE = "raw_peaks"

CHIP_NORM_BIGWIG_FOLDER = f"{PROCESSED_DATA}/ChIP_bigwig"
MERGED_BIGWIG_FOLDER    = f"{CHIP_NORM_BIGWIG_FOLDER}/merged"


wildcard_constraints:
    condition = "|".join(CONDITION),
    AB = "|".join(AB)


# =============================================================================
# Main output list
# =============================================================================

rule_ChIP_input_list = [
    expand(f"{OUTPUT_TABLES}/peaks/{{condition}}_{{AB}}_consensus_peaks.txt", AB=AB, condition=CONDITION),
    expand(f"{OUTPUT_TABLES}/peaks/{{condition}}_{{AB}}_consensus_peaks.bed", AB=AB, condition=CONDITION),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_raw_counts.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_norm_counts.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_samplesheet.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_size_factors.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_all_regions.bed", AB=AB),
    expand(f"{MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt", AB=AB),
]


# =============================================================================
# Peak processing
# =============================================================================

rule ChIP_Get_consensus_regions_per_condition:
    input:
        raw_peaks = RAW_PEAKS,
        genome    = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        condition   = "{condition}_{AB}",
        promoter    = config["CR_promoter_distance"],
        minimum_rep = config["ChIP_min_peak_rep"],
        tmp         = TMP
    output:
        consensus_peak_table = f"{OUTPUT_TABLES}/peaks/{{condition}}_{{AB}}_consensus_peaks.txt",
        consensus_peak_bed   = f"{OUTPUT_TABLES}/peaks/{{condition}}_{{AB}}_consensus_peaks.bed"
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/ChIP/ChIP_get_consensus_per_condition.R"
# =============================================================================
# Quantification matrices
# =============================================================================

rule ChIP_Filter_matrices:
    input:
        samplesheet = SAMPLESHEET,
        counts      = f"{RAW_PEAKS}/{{AB}}/{{AB}}.consensus_peaks.featureCounts.txt"
    params:
        minReads = config["ChIP_minReads"]
    output:
        counts       = f"{OUTPUT_TABLES}/matrices/{{AB}}_raw_counts.csv",
        norm_counts  = f"{OUTPUT_TABLES}/matrices/{{AB}}_norm_counts.csv",
        samplesheet  = f"{OUTPUT_TABLES}/matrices/{{AB}}_samplesheet.csv",
        size_factors = f"{OUTPUT_TABLES}/matrices/{{AB}}_size_factors.csv",
        bed          = f"{OUTPUT_TABLES}/matrices/{{AB}}_all_regions.bed"
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/ChIP/ChIP_prepare_matrices.R"

# =============================================================================
# BigWig normalisation and merging
# =============================================================================

rule ChIP_Normalize_bigwig:
    input:
        size_factors = f"{OUTPUT_TABLES}/matrices/{{AB}}_size_factors.csv"
    params:
        AB                = "{AB}",
        bigwig_folder     = BIGWIG,
        new_bigwig_folder = CHIP_NORM_BIGWIG_FOLDER
    output:
        output_file = f"{CHIP_NORM_BIGWIG_FOLDER}/{{AB}}_size_factors.csv"
    threads: 2
    resources:
        mem_mb = 24000
    script:
        "../scripts/ChIP/ChIP_norm_bigwig.R"


rule ChIP_Merge_norm_bigwig:
    input:
        output_file = f"{CHIP_NORM_BIGWIG_FOLDER}/{{AB}}_size_factors.csv"
    params:
        AB            = "{AB}",
        bigwig_folder = CHIP_NORM_BIGWIG_FOLDER,
        merged_folder = MERGED_BIGWIG_FOLDER,
        chom_size     = config["chom_size"]
    output:
        output_file = f"{MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt"
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        sh workflow/scripts/ChIP/ChIP_merge_bigwig.sh \
            {params.bigwig_folder} \
            {params.AB} \
            {params.merged_folder} \
            {params.chom_size} \
            {output.output_file}
        """
