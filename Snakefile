'''
Author: Isabelle Stévant
Affiliation: CNRS UMR 9002
Date: 26/08/2026
Licence: MIT

Snakemake rules to run the full pipeline.
The parameters of the analysis are defined in the analysis_parameters.yaml configuration file.

Pipeline created to analyse the ATAC-seq data from control and TRIM28 cKO adult mouse ovaries.
'''


# Import the different pipeline modules
include: "workflow/rules/0_prepare_CUTandRUN_data.smk"
include: "workflow/rules/0_prepare_ATAC_data.smk"
include: "workflow/rules/0_prepare_ChIP_SUMO_data.smk"
include: "workflow/rules/1_H3K9me3_landscape_Ctrl_ovaries.smk"
include: "workflow/rules/2_H3K9me3_Ctrl_vs_cKO_ovaries.smk"
include: "workflow/rules/3_ATAC_Ctrl_vs_cKO_ovaries.smk"
include: "workflow/rules/4_TRIM28_ovTF_colocalisation.smk"
include: "workflow/rules/5_TRIM28_TF_destabilisation.smk"


###########################################
#                                         #
#                  Rules                  #
#                                         #
###########################################

annotation_list = [
    "results/data/gencode.vM25.annotation.gtf.gz",
    "results/data/mm10_enhancers_Encode_Fantom5.bed",
    "results/data/mm10_repeats.bed"
]

# Run the whole pipeline
rule all:
    input:
        (
            annotation_list 
            + rule_CUTandRUN_input_list
            + rule_ATAC_input_list
            + rule_ChIP_SUMO_input_list
            + rule_Result1_input_list
            + rule_Result2_input_list
            + rule_Result3_input_list
            + rule_Result4_input_list
            + rule_Result5_input_list
        )

# Install the necessary R packages using Renv
rule install_packages:
    script:
        "renv/restore.R"

# Get mm10 annotation from Gencode
rule Get_genome:
    output:
        genome = "results/data/gencode.vM25.annotation.gtf.gz"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/Get_genome_annotation.R"

# Get enhancer coordinates from Encode and Fantom5 for mm10
rule Get_enhancer_catalogue:
    output:
        enhancers = "results/data/mm10_enhancers_Encode_Fantom5.bed"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/Get_enhancer_catalogue.R"

# Get repeatmasker for mm10
rule Get_repeatmasker:
    output:
        repeats = "results/data/mm10_repeats.bed"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "workflow/scripts/Get_repeatmasker.R"

