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

PROCESSED_DATA = "results/CR/processed_data"
OUTPUT_TABLES  = "results/CR/tables"
OUTPUT_PNG     = "results/CR/graphs/PNG"
OUTPUT_PDF     = "results/CR/graphs/PDF"
TMP            = "results/CR/tmp"

AB             = config["CR_nf-core_AB"]
CONDITION      = config["CR_nf-core_cond"]
MAPPING_FOLDER = config["CR_nf-core_path"]
SAMPLESHEET    = config["CR_nf-core_samplesheet"]
CR_GROUPS      = expand("{condition}_{AB}", condition=CONDITION, AB=AB)


MULTIQC        = f"{MAPPING_FOLDER}/{config['CR_nf-core_multiQC']}"
RAW_PEAKS      = f"{MAPPING_FOLDER}/{config['CR_nf-core_rawPeaks']}"
BAM            = f"{MAPPING_FOLDER}/{config['CR_nf-core_bam']}"
BIGWIG         = f"{MAPPING_FOLDER}/{config['CR_nf-core_bigwig']}"

CR_NORM_BIGWIG_FOLDER = f"{PROCESSED_DATA}/CR_bigwig"
MERGED_BIGWIG_FOLDER  = f"{CR_NORM_BIGWIG_FOLDER}/merged"

PEAK_TYPES = ["raw_peaks", "domains"] if config["merge_peaks"] == "yes" else ["raw_peaks"]
TFBS_BACKGROUND = ["genome", "conditions"]


# =============================================================================
# Main output list
# =============================================================================

rule_CUTandRUN_input_list = [
    f"{OUTPUT_PDF}/CR_QC_global.pdf",
    expand(f"{OUTPUT_PDF}/CR_QC_detailed_{{group}}.pdf", group=CR_GROUPS),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_raw_counts.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_norm_counts.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_samplesheet.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_size_factors.csv", AB=AB),
    expand(f"{OUTPUT_TABLES}/matrices/{{AB}}_all_regions.bed", AB=AB),
    expand(f"{MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt", AB=AB)
]

# =============================================================================
# Peak and domain processing
# =============================================================================

rule CR_Merge_adjacent_peaks:
    params:
        raw_peaks = RAW_PEAKS,
        distance  = config["CR_merge_peaks"]
    output:
        output_folder = directory(f"{OUTPUT_TABLES}/peaks/domains"),
        peak_counts   = f"{OUTPUT_TABLES}/peaks/domains/domain_count.csv"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/CR/CR_merge_adjacent_peaks.R"


rule CR_Get_consensus_regions_per_condition:
    input:
        raw_peaks = RAW_PEAKS,
        domains   = rules.CR_Merge_adjacent_peaks.output.output_folder,
        genome    = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        condition   = "{condition}_{AB}",
        promoter    = config["CR_promoter_distance"],
        distance    = config["CR_merge_peaks"],
        minimum_rep = config["CR_min_peak_rep"],
        tmp         = TMP
    output:
        consensus_peak_table   = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_peaks.txt",
        consensus_domain_table = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_domains.txt",
        consensus_peak_bed     = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_peaks.bed",
        consensus_domain_bed   = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_domains.bed"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/CR/CR_get_consensus_per_condition.R"


rule CR_Get_peaks_domains_count_per_condition:
    input:
        peaks   = expand(f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_peaks.bed", AB=AB, condition=CONDITION),
        domains = expand(f"{OUTPUT_TABLES}/peaks/consensus_per_condition/{{condition}}_{{AB}}_consensus_domains.bed", AB=AB, condition=CONDITION)
    output:
        peak_count   = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/peak_count.csv",
        domain_count = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/domain_count.csv"
    shell:
        """
        wc -l {input.peaks} |
            awk 'BEGIN{{OFS=","}} $NF!="total"{{f=$NF; sub(".*/","",f); print f,$1}}' \
            > {output.peak_count}

        wc -l {input.domains} |
            awk 'BEGIN{{OFS=","}} $NF!="total"{{f=$NF; sub(".*/","",f); print f,$1}}' \
            > {output.domain_count}
        """


rule CR_Get_consensus_regions_per_AB:
    input:
        raw_peaks = RAW_PEAKS,
        domains   = rules.CR_Merge_adjacent_peaks.output.output_folder,
        genome    = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        AB          = "{AB}",
        promoter    = config["CR_promoter_distance"],
        distance    = config["CR_merge_peaks"],
        minimum_rep = config["CR_min_peak_rep"],
        tmp         = TMP
    output:
        consensus_peak_table   = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_peaks.txt",
        consensus_domain_table = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_domains.txt",
        consensus_peak_bed     = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_peaks.bed",
        consensus_domain_bed   = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_domains.bed",
        gtf_files              = f"{TMP}/{{AB}}_consensus_peaks.gtf"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/CR/CR_get_consensus_per_AB.R"


rule CR_Get_peaks_domains_count_per_AB:
    input:
        peaks   = expand(f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_peaks.bed", AB=AB),
        domains = expand(f"{OUTPUT_TABLES}/peaks/consensus_per_AB/{{AB}}_consensus_domains.bed", AB=AB)
    output:
        peak_count   = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/peak_count.csv",
        domain_count = f"{OUTPUT_TABLES}/peaks/consensus_per_AB/domain_count.csv"
    shell:
        """
        wc -l {input.peaks} |
            awk 'BEGIN{{OFS=","}} $NF!="total"{{f=$NF; sub(".*/","",f); print f,$1}}' \
            > {output.peak_count}

        wc -l {input.domains} |
            awk 'BEGIN{{OFS=","}} $NF!="total"{{f=$NF; sub(".*/","",f); print f,$1}}' \
            > {output.domain_count}
        """


# =============================================================================
# Mapping and quality-control reports
# =============================================================================

rule CR_Mapping_summary:
    input:
        samplesheet    = SAMPLESHEET,
        mapping_folder = MAPPING_FOLDER,
        multiQC        = MULTIQC
    output:
        counts      = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv",
        percentages = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_percent.csv"
    threads: 1
    resources:
        mem_mb = 12000
    shell:
        """
        sh workflow/scripts/CR/CR_get_mapping_summary.sh \
            {input.samplesheet} \
            {input.multiQC} \
            {input.mapping_folder} \
            {output.counts} \
            {output.percentages}
        """


rule CR_Global_mapping_report:
    input:
        summary_counts = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv",
        merged_peaks   = f"{OUTPUT_TABLES}/peaks/domains/domain_count.csv",
        peak_count     = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/peak_count.csv",
        domain_count   = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/domain_count.csv"
    output:
        pdf_report = f"{OUTPUT_PDF}/CR_QC_global.pdf",
        png_report = f"{OUTPUT_PNG}/CR_QC_global.png"
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/CR/CR_QC_raw_data_global.R"


rule CR_Detailed_mapping_report:
    input:
        summary_counts  = f"{OUTPUT_TABLES}/mapping_summary/mapping_summary_count.csv",
        raw_peaks       = RAW_PEAKS,
        consensus_peaks = f"{OUTPUT_TABLES}/peaks/consensus_per_condition/peak_count.csv",
        repeatMasker    = "results/data/mm10_repeats.bed",
        genome          = "results/data/gencode.vM25.annotation.gtf.gz"
    params:
        mapping_folder = MAPPING_FOLDER,
        ctrl_AB        = config["CR_nf-core_AB_Ctrl"],
        min_peak_rep   = config["CR_min_peak_rep"],
        promoter       = config["CR_promoter_distance"]
    output:
        pdf = expand(f"{OUTPUT_PDF}/CR_QC_detailed_{{group}}.pdf", group=CR_GROUPS),
        png = expand(f"{OUTPUT_PNG}/CR_QC_detailed_{{group}}.png", group=CR_GROUPS)
    threads: 1
    resources:
        mem_mb = 12000
    script:
        "../scripts/CR/CR_QC_raw_data_detailed.R"


# =============================================================================
# Quantification matrices
# =============================================================================

rule CR_Get_read_coverage:
    input:
        consensus_gtf = f"{TMP}/{{AB}}_consensus_peaks.gtf"
    params:
        AB        = "{AB}",
        bam_files = f"{BAM}/*{{AB}}*.bam"
    output:
        counts = f"{OUTPUT_TABLES}/featureCount/{{AB}}_featureCounts.csv"
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        featureCounts \
            -p \
            --countReadPairs \
            -a {input.consensus_gtf} \
            -F SAF \
            -o {output.counts} \
            {params.bam_files}
        """


rule CR_Filter_matrices:
    input:
        samplesheet = SAMPLESHEET,
        counts      = f"{OUTPUT_TABLES}/featureCount/{{AB}}_featureCounts.csv"
    params:
        minReads = config["CR_minReads"]
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
        "../scripts/CR/CR_prepare_matrices.R"


# =============================================================================
# BigWig normalisation and merging
# =============================================================================

rule CR_Normalize_bigwig:
    input:
        size_factors = f"{OUTPUT_TABLES}/matrices/{{AB}}_size_factors.csv"
    params:
        AB                = "{AB}",
        bigwig_folder     = BIGWIG,
        new_bigwig_folder = CR_NORM_BIGWIG_FOLDER
    output:
        output_file = f"{CR_NORM_BIGWIG_FOLDER}/{{AB}}_size_factors.csv"
    threads: 2
    resources:
        mem_mb = 24000
    script:
        "../scripts/CR/CR_norm_bigwig.R"


rule CR_Merge_norm_bigwig:
    input:
        output_file = f"{CR_NORM_BIGWIG_FOLDER}/{{AB}}_size_factors.csv"
    params:
        AB            = "{AB}",
        bigwig_folder = CR_NORM_BIGWIG_FOLDER,
        merged_folder = MERGED_BIGWIG_FOLDER,
        chom_size     = config["chom_size"]
    output:
        output_file = f"{MERGED_BIGWIG_FOLDER}/{{AB}}_merge_bw.txt"
    threads: 1
    resources:
        mem_mb = 16000
    shell:
        """
        sh workflow/scripts/CR/CR_merge_bigwig.sh \
            {params.bigwig_folder} \
            {params.AB} \
            {params.merged_folder} \
            {params.chom_size} \
            {output.output_file}
        """