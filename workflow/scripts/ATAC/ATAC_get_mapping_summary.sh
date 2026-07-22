#!/bin/bash


# Usage: sh mapping_summary.sh samplesheet.csv


# # Check if argument given
# if [[ "$#" -ne 3 ]] ; then
#     echo -e 'Please provide the samplesheet and the two output files.\nUsage: sh mapping_summary.sh samplesheet.csv'
#     exit 1
# fi

# Get sample name and final result file name
samples=(`awk -F "\"*,\"*" 'NR>1{print $1}' $1 | sed 's/_T1//g' | uniq`)
# echo ${samples[*]}

results=$2
# echo $results
res_count_name=$3
res_pct_name=$4

# Column names
echo -e "Sample,Total,Kept after trimming,Mapped,Duplicated,Mitochondrial,Kept after filtering,Peaks,FRiP" > $res_count_name
echo -e "Sample,% kept after trimming*,% mapped**,% duplicated***,% mitochondrial***,% kept after filtering***" > $res_pct_name

# For each sample
for sample in "${samples[@]}"; do
	# Get the total number of reads in fastq files
	tot_reads=`grep $sample $results/multiqc_fastqc.txt | awk -F'\t' '{print $5}' | sed 's/\.0//g' | jq -s 'add'`

	# Reads after trimming
	trimmed=`grep $sample $results/multiqc_samtools_stats_1.txt | awk '{print $2}' | sed 's/\.0//g'`
	pct_trimmed=`awk -v a=$trimmed -v b=$tot_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Mapped reads
	mapped=`grep $sample $results/multiqc_samtools_stats_1.txt | awk '{print $8}' | sed 's/\.0//g'`
	pct_mapped=`awk -v a=$mapped -v b=$trimmed -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Duplicated reads
	dup=`grep $sample $results/multiqc_samtools_stats_1.txt | awk '{print $13}' | sed 's/\.0//g'`
	pct_dup=`awk -v a=$dup -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Mitochondtial reads
	mt=`grep $sample $results/mqc_samtools-idxstats-mapped-reads-plot-2_Raw_Counts.txt | awk '{print $21}' | sed 's/\.0//g'`
	pct_mt=`awk -v a=$mt -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Filtered and mapped reads
	filtered=`grep $sample $results/multiqc_samtools_stats_2.txt | awk '{print $2}' | sed 's/\.0//g'`
	pct_filtered=`awk -v a=$filtered -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Peaks
	peaks=`awk 'NR>1{print $0}' $results/multiqc_mlib_peak_count-plot.txt | grep $sample | awk '{$1=""; print $0}' | sed 's/\t//g' | sed 's/\.0//g'`
	FRiP=`awk 'NR>1{print $0}' $results/multiqc_mlib_frip_score-plot.txt | grep $sample  | awk '{$1=""; print $0}' | sed 's/\t//g'`
	

	counts=`echo -e "$sample,$tot_reads,$trimmed,$mapped,$dup,$mt,$filtered,$peaks,$FRiP" | sed 's/ //g'`
	percentages=`echo -e "$sample,$pct_trimmed,$pct_mapped,$pct_dup,$pct_mt,$pct_filtered"`

	echo $counts >> $res_count_name
	echo $percentages >> $res_pct_name

done
