#!/bin/bash


# Usage: sh mapping_summary.sh samplesheet.csv


# # Check if argument given
# if [[ "$#" -ne 3 ]] ; then
#     echo -e 'Please provide the samplesheet and the two output files.\nUsage: sh mapping_summary.sh samplesheet.csv'
#     exit 1
# fi

# Get sample name and final result file name
samples=(`awk -F "\"*,\"*" 'NR>1{print $1"_R"$2}' $1`)
# echo ${samples[*]}

results=$2

mapping=$3
trim=`echo $mapping/01_prealign/trimgalore`
mapping_res=`echo $mapping/02_alignment/bowtie2/target/markdup`

res_count_name=$4
res_pct_name=$5

# Column names
echo -e "Sample,Total,Trimmed,Mapped,Multimapped,Duplicated,Peaks" > $res_count_name
echo -e "Sample,% Trimmed, % Mapped*,% Multimapped**,% Duplicated**" > $res_pct_name

# For each sample
for sample in "${samples[@]}"; do
	# Get the total number of read paires in fastq files
	tot_reads=`grep "sequences processed in total" $trim/$sample\_2*report.txt | awk '{print $1}'`

	# Trimmed reads
	trimmed_out=`grep "Number of sequence pairs" $trim/$sample\_2*report.txt | awk '{print $19}'`
	trimmed=$(expr $tot_reads - $trimmed_out)
	pct_trimmed=`awk -v a=$trimmed -v b=$tot_reads -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Mapped reads
	unmapped=`grep $sample $results/multiqc_bowtie2.txt | awk '{print $4}'`
	mapped=$(expr $trimmed - $unmapped)
	pct_mapped=`awk -v a=$mapped -v b=$trimmed -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Duplicated reads
	dup=`grep "reads duplicated" $mapping_res/$sample.stats | awk '{print $4}'`
	pct_dup=`awk -v a=$dup -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Multimapped reads
	mult=`grep $sample $results/mqc_bowtie2_pe_plot_1.txt | awk '{print $3}' | sed 's/\.0//g'`
	pct_mult=`awk -v a=$mult -v b=$mapped -v OFMT="%.2f" 'BEGIN { print a*100/b }'`

	# Peaks
	peaks=`grep $sample $results/multiqc_primary_peakcounts_plot.txt | awk '{print $2}' | sed 's/\.0//g'`

	counts=`echo -e "$sample,$tot_reads,$trimmed,$mapped,$mult,$dup,$peaks" | sed 's/ //g'`
	percentages=`echo -e "$sample,$pct_trimmed,$pct_mapped,$pct_mult,$pct_dup"`

	echo $counts >> $res_count_name
	echo $percentages >> $res_pct_name

done
