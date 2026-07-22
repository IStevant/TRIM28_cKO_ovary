#!/bin/bash

bigwig_folder=$1
AB=$2
merge_folder=$3
chromSize=$4
output_file=$5


function prefix
{
	for file in $bigwig_folder/*${AB}_R*.bw
	do
		echo `basename ${file%_R*.bw}`
	done | sort -u
}

prefix | while read prefix
do
	echo ${prefix}
	rm $output_file
	wiggletools median $bigwig_folder/${prefix}_R*.bw  > $merge_folder/${prefix}.wig && \
	wigToBigWig -clip $merge_folder/${prefix}.wig $chromSize $merge_folder/${prefix}.bw && \
	rm $merge_folder/${prefix}.wig && \
	echo "OK" > $output_file
done
