#!/bin/bash

bigwig_folder=$1
merge_folder=$2
chromSize=$3
output_file=$4

mkdir $merge_folder

function prefix
{
    for file in $bigwig_folder/*.bw
    do
        echo `basename ${file%_*.bw}`
    done | sort -u
}

prefix | while read prefix
do
	rm $output_file
	wiggletools median $bigwig_folder/${prefix}_*.bw  > $merge_folder/${prefix}.wig && \
	wigToBigWig -clip $merge_folder/${prefix}.wig $chromSize $merge_folder/${prefix}.bw && \
	rm $merge_folder/${prefix}.wig && \
	echo "OK" > $output_file
done
