#!/bin/bash

# Parse Mature
count=$(grep '>' $(find GenomeMapping/mature_hsapiens/* | grep '\.fasta') | awk '{sum += $4} END {print sum}')
echo -e "miRNA\t$count"

# Parse Others
for file in $(find GenomeMapping/bam/*)
do 
	count=$((count+1))
	name=$(basename ${file%.*})
	count=$(samtools view -F4 $file | grep 'NM:i:0' | cut -f1 | sort | uniq | ~/git-clones/miRNA-Tools/Utility/getFasta.py RawSequenceFiltering/*.fasta match | grep '>' | awk '{sum += $4} END { print sum }')
	if [ "$count" == "" ]
	then
		count=0
	fi
	echo -e "$name\t$count"
done
