#!/bin/bash

fastq_file=$1
adaptor_sequence=$2
scripts_path=$3

trim_output=$(basename $fastq_file)
trim_output="RawSequenceFiltering/"${trim_output%.*}".trimmed.fastq"
qual_output=${trim_output%.*}".qual.fastq"
copy_output=${trim_output%.*}".copy.fastq"

mkdir -p RawSequenceFiltering

cat $fastq_file | $scripts_path/AdaptorFinder2.py $adaptor_sequence 4 3 > $trim_output
cat $trim_output | $scripts_path/AdaptorQualFilter.py 9 > $qual_output
$scripts_path/CopyKeeper.py $qual_output 4 fasta > $copy_output

count=$(grep '>' $copy_output | wc -l)
echo "$count unique sequences!"
