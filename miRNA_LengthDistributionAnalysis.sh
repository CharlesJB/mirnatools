#!/bin/bash

id=$1
path=$2

echo id:$id
for file in $(ls | grep $id | grep '_ID' | grep txt) 
do
	echo file:$file
	tmpfile=${file%.*}".tmp"
	cut -f2- $file | sed 's/\t/\n/g' | sort | uniq > $tmpfile
	name=${file%_*}".fa"
	$path/getFasta.py $file $id"_UsableSeq.fasta" > $name
	lengthDist="lengthDist_"${name%.*}".txt"
	cat $name | $path/LengthDistribution.py 60 > $lengthDist
	Rscript $path/PlotLengthDistribution.R $lengthDist ${file%.*}
done

rm *.tmp
