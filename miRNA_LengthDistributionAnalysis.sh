#!/bin/bash

id=$1
path=$2

for file in $(ls | grep $id | grep '_ID' | grep txt) 
do
	tmpfile=${file%.*}".tmp"
	cut -f2- $file | sed 's/\t/\n/g' | sort | uniq > $tmpfile
	name=${file%_*}".fa"
	$path/getFasta.py $tmpfile $id"_UsableSeq.fasta" > $name
	lengthDist="lengthDist_"${name%.*}".txt"
	cat $name | $path/LengthDistribution.py 60 > $lengthDist
	Rscript $path/PlotLengthDistribution.R $lengthDist ${file%.*} | sed '/null device/d' | sed -e '/          1/d'
done

rm -f *.tmp
