#!/bin/bash

id=$1

for file in $(ls | grep $id | grep '_ID' | grep txt) 
do
	tmpfile=${file%.*}".tmp"
	cut -f2- $file | sed 's/\t/\n/g' | sort | uniq > $tmpfile
	name=${file%_*}".fa"
	$MIRNA_TOOLS_PATH/getFasta.py $tmpfile $id"_UsableSeq.fasta" > $name
	lengthDist="lengthDist_"${name%.*}".txt"
	cat $name | $MIRNA_TOOLS_PATH/LengthDistribution.py 60 > $lengthDist
	Rscript $MIRNA_TOOLS_PATH/Rscripts/PlotLengthDistribution.R $lengthDist ${file%.*} | sed '/null device/d' | sed -e '/          1/d'
done

rm -f *.tmp
