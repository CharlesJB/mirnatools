#!/bin/bash

id=$1
path=$2

cut -f2- "mature_"$id"_looseSeqID.txt" | sed 's/\t/\n/g' | sort | uniq > $id"_looseSeqID.tmp"
cut -f2- "mature_"$id"_noMatchGaps_ID.txt" | sed 's/\t/\n/g' | sort | uniq > $id"_noMatchGaps_ID.tmp"
cut -f2- "mature_"$id"_noMatchLength_ID.txt" | sed 's/\t/\n/g' | sort | uniq > $id"_noMatchLength_ID.tmp"
cut -f2- "mature_"$id"_noMatchMismatches_ID.txt" | sed 's/\t/\n/g' | sort | uniq > $id"_noMatchMismatches_ID.tmp"
cut -f2- "mature_"$id"_perfectSeqID.txt" | sed 's/\t/\n/g' | sort | uniq > $id"_perfectSeqID.tmp"

for file in $(ls | grep 'ID\.tmp')
do
	name=${file%_*}
	name=${nameEnd%.*}
	name="mature_"$id$name".fa"
	$path/getFasta.py $file $id"_UsableSeq.fasta" > $name
	lengthDist="lengthDist"${name%.*}".txt"
	cat $name | $path/LengthDistribution.py 60 > $lengthDist
	Rscript $path/PlotLengthDistribution.R $lengthDist ${file%.*}
done
