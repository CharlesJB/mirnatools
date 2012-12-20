#!/bin/bash

reference=$1
usableSequences=$2
evalue=$3

if [ "$evalue" == "" ]
then
	evalue="0.01"
fi

if [ -e "$reference" ]
then
	# 1. Create Reference directory
	echo "	Creating reference directory."
	if [ -d "Reference" ]
	then
		rm -rf Reference
	fi
	mkdir Reference

	# 2. Link reference
	echo "	Linking reference."
	ln -s ../$reference Reference/

	# 3. Format database
	echo "	Formating database."
	input="Reference/$(basename $reference)"
	out=${input%.*}
	makeblastdb -in $input -dbtype 'nucl' -out $out

	# 4. Do the actual blast
	echo ""
	echo "	Doing the actual blast."
	out=$(basename ${reference%.*}).out
	blastn -task blastn-short -db ${input%.*} -query $usableSequences -out $out -word_size 4 -evalue $evalue -num_threads 4
else
	echo "Cannot find reference data: $reference"
	echo "Usage:"
	echo "BlastAnalysis.sh <referenceData> <usableSequences> <evalue>"
	echo "referenceData: Complete path to reference data for the blast analysis. Must be in fasta format."
	echo "usableSequences: List of valid RNA sequence (after trimming)."
	echo "evalue: threshold score for blast analysis (default = 0.01)."
fi
