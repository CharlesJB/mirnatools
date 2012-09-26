#!/bin/bash

reference=$1
usableSequences=$2

if [ -e "$reference" ]
then
	# 1. Create Reference directory
    echo "    Creating reference directory."
	if [ -d "Reference" ]
	then
		rm -rf Reference
	fi
	mkdir Reference

	# 2. Link reference
    echo "    Linking reference."
	ln -s ../$reference Reference/

	# 3. Format database
    echo "    Formating database."
	input="Reference/$(basename $reference)"
	out=${input%.*}
	makeblastdb -in $input -dbtype 'nucl' -out $out

	# 4. Do the actual blast
    echo "    Doing the actual blast."
	out=$(basename ${reference%.*}).out
	blastn -task blastn-short -db ${input%.*} -query $usableSequences -out $out -word_size 4 -evalue 0.01 -num_threads 2
else
	echo "Cannot find reference data: $reference"
	echo "Usage:"
	echo "BlastAnalysis.sh <referenceData>"
	echo "referenceData: Complete path to reference data for the blast analysis. Must be in fasta format."
fi
