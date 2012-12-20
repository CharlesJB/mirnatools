#!/bin/bash

Usage() {
	echo ""
	echo "Usage:"
	echo "RawSequenceFiltering.sh <fastq_file> <adaptor_sequence>" 
	echo "    fastq_file: Raw data in fastq format."
	echo "    adaptor_sequence: Sequence of the adaptor to remove."
	echo ""
	echo "Note: You must set the MIRNA_TOOLS_PATH environment variable to miRNA-Tools folder"
	echo "      This can be done either manually or by using LoadModules, i.e.:"
	echo "      source ~/git-clones/miRNA-Tools LoadModules"
	echo ""
}

# Set variables
source $MIRNA_TOOLS_PATH/Utility/utility.sh
EXPECTED_ARGS=2

# Load parameters
fastq_file=$1
adaptor_sequence=$2

if [ $# -eq $EXPECTED_ARGS ]
then
	# Validate files
	ValidateFile $fastq_file
	ValidateFile $MIRNA_TOOLS_PATH/AdaptorFinder2.py
	ValidateFile $MIRNA_TOOLS_PATH/AdaptorQualFilter.py
	ValidateFile $MIRNA_TOOLS_PATH/CopyKeeper.py

	# Set output filenames
	trim_output=$(basename $fastq_file)
	trim_output="RawSequenceFiltering/"${trim_output%.*}".trimmed.fastq"
	qual_output=${trim_output%.*}".qual.fastq"
	copy_output=${trim_output%.*}".copy.fasta"

	# Do the actual analysis
	mkdir -p RawSequenceFiltering

	cat $fastq_file | AdaptorFinder2.py $adaptor_sequence 4 3 > $trim_output
	cat $trim_output | AdaptorQualFilter.py 9 > $qual_output
	CopyKeeper.py $qual_output 4 fasta > $copy_output

	# Show some metadata
	count=$(grep '>' $copy_output | wc -l)
	echo "$count unique sequences!"
else
	Usage
fi
