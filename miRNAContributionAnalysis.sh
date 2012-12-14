#!/bin/bash

Usage() {
	echo "Usage:"
	echo "ContributionAnalysis.sh <baseRef> <usableSequences> <scriptsPath>"
	echo "baseRef: basename of the analysis file."
	echo "usableSequences: List of valid RNA sequence (after trimming)."
	echo "scriptsPath: path containing the miRNA-Tools git clone."
	echo "			 (https://github.com/CharlesJB/miRNA-Tools)"
}

ValidateFile() {
	fileToCheck=$1
	if [ ! -e "$fileToCheck" ]
	then
		echo "Invalid file: $fileToCheck"
		exit
	fi
}

id=$1
usableSequences=$2
scriptsPath=$3

# 1. Calculate total number of sequences
echo "	Caculating total number of usable sequences."
numberOfSequences=$(grep '>' $usableSequences | awk '{ SUM += $4 } END { print SUM }')
echo "		$numberOfSequences sequences."

# 2. Calulate relative contribution based on total number of sequences
for file in $(ls | grep mature | grep fa)
do
	echo "	Calculating relative contribution for $file."
#	uniqueFastaOutput=$baseRef"_blast_ID_unique.fa"
	relativeOutput="lengthDist_"$id"_relative.txt"
	cat $file | $scriptsPath/LengthDistribution.py 60 $numberOfSequences > $relativeOutput

	# 3. Plot the relative distribution
	echo "	Plotting the relative distribution."
	cat $usableSequences | $scriptsPath/LengthDistribution.py 60 > lengthDist.tmp
	Rscript $scriptsPath/PlotRelativeLengthDistribution.R $relativeOutput lengthDist.tmp $baseRef | sed '/null device/d' | sed -e '/          1/d'
done

rm *.tmp
