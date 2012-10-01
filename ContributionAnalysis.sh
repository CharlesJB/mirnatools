#!/bin/bash

Usage() {
	echo "Usage:"
	echo "ContributionAnalysis.sh <blastOutput> <idToRemove> <usableSequences> <scriptsPath>"
	echo "blastOutput: Relative path to blast blastIdOutput."
	echo "idToRemove: List of IDs that match known miRNA. We need to remove them from the results."
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

baseRef=$1
usableSequences=$2
scriptsPath=$3

# 1. Calculate total number of sequences
echo "	Caculating total number of usable sequences."
numberOfSequences=$(grep '>' $usableSequences | awk '{ SUM += $4 } END { print SUM }')
echo "		$numberOfSequences sequences."

# 2. Calulate relative contribution based on total number of sequences
echo "	Calculating relative contribution."
uniqueFastaOutput=$baseRef"_blast_ID_unique.fa"
relativeOutput="lengthDist_"$baseRef"_unique_relative.txt"
cat $uniqueFastaOutput | $scriptsPath/LengthDistribution.py 60 $numberOfSequences > $relativeOutput

# 3. Plot the relative distribution
echo "	Plotting the relative distribution."
cat $usableSequences | $scriptsPath/LengthDistribution.py 60 > lengthDist.tmp
Rscript $scriptsPath/PlotRelativeLengthDistribution.R $relativeOutput lengthDist.tmp $baseRef | sed '/null device/d' | sed -e '/          1/d'
