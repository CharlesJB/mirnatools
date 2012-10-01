#!/bin/bash

Usage() {
	echo "Usage:"
	echo "LengthDistribtionAnalysis.sh <blastOutput> <idToRemove> <usableSequences> <scriptsPath>"
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

blastOutput=$1
idToRemove=$2
usableSequences=$3
scriptsPath=$4

# File validation
ValidateFile $blastOutput
ValidateFile $idToRemove
ValidateFile $usableSequences
ValidateFile $scriptsPath/getFasta.py
ValidateFile $scriptsPath/LengthDistributionAnalysis.sh

# 1. Extract sequence IDs
echo "	Extracting sequence IDs."
blastIdOutput=${blastOutput%.*}_blast_ID.txt
grep 'Sequences producing' $blastOutput -B4 | grep 'Query=' | cut -f2 -d' ' \
| sed 's/PP_//g' | sort -n | sed 's/^/PP_/g' > $blastIdOutput

# 2. Remove miRNA sequence IDs
echo "	Remove known miRNA IDs."
uniqueOutput=${blastIdOutput%.*}'_unique.txt'
sort $blastIdOutput > out.tmp
sort $idToRemove > toRemove.tmp
comm -23 out.tmp toRemove.tmp | sed 's/PP_//g' | sort -n | sed 's/^/PP_/g' > $uniqueOutput
rm *.tmp

# 3. Fetch fasta sequences corresponding to IDs
echo "	Converting IDs to fasta."
uniqueFastaOutput=${uniqueOutput%.*}".fa"
$scriptsPath/getFasta.py $uniqueOutput $usableSequences > $uniqueFastaOutput

# 4. Calculate length distribution
echo "	Calculating length distribution."
lengthDistOutput="lengthDist_$uniqueOutput"
cat $uniqueFastaOutput | $scriptsPath/LengthDistribution.py 60 > $lengthDistOutput

# 5. Plot the distribution
echo "	Plotting the distribution."
Rscript $scriptsPath/PlotLengthDistribution.R $lengthDistOutput ${blastOutput%.*} | sed '/null device/d' | sed -e '/          1/d'
