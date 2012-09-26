#/bin/bash

Usage() {
	echo ""
	echo "Usage:"
	echo "GenomeMapping.sh <referenceData> <usableSequences> <scriptsPath>"
	echo "referenceData: For the blast analysis. Must be in fasta format."
	echo "usableSequences: List of valid RNA sequence (after trimming)."
	echo "scriptsPath: path containing the miRNA-Tools git clone."
	echo "             (https://github.com/CharlesJB/miRNA-Tools)"
}

ValidateFile() {
	fileToCheck=$1
	if [ ! -e "$fileToCheck" ]
	then
		echo "Invalid file: $fileToCheck"
		Usage
		exit
	fi
}

reference=$1
usableSequences=$2
idToRemove=$3
scriptsPath=$4

ValidateFile $reference
ValidateFile $usableSequences
ValidateFile $idToRemove
ValidateFile $scriptsPath/BlastAnalysis.sh
ValidateFile $scriptsPath/LengthDistributionAnalysis.sh

# 1. Blast
echo ""
echo "Beginning blast analysis..."
$scriptsPath/BlastAnalysis.sh $reference $usableSequences
echo "Done!"

echo ""
echo "Beginning length distribution analysis..."
blastOutput="$(basename ${reference%.*}).out"
$scriptsPath/LengthDistributionAnalysis.sh $blastOutput $idToRemove $usableSequences $scriptsPath
echo "Done!"

# 3. Contribution analysis
echo ""
echo "Beginning contribution analysis..."
baseRef=$(basename ${reference%.*})
$scriptsPath/ContributionAnalysis.sh $baseRef $usableSequences $scriptsPath
echo "Done!"

# 4. Calculation of the total of sequences
echo ""
echo "Calculating total number of sequences..."
totalOutput=$baseRef"_totalCount.txt"
grep '>' $baseRef"_blast_ID_unique.fa" | awk '{ SUM += $4 } END { print SUM }' > $totalOutput
totalCount=$(cat $totalOutput)
echo "    $totalCount sequences matched $baseRef!"
echo "Done!"

# 5. Clean up
rm *.tmp
echo ""
