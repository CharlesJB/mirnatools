#/bin/bash

Usage() {
	echo ""
	echo "Usage:"
	echo "GenomeMapping.sh <referenceData> <usableSequences> <scriptsPath>"
	echo "referenceData: For the blast analysis. Must be in fasta format."
	echo "usableSequences: List of valid RNA sequence (after trimming)."
	echo "scriptsPath: path containing the miRNA-Tools git clone."
	echo "			 (https://github.com/CharlesJB/miRNA-Tools)"
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
scriptsPath=$3

ValidateFile $reference
ValidateFile $usableSequences
ValidateFile $scriptsPath/BlastAnalysis.sh
ValidateFile $scriptsPath/LengthDistributionAnalysis.sh

# 1. Blast
echo ""
echo "Beginning blast analysis..."
$scriptsPath/BlastAnalysis.sh $reference $usableSequences "0.1"
echo "Done!"

blastOutput=$(basename ${reference%.*})".out"

echo ""
echo "Launching Preethi's algorithm..."
$scriptsPath/preethi_blast_analysis.py $blastOutput ${blastOutput%.*}
echo "Done!"

echo ""
echo "Beginning length distribution analysis..."
id=$(echo $(basename $usableSequences) | sed 's/_UsableSeq.fasta//g')
$scriptsPath/miRNA_LengthDistributionAnalysis.sh $id $scriptsPath
echo "Done!"

# 3. Contribution analysis
echo ""
echo "Beginning contribution analysis..."
baseRef=$(basename ${reference%.*})
$scriptsPath/miRNAContributionAnalysis.sh $id $usableSequences $scriptsPath
echo "Done!"

# 5. Clean up
rm -f *.tmp
echo ""
