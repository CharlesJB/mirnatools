#/bin/bash

Usage() {
	echo ""
	echo "Usage:"
	echo "GenomeMapping.sh <referenceData> <usableSequences>"
	echo "referenceData: For the blast analysis. Must be in fasta format."
	echo "usableSequences: List of valid RNA sequence (after trimming)."
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
reference=$1
usableSequences=$2

if [ $# -eq $EXPECTED_ARGS ]
then
	# Validate files
	ValidateFile $reference
	ValidateFile $usableSequences
	ValidateFile $MIRNA_TOOLS_PATH/BlastAnalysis.sh
	ValidateFile $MIRNA_TOOLS_PATH/miRNA_LengthDistributionAnalysis.sh
	ValidateFile $MIRNA_TOOLS_PATH/miRNAContributionAnalysis.sh

	# Do the actual analysis
	echo ""
	echo "Beginning blast analysis..."
	BlastAnalysis.sh $reference $usableSequences "0.1"
	echo "Done!"

	blastOutput=$(basename ${reference%.*})".out"

	echo ""
	echo "Launching Preethi's algorithm..."
	preethi_blast_analysis.py $blastOutput ${blastOutput%.*}
	echo "Done!"

	echo ""
	echo "Beginning length distribution analysis..."
	id=$(echo $(basename $usableSequences) | sed 's/_UsableSeq.fasta//g')
	miRNA_LengthDistributionAnalysis.sh $id $scriptsPath
	echo "Done!"

	echo ""
	echo "Beginning contribution analysis..."
	baseRef=$(basename ${reference%.*})
	miRNAContributionAnalysis.sh $id $usableSequences $scriptsPath
	echo "Done!"

	# Clean up
	rm -f *.tmp
	echo ""
else
	Usage
fi
