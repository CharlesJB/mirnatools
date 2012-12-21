#!/bin/bash

Usage() {
	echo "Usage:"
	echo "ExpressionAnalysis.sh <miRNApath> <toGrep> <scriptsPath>"
	echo "    miRNApath: path to miRNA analysis (relative)"
	echo "    toGrep: extra filter for miRNA analysis files"
	echo "    scriptsPath: path containing the miRNA-Tools git clone."
	echo "		     ( https://github.com/CharlesJB/miRNA-Tools )"
}

path=$1
toGrep=$2
scriptsPath=$3

mkdir -p ExpressionLevels

# Get the sum of all counts from GenomeMapping miRNA analysis
file=$path/$(ls $path | grep mature | grep combined | grep $toGrep | head -n1)
count=$(awk '{ sum += $2 } END {printf "%2.2f\n", sum}' $file)

# Convert the absolute sequence count by the sum of all absolute count
toAwk="awk '{ printf \"%4.8f\\n\", (\$2/$count) }' $file > ExpressionLevels/percent.tmp"
eval $toAwk

# Merge the percentage with the correct names and clean up intermediate files
cut -f1 $file > ExpressionLevels/names.tmp
output="ExpressionLevels/mature_combinedMatches_relativeCount_"$toGrep".txt"
paste ExpressionLevels/names.tmp ExpressionLevels/percent.tmp | sort -k2 -nr > $output
rm ExpressionLevels/*.tmp

# Generate the data file to be used by the graph production scripts
distributionFile=${output%.*}"_distribution.txt"
#   pie chart graph
sort -k2 -nr $output | head -n15 > $distributionFile
count=$(wc -l $distributionFile | awk '{print $1}')
if [ "$count" -ge 15 ]
then
	echo -e "others\t$(sort -k2 -nr $output | tail -n+16 | awk '{ sum += $2 } END { print sum}')" >> $distributionFile
fi

pieChartBasename="ExpressionLevels/"$toGrep"_pieChart"
Rscript $scriptsPath/Rscripts/percent_pie_chart.R jpeg $distributionFile $pieChartBasename | sed '/null device/d' | sed -e '/          1/d'
Rscript $scriptsPath/Rscripts/percent_pie_chart.R tiff $distributionFile $pieChartBasename | sed '/null device/d' | sed -e '/          1/d'

# miRNA expression level graph
absoluteCountFile="ExpressionLevels/"$(basename $file)
expressionBasename="ExpressionLevels/"$toGrep"_expressionLevels"
sort -k2 -nr $file > $absoluteCountFile
Rscript $scriptsPath/Rscripts/profiling.R jpeg $absoluteCountFile $expressionBasename | sed '/null device/d' | sed -e '/          1/d' 
Rscript $scriptsPath/Rscripts/profiling.R tiff $absoluteCountFile $expressionBasename | sed '/null device/d' | sed -e '/          1/d'
