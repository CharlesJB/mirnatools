#!/bin/bash

Usage() {
	echo "Usage:"
	echo "pieChartProducter.sh <usableSequences> <miRNApath> <toGrep>"
	echo "    usableSequences: list of usable sequences obtained during raw sequence filtering"
	echo "    miRNApath: path to miRNA analysis (relative)"
	echo "    toGrep: extra filter for miRNA analysis files"
	echo "    scriptsPath: path containing the miRNA-Tools git clone."
	echo "		     ( https://github.com/CharlesJB/miRNA-Tools )"
}

usableSequences=$1
path=$2
toGrep=$3
scriptsPath=$4

mkdir -p ExpressionLevels

file=$path/$(ls $path | grep mature | grep combined | grep $toGrep | head -n1)

count=$(grep '>' $usableSequences | cut -f4 -d' ' | awk '{ sum += $1 } END {print sum/1000000}')
toAwk="awk '{ printf \"%4.3f\\n\", (\$2/$count) }' $file > ExpressionLevels/percent.tmp"
eval $toAwk
cut -f1 $file > ExpressionLevels/names.tmp
output="ExpressionLevels/mature_combinedMatches_relativeCount_"$toGrep".txt"
paste ExpressionLevels/names.tmp ExpressionLevels/percent.tmp | sort -k2 -nr > $output
rm ExpressionLevels/*.tmp

distributionFile=${output%.*}"_distribution.txt"
sort -k2 -nr $output | head -n15 > $distributionFile
echo -e "others\t$(sort -k2 -nr $output | tail -n+16 | awk '{ sum += $2 } END { print sum}')" >> $distributionFile

pieChartBasename=$toGrep"_pieChart
Rscript $scriptsPath/percent_pie_chart.R jpeg $distributionFile $pieChartBasename | sed '/null device/d' | sed -e '/          1/d'
Rscript $scriptsPath/percent_pie_chart.R tiff $distributionFile $pieChartBasename | sed '/null device/d' | sed -e '/          1/d'
