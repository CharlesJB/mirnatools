#!/bin/bash
id=$1
output=$id".wiki"
link=$2

# 1. Analysis name
echo "= $id - Analysis Report =" > $output

# 2. Raw Sequence Filtering
echo "== Raw Sequence Filtering ==" >> $output
folder="RawSequenceFiltering"
count_original=$(wc -l $folder/$(ls $folder | grep $id | grep txt | grep -v 'trimmed\|qual' | head -n1) | awk '{ print $1/4 }')
echo "Original count: "$count_original" sequences.<br/>" >> $output
count_qual=$(wc -l $folder/$(ls $folder | grep $id | grep qual | grep txt | head -n1) | awk '{ print $1/4 }')
percent_qual=$(echo "scale=2; $count_qual / $count_original * 100" | bc -l)
echo "After Trimming: "$count_qual" sequences ("$percent_qual" of original count).<br/>" >> $output
unique_count=$(grep '>' $folder/$(ls $folder | grep $id | grep fasta) | wc -l | awk '{ print $1 }')
echo "<b>Usable sequences: "$unique_count" unique sequences.</b><br/>" >> $output
echo "" >> $output

# 3. Length Distribution Analysis
echo "== Length Distribution Analysis ==" >> $output
figure_name=$(ls LengthDistributionAnalysis | grep $id | grep jpeg | head -n1)
echo $link"/"$figure_name >> $output
echo "" >> $output

# 4. Genome mapping - miRNAs
echo "== Genome Mapping - miRNAs ==" >> $output

#   4.1 Counts
echo "=== Sequence Count ===" >> $output
folder="GenomeMapping/miRNA"
perfect_count=$(wc -l $folder/$(ls $folder | grep $id | grep perfectMatches | head -n1) | awk '{ print $1 }')
loose_count=$(wc -l $folder/$(ls $folder| grep $id | grep looseMatches | head -n1) | awk '{ print $1 }')
echo "After using Preethi's algorithm, we have:<br/>" >> $output
echo "*"$perfect_count" miRNAs that were perfect matches<br/>" >> $output
echo "*"$loose_count" miRNAs that were loose matches<br/>" >> $output
echo "" >> $output

noMatchGaps_count=$(wc -l $folder/$(ls $folder| grep -v length | grep $id | grep 'noMatchGaps.txt' | head -n1) | awk '{ print $1 }')
noMatchLength_count=$(wc -l $folder/$(ls $folder | grep -v length | grep $id | grep 'noMatchLength.txt' | head -n1) | awk '{ print $1 }')
noMatchMismatches_count=$(wc -l $folder/$(ls $folder | grep -v length | grep $id | grep 'noMatchMismatches.txt' | head -n1) | awk '{ print $1 }')
echo "Details of the sequences that were rejected:<br/>" >> $output
echo "*"$noMatchGaps_count" miRNAs were rejected because of gaps<br/>" >> $output
echo "*"$noMatchLength_count" miRNAs were rejected because of length<br/>" >> $output
echo "*"$noMatchMismatches_count" miRNAs were rejected because of mismatches<br/>" >> $output
echo "" >> $output

#   4.2 Insert figures
echo "=== Length Distribution and Relative Contribution ===" >> $output
for file in $(ls GenomeMapping/miRNA | grep $id | grep jpeg | grep -v relative)
do
	ID=$(echo $file | sed 's/lengthDist_mature_//g' | sed 's/_ID\.jpeg//g')
	relative="lengthDist_"$ID"_relative.jpeg"
	fasta="mature_"$ID".fa"
	count=$(grep '>' "GenomeMapping/miRNA/"$fasta | awk '{ sum += $4 } END { print sum }')
	if [ "$count" == "" ]
	then
		count=0
	fi
	echo "==== "$ID" ====" >> $output
	echo "Number of sequences: "$count"<br/>">> $output
	file="miRNA_"$file
	relative="miRNA_"$relative
	echo $link"/"$file >> $output
	echo $link"/"$relative >> $output
	echo "" >> $output
done

echo "=== Expression Levels ===" >> $output
pieChart=$(ls ExpressionLevels | grep $id | grep jpeg | grep pie)
expression=$(ls ExpressionLevels | grep $id | grep jpeg | grep expressionLevels)
echo $link"/"$pieChart >> $output
echo $link"/"$expression >> $output

# 5. Genome mapping - other types of small RNAs
echo "== Genome mapping - other types of small RNAs ==" >> $output
dirpath="GenomeMapping/"$id
for dir in $(ls $dirpath | grep -v old)
do
	length=$id"_"$dir"lengthDist_"$dir".jpeg"
	relative=$id"_"$dir"lengthDist_"$dir"_relative.jpeg"
	count=$(cat $dirpath/$dir/$dir"_totalCount.txt")
	if [ "$count" == "" ]
	then
		count=0
	fi
	echo "=== "$dir" ===" >> $output
	echo "Number of sequences: "$count"<br/>" >> $output
	echo "$link/$length" >> $output
	echo "$link/$relative" >> $output
	echo "" >> $output
done
