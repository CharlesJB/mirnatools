#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-10-10

"""
Usage:

removeDuplicatedMapping.py <file_list> <miRNA_fasta>

	file_list: The list of every fasta files for non-miRNA mappings, one file per line
"""


import sys

if __name__=="__main__":
        if len(sys.argv)!=3:
                print __doc__
                sys.exit(1)

	file_list = sys.argv[1]
	miRNA_fasta = sys.argv[2]

	# Calculate total count for miRNA
	total_miRNA = 0.0
	miRNA_name = miRNA_fasta.split("/")[len(miRNA_fasta.split("/"))-1].split('.')[0]
	miRNA_entry = {}
	for line in open(miRNA_fasta):
		if line[0] == '>':
			total_miRNA += float(line.split()[3])
			name = line.split()[1]
			miRNA_entry[name] = total_miRNA

	# Fetch every collapsed read and count number of file for each one of them
	read_count = {}
	file_count = {}
	for filename in open(file_list):
		for line in open(filename.strip()):
			if line[0] == '>':
				name = line.split()[1]
				if name not in file_count:
					file_count[name] = 0
				file_count[name] += 1
				if name not in miRNA_entry:
					count = float(line.split()[3])
					miRNA_entry[name] = count

	# Calculate normalized count for every rna specie
	rna_count = {}
	for filename in open(file_list):
		rna_name = filename.split("/")[len(filename.split("/"))-1].split('.')[0]
		total_count = 0.0
		for line in open(filename.strip()):
			if line[0] == '>':
				name = line.split()[1]
				total_count += miRNA_entry[name] / float(file_count[name])
		rna_count[rna_name] = total_count

	# Calcute percentage
	total_count = total_miRNA
	rna_percentage = {}
	for rna_name in rna_count:
		total_count += rna_count[rna_name]
	for rna_name in rna_count:
		rna_percentage[rna_name] = rna_count[rna_name] / total_count
	
	# Print results
	print "RNA_specie\tCount\tPercentage"
	print miRNA_name + "\t" + str(total_miRNA) + "\t" + str(total_miRNA/total_count)
	for rna_name in rna_count:
		print rna_name + "\t" + str(rna_count[rna_name]) + "\t" + str(rna_percentage[rna_name])
