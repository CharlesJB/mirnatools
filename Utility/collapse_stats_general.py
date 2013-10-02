#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-09-30

import sys

class miRNA:
	def __init__(self):
		self.reset()

	def reset(self):
		self.strand = ""
		self.none = 0.0
		self.other = 0.0
		self.uridylation = 0.0
		self.adenylation = 0.0

#def print_line(name, strand, modification, value):
#	t = "\t"
#	print name + t + strand + t + modification + t + str(value)
#		

toPrint = sys.argv[1]
if toPrint == "":
	toPrint = "all"

miRNA_list = {}
print "strand\ttype\tcount"
print "-\turidylation_predominant\t0"
print "-\tadenylation_predominant\t0"
print "-\tsimilar\t0"
print "-\texclusive_uridylation\t0"
print "-\texclusive_adenylation\t0"
#print "name\tstrand\tnone\tother\turidylation\tadenylation"
for line in sys.stdin:
	name = line.split()[1]
	strand = line.split()[2]
	modification = line.split()[3]
	count = float(line.split()[4])

	if name not in miRNA_list:
		miRNA_list[name] = miRNA()

	miRNA_list[name].strand = strand

	if modification== "none":
		miRNA_list[name].none += count
	elif modification == "other":
		miRNA_list[name].other+= count
	elif modification == "uridylation":
		miRNA_list[name].uridylation += count
	elif modification == "adenylation":
		miRNA_list[name].adenylation += count

grand_total = 0.0
total_3p = 0.0
total_5p = 0.0

for entry in miRNA_list:
	total = miRNA_list[entry].none + miRNA_list[entry].other + miRNA_list[entry].uridylation + miRNA_list[entry].adenylation
	if total > 0.0:
		strand = miRNA_list[entry].strand
		grand_total += total
		if strand == "3p":
			total_3p += total
		if strand == "5p":
			total_5p += total

for entry in miRNA_list:
	total = miRNA_list[entry].none + miRNA_list[entry].other + miRNA_list[entry].uridylation + miRNA_list[entry].adenylation
	if total > 0.0:
		strand = miRNA_list[entry].strand
		if toPrint == "all" or strand == toPrint:
			none = miRNA_list[entry].none
			other = miRNA_list[entry].other
			uridylation = miRNA_list[entry].uridylation
			adenylation = miRNA_list[entry].adenylation

			result = ""
			if adenylation >= 0.1 and uridylation >= 0.1:
				if adenylation > uridylation*2:
					result="adenylation_predominant"
				elif uridylation > adenylation*2:
					result="uridylation_predominant"
				else:
					retult="similar"
			elif adenylation >= 0.1 and uridylation < 0.1:
				result="exclusive_adenylation"
			elif uridylation >= 0.1 and adenylation< 0.1:
				result="exclusive_uridylation"
			elif none >= 0.1:
				result="other"
			else:
				result="none"

			if result != "none":
				total_print = 0.0
				if toPrint == "all":
					total_print = total / grand_total
				elif toPrint == "3p":
					total_print = total / total_3p 
				elif toPrint == "5p":
					total_print = total / total_5p 
				t = "\t"
				print strand + t + result + t + str(total_print)
#		print_line(entry, miRNA_list[entry].strand, "none", percent_none)
#		print_line(entry, miRNA_list[entry].strand, "other", percent_other)
#		print_line(entry, miRNA_list[entry].strand, "uridylation", percent_uridylation)
#		print_line(entry, miRNA_list[entry].strand, "adenylation", percent_adenylation)
