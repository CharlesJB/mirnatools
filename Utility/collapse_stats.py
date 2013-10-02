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

def print_line(name, strand, modification, value):
	t = "\t"
	print name + t + strand + t + modification + t + str(value)
		

miRNA_list = {}
print "name\tstrand\ttype\tpercent"
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

for entry in miRNA_list:
	total = miRNA_list[entry].none + miRNA_list[entry].other + miRNA_list[entry].uridylation + miRNA_list[entry].adenylation
	if total > 0.0:
		percent_none = miRNA_list[entry].none / total
		percent_other = miRNA_list[entry].other / total
		percent_uridylation = miRNA_list[entry].uridylation / total
		percent_adenylation = miRNA_list[entry].adenylation / total

		print_line(entry, miRNA_list[entry].strand, "none", percent_none)
		print_line(entry, miRNA_list[entry].strand, "other", percent_other)
		print_line(entry, miRNA_list[entry].strand, "uridylation", percent_uridylation)
		print_line(entry, miRNA_list[entry].strand, "adenylation", percent_adenylation)
#		print entry + t + strand + t + str(percent_none) + t + str(percent_other) + t + str(percent_uridylation) + t + str(percent_adenylation)
#		print entry + "\t" + miRNA_list[entry].strand + "\t" + str(miRNA_list[entry].none) + "\t" + str(miRNA_list[entry].other) + "\t" + str(miRNA_list[entry].uridylation) + "\t" + str(miRNA_list[entry].adenylation)
