#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-06-21

"""

This script will create summary files for expression analysis.

Usage:

cat filename.fastq | BlastSummaryAnalysis.py 

"""

import sys

names = []
count = {}
total = 0.0

for line in sys.stdin:
	if line.split()[0] == "perfect_match" or line.split()[0] == "loose_match":
		name = line.strip().split()[1]
		if name not in names:
			names.append(name)
			count[name] = float(line.strip().split()[2])
			total += count[name]

for name in names:
	print name + "\t" + str('{0:.10f}'.format((count[name]/total)*100.0))
