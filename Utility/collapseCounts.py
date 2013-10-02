#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-09-25

import sys

toPrint = {}

for line in sys.stdin:
	value = line.split()[0]
	count = float(line.split()[1])
	if value not in toPrint:
		toPrint[value] = 0
	toPrint[value] += count

total = 0.0
for value in toPrint:
	total += toPrint[value]

for value in toPrint:
	print value + "\t" + str(toPrint[value])
#	print value + "\t" + "{0:.10f}".format(float(toPrint[value])/float(total))
