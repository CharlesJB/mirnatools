#!/usr/bin/python
# encoding: utf-8
# Author: Charles Joly Beauparlant
# 2013-06-14

VERSION=0.1

def Usage():
	print ""
	print "This program will convert U to T in sequences"
	print ""
	print "Usage "
	print "cat joe.fasta | ./convertRNAfastaToDNA.py > output.fasta"
	print ""

import sys

for line in sys.stdin:
	if line[0] != '>':
		line = line.upper().replace('U', 'T')
	print line.strip()
