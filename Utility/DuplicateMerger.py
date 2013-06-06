#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-07-23

"""

This script filter a fastq file to keep only sequences with at least <minCopy> copy.
Saves data as fasta: All the identical sequences are collapsed together with the copy count in the header. 

Usage:

cat filename.fastq | DuplicateMerger.py <minCopy> > out.fasta
	minCopy: Minimal number of copy to keep sequence.

"""

VERSION="2.1"

class DuplicateMerger:
	def __init__(self, minCopy):
		self.minCopy = minCopy
		self.clear()

	def clear(self):
		self.seqCounts = {}

	def processSequence(self, sequence):
		if sequence == "":
			sys.stderr.write("Warning (QualityFilter.py): sequence line is empty.\n")
			return
		if sequence in self.seqCounts:
			self.seqCounts[sequence] += 1
		else:
			self.seqCounts[sequence] = 1

	def printResults(self):
		sorted(self.seqCounts, key = self.seqCounts.get)
		i = 0
		for seq in sorted(self.seqCounts, key = self.seqCounts.get, reverse = True):
			if self.seqCounts[seq] >= self.minCopy:
				sys.stdout.write("> PP_" + str(i) + " Count: " + str(self.seqCounts[seq]) + "\n")
				sys.stdout.write(seq + "\n")
#				print "> PP_" + str(i) + " Count: " + str(self.seqCounts[seq])
#				print seq
				i += 1

# Tests
sequence_tests_valid1="GGGCATCGATGCAGTCTATCGTAGTC"
sequence_tests_valid2="GGGGATCGATGCAGTCTATCGTAGTC"
sequence_tests_empty=""

def Tests():
        print "****"
	print "Valid - Add 1 sequence" 
	d = DuplicateMerger(1)
	d.processSequence(sequence_tests_valid1)
	print ""
	print "Expecting:"
	print "> PP_0 Count: 1"
	print "GGGCATCGATGCAGTCTATCGTAGTC"
	print ""
	print "Actual:"
	d.printResults()
	print ""
        print "****"
	print "Valid - Add 2 different sequences" 
	d = DuplicateMerger(1)
	d.processSequence(sequence_tests_valid1)
	d.processSequence(sequence_tests_valid2)
	print ""
	print "Expecting:"
	print "> PP_0 Count: 1"
	print "GGGCATCGATGCAGTCTATCGTAGTC"
	print "> PP_1 Count: 1"
	print "GGGGATCGATGCAGTCTATCGTAGTC"
	print ""
	print "Actual:"
	d.printResults()
	print ""
        print "****"
	print "Valid - Add 2 identical sequences" 
	d = DuplicateMerger(1)
	d.processSequence(sequence_tests_valid1)
	d.processSequence(sequence_tests_valid1)
	print ""
	print "Expecting:"
	print "> PP_0 Count: 2"
	print "GGGCATCGATGCAGTCTATCGTAGTC"
	print ""
	print "Actual:"
	d.printResults()
	print ""
        print "****"
	print "Valid - Below treshold"
	d = DuplicateMerger(2)
	d.processSequence(sequence_tests_valid1)
	d.processSequence(sequence_tests_valid2)
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	d.printResults()
	print ""
        print "****"
	print "Invalid - Empty sequence"
	d = DuplicateMerger(1)
	d.processSequence(sequence_tests_empty)
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	d.printResults()
	print ""
		

import sys

if __name__=="__main__":
	if len(sys.argv) == 1:
		print __doc__
                sys.exit(1)

	if sys.argv[1] == "version":
		print "QualityFilter.py v." + str(VERSION)
		sys.exit()

	if sys.argv[1] == "tests":
		Tests()
		sys.exit()

        if len(sys.argv)!=2:
                print __doc__
                sys.exit(1)

	try:
		minCopy = int(sys.argv[1])
	except ValueError:
		print ""
		print "Invalid minCopy(" + sys.argv[1] + "), must be an integer value."
		print __doc__
		sys.exit()

        duplicateMerger = DuplicateMerger(minCopy)

	sequence=""
	i = 0
	for line in sys.stdin:
		if i==1:
			sequence=line.strip()
			duplicateMerger.processSequence(sequence)

		i+=1
		if i==4:
			i=0

	duplicateMerger.printResults()
