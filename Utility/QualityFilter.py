#!/usr/bin/python
# encoding: utf-8
# author: Sébastien Boisvert
# 2012-01-24
# modified by: Charles Joly Beauparlant

import sys
VERSION=2.1

"""

This program trims sequences based on their Phred Qualilty Score.
The longest stretch of sequence with score higher or equal to <minQual> will be printed.

Usage: 
	cat joe.fastq | ./AdaptorRemover.py <minQual> > joe.trimmed.fastq

"""

def process(header,sequence,dummy,quality,minQual):
	if header == "":
		sys.stderr.write("Warning (QualityFilter.py): header line is empty.\n")
		return

	if sequence == "":
		sys.stderr.write("Warning (QualityFilter.py): sequence line is empty.\n")
		return

	if dummy == "":
		sys.stderr.write("Warning (QualityFilter.py): dummy line is empty.\n")
		return

	if quality == "":
		sys.stderr.write("Warning (QualityFilter.py): quality line is empty.\n")
		return

	if len(sequence) != len(quality):
		sys.stderr.write("Warning (QualityFilter.py): sequence length does not match quality length.\n")
		return

	current=0
	longest=0
	longestLength=0
	# If it's not possible to get a better score, stop looping
	while len(sequence)-longestLength > current:
		miss=False
		i=0

		 # Stop checking when you find a score below treshold or if you are at the end of the sequence
		while miss==False and current+i < len(sequence):
			if ord(quality[current+i]) < minQual+33:
				miss=True
			else:
				i=i+1

		if i > longestLength:
			longest=current
			longestLength=i

		current=current+i+1

	# Only print if you have a sequence of at least 10 nucleotides with high enough quality score
	if longestLength >= 10:
		print header
		print sequence[longest:longest+longestLength]
		print dummy
		print quality[longest:longest+longestLength]

# Tests 

header_test_normal="@HWI-EASXXX:2:1:0:672#0/1"
header_test_empty=""
sequence_test_normal="NCGGGTGATGCGAACTGGAGTCTGAGCATCTCGTAT"
sequence_test_empty=""
sequence_test_short="NCGGGTGATGCGAACTGGAGTCTGAGCATCTCGTA"
dummy_test_normal="+HWI-EASXXX:2:1:0:672#0/1"
dummy_test_empty=""
quality_test_normal="DNVVVVVVVTSTTTTTPQSTTTTRNRTPRQQRQBBB"
quality_test_10valid="DNVVVVVVVT*TTTTTP*STTTTR*RTPRQQ*QBBB"
quality_test_10validMiddle="DNVVV*VVVTVTTTTT**STTTTR*RTPRQQ*QBBB"
quality_test_9valid="DNVVVVVVV**TTTTTP*STTTTR*RTPRQQ*QBBB"
quality_test_0valid="************************************"
quality_test_empty=""

def Tests():
        print "****"
	print "Valid Entry - The complete sequence is valid"
	print ""
	print "Expecting:"
	print header_test_normal
	print sequence_test_normal
	print dummy_test_normal
	print quality_test_normal
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_normal, 10)
	print ""
        print "****"
	print "Valid Entry - 10 nucleotides of the sequence are valid"
	print ""
	print "Expecting:"
	print header_test_normal
	print "NCGGGTGATG"
	print dummy_test_normal
	print "DNVVVVVVVT"
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_10valid, 10)
	print ""
        print "****"
	print "Valid Entry - 10 nucleotides of the sequence are valid and in middle of sequence"
	print ""
	print "Expecting:"
	print header_test_normal
	print "GATGCGAACT"
	print dummy_test_normal
	print "VVVTVTTTTT"
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_10validMiddle, 10)
	print ""
        print "****"
	print "Valid Entry - 9 nucleotides of the sequence are valid"
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_9valid, 10)
	print ""
	print "****"
	print "Valid Entry - None of the sequence is valid"
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_0valid, 10)
	print ""
	print "****"
	print "Valid Entry - minQual is 0"
	print ""
	print "Expecting:"
	print header_test_normal
	print sequence_test_normal
	print dummy_test_normal
	print quality_test_normal
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_normal, 0)
	print ""
	print "****"
	print "Valid Entry - minQual smaller than 0"
	print ""
	print "Expecting:"
	print header_test_normal
	print sequence_test_normal
	print dummy_test_normal
	print quality_test_normal
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_normal, -1)
	print ""
	print "****"
	print "Invalid Entry - Empty header"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): header line is empty."
	print ""
	print "Actual:"
	process(header_test_empty, sequence_test_normal, dummy_test_normal, quality_test_normal, 10)
	print ""
	print "****"
	print "Invalid Entry - Empty sequence"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): sequence line is empty."
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_empty, dummy_test_normal, quality_test_normal, 10)
	print ""
	print "****"
	print "Invalid Entry - Empty dummy"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): dummy line is empty."
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_empty, quality_test_normal, 10)
	print ""
	print "****"
	print "Invalid Entry - Empty quality"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): quality line is empty."
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_normal, dummy_test_normal, quality_test_empty, 10)
	print ""
	print "****"
	print "Invalid Entry - Sequence and Quality does not match."
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): sequence length does not match quality length."
	print ""
	print "Actual:"
	process(header_test_normal, sequence_test_short, dummy_test_normal, quality_test_normal, 10)
	print ""

import sys

# Parse arguments
if len(sys.argv) == 1:
	print __doc__
	sys.exit()

if sys.argv[1] == "version":
        print "QualityFilter.py v." + str(VERSION)
        sys.exit()

if sys.argv[1] == "usage":
        print __doc__
        sys.exit()

if sys.argv[1] == "tests":
        Tests()
        sys.exit()

if len(sys.argv)!=2:
        print __doc__
        sys.exit()

if sys.stdin.isatty():
        print __doc__
        sys.exit()

try:
        minQual=int(sys.argv[1])
except ValueError:
        print ""
        print "Invalid minQual: " + sys.argv[1]
        print __doc__
        sys.exit()

i=0

l0=""
l1=""
l2=""
l3=""

for line in sys.stdin:
	if i==0:
		l0=line.strip()
	elif i==1:
		l1=line.strip()
	elif i==2:
		l2=line.strip()
	elif i==3:
		l3=line.strip()

		j=0

		process(l0,l1,l2,l3,minQual)

	i+=1

	if i==4:
		i=0
