#!/usr/bin/python
# encoding: utf-8
# author: Sébastien Boisvert
# 2012-01-24
# modified by: Charles Joly Beauparlant
# 2012-11-16

VERSION=2.1

def Usage():
	print ""
	print "This program keeps only things before the adaptor."
	print ""
	print "Usage "
	print "cat joe.fastq | ./AdaptorFinder.py <adaptorSequence> <trimSize> <maxMismatches> > joe.trimmed.fastq"
	print "trimSize: If no match is found with adaptor, the script will try again removing 1 base at a time,"
	print "          until a match is found or trimSize base are removed."
	print "maxMismatches: Maximum number of mismatches that can be tolerated."
	print ""

def mismatches(s1,s2):
	i=0
	mismatches=0
	while i<len(s1):
		if s1[i]!=s2[i]:
			mismatches+=1
		i+=1
	return mismatches

def findOffset(sequence,adaptor,maxMismatches):
	i=0
	sequenceLength=len(sequence)
	adaptorLength=len(adaptor)

	#no mismatch
	while i<=sequenceLength-adaptorLength:
		observed=sequence[i:adaptorLength]
		if observed==adaptor:
			return i

		i+=1

	i=0

	# mismatch
	while i<=sequenceLength-adaptorLength:
		observed=sequence[i:i+adaptorLength]

		nonHits=mismatches(observed,adaptor)
#		if nonHits<=adaptorLength/2:
		if nonHits<=maxMismatches:
			return i

		i+=1

	return sequenceLength


def process(header,sequence,dummy,quality,adaptor,maxMismatches):
	offset=findOffset(sequence,adaptor,maxMismatches)

	lengthDiff=len(sequence)-len(adaptor)

	# If count == 0, then there is only the adaptor
	# Sequence shorter than 10 nucleotides are removed
	if offset < 10:
		return True

	if offset>=10 and offset != len(sequence):
		print header
		print sequence[0:offset]
		print dummy
		print quality[0:offset]
		return True

	return False	

def processAdaptor(header,sequence,dummy,quality,adaptor,maxMismatches,trimSize):
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
	# If the full sequence of the adaptor was not found, we truncate it and try to find it again (up to trimSize tries)
	j=0
	tmpAdaptor = adaptor
	while j <= trimSize + 1: 
		if j <= trimSize:
			if process(header,sequence,dummy,quality,tmpAdaptor,maxMismatches):
				return
			j += 1
			tmpAdaptor=adaptor[j:len(adaptor)]
		else:
			print header
			print sequence
			print dummy
			print quality
			return

# Tests

adaptor_test="ACGTACTGGT"

header_test_normal="@HWI-EASXXX:2:1:0:672#0/1"
header_test_empty=""

dummy_test_normal="+HWI-EASXXX:2:1:0:672#0/1"
dummy_test_empty=""

quality_test_normal="DNVVVVVVVTSTTTTTPQSTTTTRNRTPRQ"
quality_test_empty=""

sequence_test_noAdaptor="ACTGCATGCATGCATCGATCGATACGGTGG"
sequence_test_adaptorEnd="ACTGCATGCATGCATCACGTACTGGTTCGT"
sequence_test_adaptor10="ACTGCATGCAACGTACTGGTCTCGTATTGC"
sequence_test_adaptor9= "ACTGCATGCACGTACTGGTCAGTCATCGTC"
sequence_test_adaptorStart="ACGTACTGGTCATGCTCTATCTCTACTATC"
sequence_test_adaptorTruncated3="ACTGCATGCATGCATCTACTGGTTCGTCTT"
sequence_test_adaptorMismatches2="ACTGCATGCATGCATCACTTACTGGTTCGT"
sequence_test_empty=""
sequence_test_short="NCGGGTGATGCGAACTGGAGTCTGAGCTT"

def Tests():
        print "****"
	print "Valid Entry No Adaptor"
	print ""
	print "Expecting:"
	print header_test_normal
	print sequence_test_noAdaptor
	print dummy_test_normal
	print quality_test_normal
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_noAdaptor, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
        print "****"
	print "Valid Entry Adaptor End"
	print ""
	print "Expecting:"
	print header_test_normal
	print "ACTGCATGCATGCATC"
	print dummy_test_normal
	print "DNVVVVVVVTSTTTTT"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptorEnd, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
        print "****"
	print "Valid Entry Adaptor 10"
	print ""
	print "Expecting:"
	print header_test_normal
	print "ACTGCATGCA"
	print dummy_test_normal
	print "DNVVVVVVVT"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptor10, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
        print "****"
	print "Valid Entry Adaptor 9"
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptor9, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
        print "****"
	print "Valid Entry Adaptor Start"
	print ""
	print "Expecting:"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptorStart, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
        print "****"
	print "Valid Entry Adaptor Truncated 3"
	print ""
	print "Expecting:"
	print header_test_normal
	print "ACTGCATGCATGCATC"
	print dummy_test_normal
	print "DNVVVVVVVTSTTTTT"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptorTruncated3, dummy_test_normal, quality_test_normal, adaptor_test, 0, 3)
	print ""
        print "****"
	print "Valid Entry Adaptor Mismatches 2"
	print ""
	print "Expecting:"
	print header_test_normal
	print "ACTGCATGCATGCATC"
	print dummy_test_normal
	print "DNVVVVVVVTSTTTTT"
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_adaptorMismatches2, dummy_test_normal, quality_test_normal, adaptor_test, 2, 0)
	print ""
	print "****"
	print "Invalid Entry - Empty header"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): header line is empty."
	print ""
	print "Actual:"
	processAdaptor(header_test_empty, sequence_test_noAdaptor, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
	print "****"
	print "Invalid Entry - Empty sequence"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): sequence line is empty."
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_empty, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""
	print "****"
	print "Invalid Entry - Empty dummy"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): dummy line is empty."
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_noAdaptor, dummy_test_empty, quality_test_normal, adaptor_test, 0, 0)
	print ""
	print "****"
	print "Invalid Entry - Empty quality"
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): quality line is empty."
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_noAdaptor, dummy_test_normal, quality_test_empty, adaptor_test, 0, 0)
	print ""
	print "****"
	print "Invalid Entry - Sequence and Quality does not match."
	print ""
	print "Expecting:"
	print "Warning (QualityFilter.py): sequence length does not match quality length."
	print ""
	print "Actual:"
	processAdaptor(header_test_normal, sequence_test_short, dummy_test_normal, quality_test_normal, adaptor_test, 0, 0)
	print ""

import sys

# Parse arguments
if len(sys.argv) == 1:
	Usage()
	sys.exit()

if sys.argv[1] == "version":
        print "AdaptorFinder.py v." + str(VERSION)
        sys.exit()

if sys.argv[1] == "usage":
        Usage()
        sys.exit()

if sys.argv[1] == "tests":
        Tests()
        sys.exit()

if len(sys.argv)!=4:
	Usage()
	sys.exit()

if sys.stdin.isatty():
        Usage()
        sys.exit()

adaptor=sys.argv[1]

if len(adaptor) < 10:
	print ""
	print "Length of adaptor must be longer than 10."
	print "Current adaptor length(" + adaptor + "): " + str(len(adaptor))
	Usage()
	sys.exit()

try:
	trimSize=int(sys.argv[2])
except ValueError:
        print ""
        print "Invalid trimSize: " + sys.argv[3]
        Usage()
        sys.exit()

try:
	maxMismatches=int(sys.argv[3])
except ValueError:
        print ""
        print "Invalid maxMismatches: " + sys.argv[3]
        Usage()
        sys.exit()

if (maxMismatches>len(adaptor)):
	print "maxMismatches value bigger than adaptor size."
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

		processAdaptor(l0,l1,l2,l3,adaptor,maxMismatches, trimSize)

	i+=1

	if i==4:
		i=0
