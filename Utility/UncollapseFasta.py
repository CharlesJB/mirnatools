#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-06-21

"""
This script will create a unique fasta entry for every value of the count for each sequence.

Usage:

	cat <fasta_filename> | UncollapseFasta.py <prefix>

	fasta_filename: fasta file to uncollapse
	prefix: The prefix for the identifier of each entry produced

"""

import sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print __doc__
		sys.exit(1)

	prefix = sys.argv[1]

	print_count = 0
	sequence_count = 0
	sequence = ""

	for line in sys.stdin:
		if line[0] == '>':
			if sequence_count > 0 and len(sequence) > 0:
				for i in range(0, sequence_count):
					print ">" + prefix + "_" + str(print_count+1)
					print sequence
					print_count += 1
			sequence_count = int(line.strip().split()[len(line.strip().split())-1])
			sequence = ""
		else:
			sequence += line.strip()
