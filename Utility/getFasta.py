#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-28

"""
Return a subset of a fasta file containing only specifics ID

Usage:

	cat <ID_filename> | getFasta.py <fasta_filename> <match_status>

	ID_filename: ID file from BlastAnalysis
	fasta_filename: fasta file containing ID and count in header.
	match_status: match or not-match

"""
class CountGetter:
	def __init__(self):
		self.IDs = []
		self.fetchIDs()

	def fetchIDs(self):
		for line in sys.stdin:
#			for token in line.strip().split()[2:]:
			for token in line.strip().split():
				if token not in self.IDs:
					self.IDs.append(token)

	def fetchID(self, line):
#		print line.split()
		return line.split()[1]

	def printAll(self, fasta_filename, match_status):
		count = 0
		header = ""
		sequence = ""
		for line in open(fasta_filename):
			if count == 0:
				header = line.strip()
				count = 1
			else:
				sequence = line.strip()
				ID = self.fetchID(header)
				if match_status == "match":
#					print "ID: " + ID
					if ID in self.IDs:
						print header
						print sequence
				else:
					if ID not in self.IDs:
						print header
						print sequence
				count = 0

import sys

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(1)

	fasta_filename = sys.argv[1]
	match_status = sys.argv[2]

	countGetter = CountGetter()
	countGetter.printAll(fasta_filename, match_status)
