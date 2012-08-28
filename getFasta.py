#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-28

"""
Return a subset of a fasta file containing only specifics ID
usage:
getFasta.py <ID_filename> <fasta_filename> 
ID_filename: File containing the list of ID, one ID per file.
fasta_filename: fasta file containing ID and count in header.
"""
class CountGetter:
	def __init__(self, ID_filename):
		self.IDs = []
		self.fetchIDs(ID_filename)

	def fetchIDs(self, ID_filename):
		for line in open(ID_filename):
			tokens = line.split()
			self.IDs.append(tokens[0])

	def fetchID(self, line):
		return line.split()[0][1:]

	def printAll(self, fasta_filename):
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
				if ID in self.IDs:
					print header
					print sequence
				count = 0

import sys

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(1)

	ID_filename = sys.argv[1]
	fasta_filename = sys.argv[2]
		
	countGetter = CountGetter(ID_filename)
	countGetter.printAll(fasta_filename)
