#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-27

"""
Parse nt blast output. Keep only sequences that do not have Homo Sapiens in highest score.
Usage:
nt_blast_analysis.py <blast.out> <output> <seq_count>
    blast.out: Output of blast against nt database.
    output: Prefix that will be added to output files.
    seq_count: The number of sequences after quality trimming.
"""

from Parser_BlastOutput import *

class BlastAnalyzer:
	def __init__(self, filename, output, total_count):
		self.parser = Parser(filename)
		self.output = output
		self.total_count = float(total_count)
		self.perfectSeqIDs = {}
		self.perfectCounts = {}
		self.perfectSequences = {}
		self.perfectFullName = {}

	def addSeqName(self, key, result, name):
		if result == "perfectMatch":
			if key in self.perfectSeqIDs:
				self.perfectSeqIDs[key].append(name)
			else:
				self.perfectSeqIDs[key] = []
				self.perfectSeqIDs[key].append(name)

	def addCount(self, key, result, count):
		if result == "perfectMatch":
			if key in self.perfectCounts:
				self.perfectCounts[key] += count
			else:
				self.perfectCounts[key] = count

	def addQuerySequence(self, key, result, sequence):
		if result == "perfectMatch":
			if key not in self.perfectSequences:
				self.perfectSequences[key] = sequence

	def addQueryFullName(self, key, result, fullName):
		if result == "perfectMatch":
			if key not in self.perfectFullName:
				self.perfectFullName[key] = fullName

	def getQuerySequence(self, key):
		return self.perfectSequences[key]

	def checkHuman(self, query_fullName):
		if "Homo" in query_fullName:
			return True
		if "Homo sapiens" in query_fullName:
			return True
		if "Sapiens" in query_fullName:
			return True
		if "sapiens" in query_fullName:
			return True
		if "Human" in query_fullName:
			return True
		if "human" in query_fullName:
			return True
		return False

	def processMatch(self, token, i):
		name = token.getNames()[i]
		query_fullName = token.getData(name, "query_fullName")
		query_start = token.getData(name, "query_start")
		query_end = token.getData(name, "query_end")
		query_length = token.getData("Query", "length")
		name_start = token.getData(name, "start")
		name_end = token.getData(name, "end")
		name_length = token.getData(name, "length")
		mismatches = token.getData(name, "mismatches")
		gaps = token.getData(name, "gaps")

		# Check if 
		if self.checkHuman(query_fullName) == True:
			return "noMatch, hsapiens"

		# Check for no match
		if query_start > 1 or query_end <= query_length - 1:
			return "noMatch, length"
		if gaps > 0:
			return "noMatch, gaps"
		if mismatches > 3:
			return "noMatch, mismatches"

		# Check for perfect match
		if mismatches == 0:
			return "perfectMatch"

		# Check for loose match
		else:
			return "looseMatch"

	def getShortName(self, name):
		return name.split()[0]

	def addResults(self, result, token, i):
		name = token.getNames()[i]
		sequence = token.getData(name, "query_sequence")
		fullName = token.getData(name, "query_fullName")
		count = int(token.getID().split()[2])
		self.addSeqName(name, result, token.getID())
		self.addQuerySequence(name, result, sequence)
		self.addQueryFullName(name, result, fullName)
		self.addCount(name, result, count)

	def processToken(self, token):
		# We keep only miRNA having the best score
		token.processScores()

		# Check if we have a perfect match, a loose match or no match at all
		for i in range (0, token.getNumberOfResult()):
			result = self.processMatch(token, i)
			if result == "perfectMatch":
				self.addResults(result, token, i)
				
	def parseFile(self):
		done = False
		count = 0
		while done == False:
			count += 1
			self.parser.createNextToken()
			if count % 1000 == 0:
				print count, " blast hits processed."
			if self.parser.isEOF() != True:
				token = self.parser.getToken()
				self.processToken(token)
			else:
				print count, " blast hits processed."
				done = True

	def printReport(self):
		filename = self.output + "_perfectMatches_summary.txt"
		f = open(filename, 'w')
		f.write("ID\tSequence\tcount\t%_of_total\n")
		for ID in self.perfectCounts:
			toPrint = self.perfectFullName[ID]
			toPrint = toPrint + "\t" + self.getQuerySequence(ID)
			toPrint = toPrint + "\t" + str(self.perfectCounts[ID])
			toPrint = toPrint + "\t" + str((self.perfectCounts[ID] / self.total_count) * 100)
			toPrint = toPrint + '\n'
			f.write(toPrint)
		f.close()

import sys

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print __doc__
		sys.exit(1)

	filename = sys.argv[1]
	output = sys.argv[2]
	total_count = sys.argv[3]
	blastAnalyzer = BlastAnalyzer(filename, output, total_count)
	blastAnalyzer.parseFile()
	blastAnalyzer.printReport()
