#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-27

"""

Based on the algorithm presented in: Creighton et al 2009

Usage:
    BlastAnalysis.py <blast.out> <output>

    blast.out: Output of blast against miRBase (tested with version 19).
    output: summary or id

"""

from BlastParser import *

class BlastAnalyzer:
	def __init__(self, filename):
		self.parser = Parser(filename)
		self.output = output
		self.perfectCounts = {}
		self.perfectSpecies = {}
		self.perfectSequences = {}
		self.perfectSeqIDs = {} # The name of the sequences that blast to each miRNA
		self.looseCounts = {}
		self.looseSeqIDs = {} # The name of the sequences that blast to each miRNA
		self.combinedCounts = {}
		# Below are for the analysis of sequences that didn't pass the filtering
		self.noMatchLength = {}
		self.noMatchMismatches = {}
		self.noMatchGaps = {}
		self.noMatchLength_IDs = {}
		self.noMatchMismatches_IDs = {}
		self.noMatchGaps_IDs = {}

	def addSeqName(self, key, result, name):
		if result == "perfectMatch":
			if key in self.perfectSeqIDs:
				self.perfectSeqIDs[key].append(name)
			else:
				self.perfectSeqIDs[key] = []
				self.perfectSeqIDs[key].append(name)
		if result == "looseMatch":
			if key in self.looseSeqIDs:
				self.looseSeqIDs[key].append(name)
			else:
				self.looseSeqIDs[key] = []
				self.looseSeqIDs[key].append(name)

	def addCount(self, key, result, count):
		if result == "perfectMatch":
			if key in self.perfectCounts:
				self.perfectCounts[key] += count
			else:
				self.perfectCounts[key] = count
		if result == "looseMatch":
			if key in self.looseCounts:
				self.looseCounts[key] += count
			else:
				self.looseCounts[key] = count
		if key in self.combinedCounts:
			self.combinedCounts[key] += count
		else:
			self.combinedCounts[key] = count

	def addSpecieName(self, key, result, specie):
		if result == "perfectMatch":
			if key not in self.perfectSpecies:
				self.perfectSpecies[key] = specie

	def addQuerySequence(self, key, result, sequence):
		if result == "perfectMatch":
			if key not in self.perfectSequences:
				self.perfectSequences[key] = sequence

	def getSpecieName(self, key):
		return self.perfectSpecies[key]

	def getQuerySequence(self, key):
		return self.perfectSequences[key]

	def getCount(self, key, result):
		if result == "perfectMatch":
			return self.perfectCounts[key]
		if result == "looseMatch":
			return self.looseCounts[key]

	def processMatch(self, token, i):
		name = token.getNames()[i]
		query_start = token.getData(name, "query_start")
		query_end = token.getData(name, "query_end")
		query_length = token.getData("Query", "length")
		name_start = token.getData(name, "start")
		name_end = token.getData(name, "end")
		name_length = token.getData(name, "length")
		mismatches = token.getData(name, "mismatches")
		gaps = token.getData(name, "gaps")

		# Check for no match
		if query_start > 4 or query_end <= query_length - 4:
			return "noMatch, length"
		if name_start > 4 or name_end <= name_length - 4:
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

	def getSpecie(self, fullName):
		tokens = fullName.split()
		return ' '.join(tokens[2:len(tokens)-1])

	def addResults(self, result, token, i):
		count = float(token.getCount()) / float(token.getNumberOfResult())
		name = token.getNames()[i]
		fullName = token.getData(name, "query_fullName")
		specie = self.getSpecie(fullName)
		sequence = token.getData(name, "query_sequence")
		self.addCount(name, result, count)
		self.addSeqName(name, result, token.getID())
		self.addSpecieName(name, result, specie)
		self.addQuerySequence(name, result, sequence)

	def addNoMatch(self, container, containerID, token, i):
		count = float(token.getCount()) / float(token.getNumberOfResult())
		name = self.getShortName(token.getNames()[i])
		ID = token.getID()
		if name in container:
			container[name] += count
		else:
			container[name] = count			
		if name not in containerID:
			containerID[name] = []
		containerID[name].append(ID)

	def processToken(self, token):
		# We keep only miRNA having the best score
		token.processScores()

		# Check if we have a perfect match, a loose match or no match at all
		for i in range (0, token.getNumberOfResult()):
			result = self.processMatch(token, i)
			if result == "looseMatch" or result == "perfectMatch":
				self.addResults(result, token, i)
			else:
				if result == "noMatch, length":
					self.addNoMatch(self.noMatchLength, self.noMatchLength_IDs, token, i)
				elif result == "noMatch, gaps":
					self.addNoMatch(self.noMatchGaps, self.noMatchGaps_IDs, token, i)
				elif result == "noMatch, mismatches":
					self.addNoMatch(self.noMatchMismatches, self.noMatchMismatches_IDs, token, i)
				
	def parseFile(self):
		done = False
		count = 0
		while done == False:
			count += 1
			self.parser.createNextToken()
			if count % 1000 == 0:
				sys.stderr.write(str(count) + " blast hits processed.\n")
			if self.parser.isEOF() != True:
				token = self.parser.token
				self.processToken(token)
			else:
				sys.stderr.write(str(count) + " blast hits processed.\n")
				done = True

	def printCount(self, comment, container):
		for miRNA in container:
			print comment + '\t' + miRNA + '\t' + str(container[miRNA])

	def printID(self, comment, container):
		for miRNA in container:
			toPrint = comment + '\t' + miRNA
			for i in range(0, len(container[miRNA])):
				toPrint = toPrint + '\t' + container[miRNA][i]
			print toPrint

	def printSummary(self):
		# Print perfect matches
		self.printCount("perfect_match", self.perfectCounts)
		# Print loose matches
		self.printCount("loose_match", self.looseCounts)
		# Print no matches
		self.printCount("no_match_length", self.noMatchLength)
		self.printCount("no_match_gap", self.noMatchGaps)
		self.printCount("no_match_mismatch", self.noMatchMismatches)

	def printAllIDs(self):
		self.printID("perfect_match", self.perfectSeqIDs)
		self.printID("loose_match", self.looseSeqIDs)

import sys

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(1)

	filename = sys.argv[1]
	output = sys.argv[2]
	blastAnalyzer = BlastAnalyzer(filename)
	blastAnalyzer.parseFile()
	if output != "id":
		blastAnalyzer.printSummary()
	else:
		blastAnalyzer.printAllIDs()
#	blastAnalyzer.printReport()
