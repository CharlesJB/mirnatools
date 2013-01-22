#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-27

"""
Based on the algorithm presented in: Creighton et al 2009
Usage:
preethi_blast_analysis.py <blast.out> <output> <seq_count>
    blast.out: Output of blast against miRBase (tested with version 19).
    output: Prefix that will be added to output files.
    seq_count: The number of sequences after quality trimming.
"""

from Library/Parser_BlastOutput import *
from Library/SpecieConverter.py import *

class BlastAnalyzer:
	def __init__(self, filename, output, total_count):
		self.parser = Parser(filename)
		self.output = output
		self.total_count = float(total_count)
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

	def processScores(self, token):
		maximum = 0.0
		# Get max value
		keys = token.getNames()
		for key in keys:
			current_score = token.getData(key, "score")
			if current_score > maximum:
				maximum = current_score

		# Remove keys that are below max score
		toRemove = []
		for key in keys:
			if token.getData(key, "score") < maximum:
				toRemove.append(key)
		for key in toRemove:
			token.removeName(key)


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
		self.processScores(token)

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
				print count, " blast hits processed."
			if self.parser.isEOF() != True:
				token = self.parser.getToken()
				self.processToken(token)
			else:
				print count, " blast hits processed."
				done = True

	def printCount(self, filename, container):
		f = open(filename, 'w')
		for miRNA in container:
			f.write(miRNA + '\t' + str(container[miRNA]) + '\n')
		f.close()

	def printID(self, filename, container):
		f = open(filename, 'w')
		for miRNA in container:
			toPrint = miRNA
			for i in range(0, len(container[miRNA])):
				toPrint = toPrint + '\t' + container[miRNA][i]
			toPrint += '\n'
			f.write(toPrint)
		f.close()

	def printReport(self):
		specieConverter = SpecieConverter()
		filename = self.output + "_perfectMatches_summary.txt"
		f = open(filename, 'w')
		f.write("miRNA_ID\tSpecie\tmiRNA_Sequence\tSequence_Count\t%_of_total\n")
		for miRNA in self.perfectCounts:
			toPrint = miRNA
			toPrint = toPrint + "\t" + self.getSpecieName(miRNA)
			toPrint = toPrint + "\t" + specieConverter(self.getSpecieName(miRNA))
			toPrint = toPrint + "\t" + self.getQuerySequence(miRNA)
			toPrint = toPrint + "\t" + str(self.perfectCounts[miRNA])
			toPrint = toPrint + "\t" + str((self.perfectCounts[miRNA] / self.total_count) * 100)
			toPrint = toPrint + '\n'
			f.write(toPrint)
		f.close()

	def printAll(self):
		# Print perfect matches
		filename = output + "_perfectMatches.txt"
		self.printCount(filename, self.perfectCounts)
		filename = output + "_perfectSeq_ID.txt"
		self.printID(filename, self.perfectSeqIDs)
		# Print loose matches
		filename = output + "_looseMatches.txt"
		self.printCount(filename, self.looseCounts)
		filename = output + "_looseSeq_ID.txt"
		self.printID(filename, self.looseSeqIDs)
		# Print combined (loose + perfect) matches
		filename = output + "_combinedMatches.txt"
		self.printCount(filename, self.combinedCounts)
		# Print no matches
		filename = output +  "_noMatchLength.txt"
		self.printCount(filename, self.noMatchLength)
		filename = output +  "_noMatchGaps.txt"
		self.printCount(filename, self.noMatchGaps)
		filename = output +  "_noMatchMismatches.txt"
		self.printCount(filename, self.noMatchMismatches)
		filename = output +  "_noMatchLength_ID.txt"
		self.printID(filename, self.noMatchLength_IDs)
		filename = output +  "_noMatchGaps_ID.txt"
		self.printID(filename, self.noMatchGaps_IDs)
		filename = output +  "_noMatchMismatches_ID.txt"
		self.printID(filename, self.noMatchMismatches_IDs)

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
	blastAnalyzer.printAll()
	blastAnalyzer.printReport()
