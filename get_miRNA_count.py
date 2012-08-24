#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-24

"""
usage:
get_miRNA_count.py blast.out
"""

class Parser:
	def __init__ (self):
		self.inEntry = False
		self.hasHit = False
		self.scoresFetched = False
		self.reset()
	
	def reset(self):
		self.count = {}

	def newEntry(self):
		self.inEntry = True
		self.hasHit = False
		self.scoresFetched = False
		self.entryID = ""
		self.entryCount = 0
		self.entryScores = {}

	def fetchLineInfos(self, line):
		tokens = line.split()
		self.entryID = tokens[1] # Is it useful?
		self.entryCount = int(tokens[3])

	def fetchScore(self, line):
		tokens = line.split()
		name = tokens[0]
		i = 1
		while tokens[i] != "stem-loop":
			i += 1
		score = float(tokens[i+1])
		self.entryScores[name] = score

	def addCount(self, key, count):
		if key in self.count:
			self.count[key] += count
		else:
			self.count[key] = count

	def processScores(self):
		max = 0.0
		# Get max value
		for key in self.entryScores:
			if self.entryScores[key] > max:
				max = self.entryScores[key]

		# Get number of miRNA with max value
		tie = 0
		for key in self.entryScores:
			if self.entryScores[key] >= max:
				tie += 1	

		# Calculate the count considering the number of tie
		count = self.entryCount / float(tie)

		# Add the count to the miRNA
		for key in self.entryScores:
			if self.entryScores[key] >= max:
				self.addCount(key, count)

	def parseLine(self, line):
		if self.inEntry == False:
			if "Query=" in line:
				self.newEntry()
				self.fetchLineInfos(line)   
		elif self.hasHit == False:
			if "No hits found" in line:
				self.inEntry = False

			elif "Sequences producing" in line:
				self.hasHit = True
		elif self.hasHit == True:
			if '>' in line:
				self.scoresFetched = True
				self.processScores()
				self.inEntry = False
				self.hasHit = False
			elif self.scoresFetched == False:
				if "hsa" in line:
					self.fetchScore(line)

	def printAll(self):
		for key in self.count:
			print key + '\t' + str(self.count[key])
				
   

import sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print __doc__
		sys.exit(1)

	parser = Parser()
	filename=sys.argv[1]
	for line in open(filename):		
		parser.parseLine(line)
	parser.printAll()
