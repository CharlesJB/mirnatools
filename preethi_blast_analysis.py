#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-27

"""
Based on the algorithm presented in: Creighton et all 2009
usage:
preethi_blast_analysis.py <blast.out> <output>
blast.out: Output of blast
output: Prefix that will be added to output files
"""

class Entry:
	def __init__(self):
		self.reset()

	def reset(self):
		self.data = {}
		self.data["length"] = 0
		self.data["query_start"] = 0
		self.data["query_end"] = 0
		self.data["start"] = 0
		self.data["end"] = 0
		self.data["score"] = 0
		self.data["mismatches"] = 0
		self.data["gaps"] = 0

	def setData(self, name, value):
		self.data[name] = value

	def getData(self, name):
		return self.data[name]

class Token:
	def __init__(self):
		self.reset()

	def reset(self):
		self.ID = ""
		self.count = 0
		self.data = {}
		self.subjcts_names = []
		self.query = Entry()

	def setID(self, ID):
		self.ID = ID

	def setCount(self, count):
		self.count = count
	
	def getCount(self):
		return self.count

	def addName(self, name):
		if name not in self.subjcts_names:
			self.subjcts_names.append(name)

	def setData(self, name, datatype, value):
		if name == "Query":
			self.query.setData(datatype, value)
		else:
			if name not in self.data:
				self.data[name] = Entry()
				self.addName(name)
			self.data[name].setData(datatype, value)

	def removeName(self, name):
		self.subjcts_names.remove(name)
		del self.data[name]

	def getData(self, name, datatype):
		if name == "Query":
			return self.query.getData(datatype)
		else:
			return self.data[name].getData(datatype)
	
	def getID(self):
		return self.ID
	
	def getNames(self):
		return self.subjcts_names
	
	def getNumberOfResult(self):
		return len(self.subjcts_names)

class Parser:
	def __init__ (self, filename):
		self.f = open(filename) 
		self.token = Token()
		self.eof = False
		self.state = "noEntry"
		self.queryState = "newQuery"
		self.queryCount = 0

	def setQueryState(self, state):
		self.queryState = state
	
	def getQueryState(self):
		return self.queryState
	
	def setState(self, state):
		self.state = state

	def getState(self):
		return self.state

	def setData(self, name, datatype, value):
		self.token.setData(name, datatype, value)

	def newEntry(self):
		self.token.reset()
		self.setState("noEntry")
		self.setQueryState("newQuery")
		self.queryCount = 0

	def fetchQueryInfos(self, line):
		tokens = line.split()
		self.token.setID(tokens[1])
		self.token.setCount(int(tokens[3]))

	def fetchLength(self, name, line):
		tokens = line.split('=')
		self.setData(name, "length", int(tokens[1]))

	def fetchScore(self, line):
		tokens = line.split()
		name = tokens[0]
#		score = float(tokens[5])
		score = float(tokens[len(tokens)-2])
		self.token.setData(name, "score", score)

	def fetchIdentitiesAndGaps(self, name, line):
		tokens = line.split()
		identities = str(tokens[2]).split('/')
		self.setData(name, "mismatches", int(identities[1]) - int(identities[0]))
		gaps = str(tokens[6]).split('/')
		self.setData(name, "gaps", int(gaps[0]))

	def fetchStartEnd(self, name, line):
		toAddStart = "start"
		toAddEnd = "end"
		if "Query" in line:
			toAddStart = "query_start"
			toAddEnd = "query_end"
		tokens = line.split()
		start = int(tokens[1])
		end = int(tokens[3])
		if end > start:
			self.setData(name, toAddStart, start)
			self.setData(name, toAddEnd, end)
		else:
			self.setData(name, toAddStart, end)
			self.setData(name, toAddEnd, start)

	def addCount(self, count):
		self.token.addCount(count)

	def parseLine(self, line):
		state = self.getState()
		if state == "newEntry" or state == "noEntry":
			if state == "noEntry":
				self.newEntry()
				self.setState("newEntry")
			if "Query=" in line:
				self.fetchQueryInfos(line)   
				self.setState("inEntry")
		elif state == "inEntry":
			if "Length=" in line:
				self.fetchLength("Query", line)
				self.setState("queryLengthFetched")

		elif state == "queryLengthFetched":
			if "No hits found" in line:
				self.newEntry()

			elif "Sequences producing" in line:
				self.setState("hasHit")

		elif state == "hasHit":
			if '>' in line:
				self.setState("scoresFetched")
			else:
				if "hsa" in line:
					self.fetchScore(line)

		elif state == "scoresFetched":
			queryState = self.getQueryState()
			if queryState == "newQuery":
				if "Lambda" in line: # This means end of Entry
					self.setState("noEntry")
				elif "Length=" in line:
					name = self.token.getNames()[self.queryCount]
					self.fetchLength(name, line)
					self.setQueryState("subjctsLengthFetched")
			elif queryState == "subjctsLengthFetched":
				if "Identities" in line:
					name = self.token.getNames()[self.queryCount]
					self.fetchIdentitiesAndGaps(name, line)
					self.setQueryState("identFetched")
			elif queryState == "identFetched":
				if "Query" in line:
					name = self.token.getNames()[self.queryCount]
					self.fetchStartEnd(name, line)
					self.setQueryState("queryStartEndFetched")
					
			elif queryState == "queryStartEndFetched":
				if "Sbjct" in line:
					name = self.token.getNames()[self.queryCount]
					self.fetchStartEnd(name, line)
					self.queryCount += 1
					self.setQueryState("newQuery")
					
	def getToken(self):
		return self.token

	def isEOF(self):
		return self.eof

	def createNextToken(self):
		done = False
		while done == False and self.isEOF() == False:
			line = self.f.readline()
			self.parseLine(line)
			if not line:
				self.eof = True
			elif self.getState() == "noEntry":
				done = True
   
class BlastAnalyzer:
	def __init__(self, filename, output):
		self.parser = Parser(filename)
		self.output = output
		self.perfectCounts = {}
		self.perfectSeqIDs = {} # The name of the sequences that blast to each miRNA
		self.looseCounts = {}
		self.looseSeqIDs = {} # The name of the sequences that blast to each miRNA
		self.noMatchLength = {}
		self.noMatchMismatches = {}
		self.noMatchGaps = {}

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

	def addResults(self, result, token, i):
		count = float(token.getCount()) / float(token.getNumberOfResult())
		name = token.getNames()[i]
		self.addCount(name, result, count)
		self.addSeqName(name, result, token.getID())

	def addNoMatch(self, result, token, i):
		count = float(token.getCount()) / float(token.getNumberOfResult())
		name = token.getNames()[i]
		if result == "noMatch, length":
			if name in self.noMatchLength:
				self.noMatchLength[name] += count
			else:
				self.noMatchLength[name] = count
		elif result == "noMatch, gaps":
			if name in self.noMatchGaps:
				self.noMatchGaps[name] += count
			else:
				self.noMatchGaps[name] = count
		elif result == "noMatch, mismatches":
			if name in self.noMatchMismatches:
				self.noMatchMismatches[name] += count
			else:
				self.noMatchMismatches[name] = count

	def processToken(self, token):
		# We keep only miRNA having the best score
		self.processScores(token)

		# Check if we have a perfect match, a loose match or no match at all
		for i in range (0, token.getNumberOfResult()):
			result = self.processMatch(token, i)
			if result == "looseMatch" or result == "perfectMatch":
				self.addResults(result, token, i)
			else:
				self.addNoMatch(result, token, i)
				
	def parseFile(self):
		done = False
		count = 0
		while done == False:
#			print "====== Token: " + str(count+1) + " ======"
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

	def printAll(self):
		# Print perfect matches
		filename = output + "_perfectMatches.txt"
		f = open(filename, 'w')
		for miRNA in self.perfectCounts:
			f.write(miRNA + '\t' + str(self.perfectCounts[miRNA]) + '\n')
		f.close()
		filename = output + "_perfectSeqID.txt"
		f = open(filename, 'w')
		for miRNA in self.perfectSeqIDs:
			toPrint = miRNA
			for i in range(0, len(self.perfectSeqIDs[miRNA])):
				toPrint = toPrint + '\t' + self.perfectSeqIDs[miRNA][i]
			toPrint += '\n'
			f.write(toPrint)
		f.close()

		# Print loose matches
		filename = output + "_looseMatches.txt"
		f = open(filename, 'w')
		for miRNA in self.looseCounts:
			f.write(miRNA + '\t' + str(self.looseCounts[miRNA]) + '\n')
		f.close()
		filename = output + "_looseSeqID.txt"
		f = open(filename, 'w')
		for miRNA in self.looseSeqIDs:
			toPrint = miRNA
			for i in range(0, len(self.looseSeqIDs[miRNA])):
				toPrint = toPrint + '\t' + self.looseSeqIDs[miRNA][i]
			toPrint += '\n'
			f.write(toPrint)
		f.close()

		# Print sequences that didn't made it
		filename = output +  "_noMatchLength.txt"
		f = open(filename, 'w')
		for miRNA in self.noMatchLength:
			f.write(miRNA + '\t' + str(self.noMatchLength[miRNA]) + '\n')
		f.close()
		filename = output +  "_noMatchGaps.txt"
		f = open(filename, 'w')
		for miRNA in self.noMatchGaps:
			f.write(miRNA + '\t' + str(self.noMatchGaps[miRNA]) + '\n')
		f.close()
		filename = output +  "_noMatchMismatches.txt"
		f = open(filename, 'w')
		for miRNA in self.noMatchMismatches:
			f.write(miRNA + '\t' + str(self.noMatchMismatches[miRNA]) + '\n')
		f.close()

import sys

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(1)

	filename = sys.argv[1]
	output = sys.argv[2]
	blastAnalyzer = BlastAnalyzer(filename, output)
	blastAnalyzer.parseFile()
	blastAnalyzer.printAll()
