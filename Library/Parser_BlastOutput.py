#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-27

"""
Create tokens from blast output file (this is a library)
"""

class Entry:
	def __init__(self):
		self.reset()

	def reset(self):
		self.data = {}
		self.data["length"] = 0
		self.data["query_fullName"] = ""
		self.data["query_start"] = 0
		self.data["query_end"] = 0
		self.data["query_sequence"] = ""
		self.data["start"] = 0
		self.data["end"] = 0
		self.data["sequence"] = ""
		self.data["score"] = 0
		self.data["e_value"] = 0
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
		self.species = {}
		self.sequences = {}
		self.query = Entry()

	def setID(self, ID):
		self.ID = ID

	def setCount(self, count):
		self.count = count
	
	def getCount(self):
		return self.count

	def getShortName(self, name):
		return str(name).split()[0]

	def addName(self, name):
		self.subjcts_names.append(name)
#		if name not in self.subjcts_names:
#			self.subjcts_names.append(name)

	def setData(self, name, datatype, value):
		if name == "Query":
			self.query.setData(datatype, value)
		else:
			shortName = self.getShortName(name)
			if shortName not in self.data:
				self.data[shortName] = Entry()
				self.addName(shortName)
			self.data[shortName].setData(datatype, value)

	def removeName(self, name):
		if self.subjcts_names.count(name) == 1:
			self.subjcts_names.remove(name)
			del self.data[name]
		else:
			isRemoved = False
			for i in range(len(self.subjcts_names)-1, -1, -1):
				if isRemoved == False:
					if self.subjcts_names[i] == name:
						self.subjcts_names.pop(i)
						isRemoved = True

	def getData(self, name, datatype):
		if name == "Query":
			return self.query.getData(datatype)
		else:
			return self.data[name].getData(datatype)
	
	def getID(self):
		return self.ID
	
	def getNames(self):
		return self.subjcts_names

	def getSequence(self, name):
		return self.sequences[name]

	def getNumberOfResult(self):
		return len(self.subjcts_names)

	def processScores(self): # Keep only value with the top score
		maximum = 0.0
		# Get max value
		keys = self.getNames()
		for key in keys:
			current_score = self.getData(key, "score")
			if current_score > maximum:
				maximum = current_score

		# Remove keys that are below max score
		toRemove = []
		for key in keys:
			if self.getData(key, "score") < maximum:
				toRemove.append(key)
		for key in toRemove:
			self.removeName(key)

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
#		self.token.setID(tokens[1])
		self.token.setID(' '.join(line.split()[1:len(line.split())]))
		self.token.setCount(int(tokens[3]))

	def fetchSubjctName(self, line):
#		tokens = line.split()
#		name = tokens[0]
#		fullName = ' '.join(tokens[:len(tokens)-2])
		fullName = line[1:].strip()
		name = fullName.split()[0]
		self.token.addName(name)
		self.token.setData(name, "query_fullName", fullName) 

	def fetchScores(self, line):
		tokens = line.split()
		name = line.split()[0]
		score = float(tokens[len(tokens)-2])
		evalue = float(tokens[len(tokens)-1])
		self.setData(name, "score", score)
		self.setData(name, "e_value", evalue)

	def parseScoreLine(self, line):
		self.fetchScores(line)

	def fetchLength(self, name, line):
		tokens = line.split('=')
		self.setData(name, "length", int(tokens[1]))

	def fetchIdentitiesAndGaps(self, line):
		name = self.token.getNames()[self.queryCount]
		tokens = line.split()
		identities = str(tokens[2]).split('/')
		gaps = str(tokens[6]).split('/')
		self.setData(name, "mismatches", int(identities[1]) - int(identities[0]))
		self.setData(name, "gaps", int(gaps[0]))

	def fetchStartEnd(self, name, line):
		if "Query" in line:
			toAddStart = "query_start"
			toAddEnd = "query_end"
		if "Sbjct" in line:
			toAddStart = "start"
			toAddEnd = "end"
		tokens = line.split()
		start = int(tokens[1])
		end = int(tokens[3])
		if end > start:
			self.setData(name, toAddStart, start)
			self.setData(name, toAddEnd, end)
		else:
			self.setData(name, toAddStart, end)
			self.setData(name, toAddEnd, start)

	def fetchSequence(self, name, line):
		tokens = line.split()
		sequence = tokens[2]
		if "Query" in line:
			self.token.setData(name, "query_sequence", sequence)
		if "Sbjct" in line:
			self.token.setData(name, "sequence", sequence)

	def parseSequenceLine(self, line):
		name = self.token.getNames()[self.queryCount]
		self.fetchStartEnd(name, line)
		self.fetchSequence(name, line)

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
#			if '>' in line:
			if line[0] == '>':
				self.setState("scoresFetched")
				self.fetchSubjctName(line)
			else:
				if len(line.strip()) > 0:
					self.parseScoreLine(line)

		elif state == "scoresFetched":
			queryState = self.getQueryState()
			if queryState == "newQuery":
				if "Lambda" in line: # This means end of Entry
					self.setState("noEntry")
				if len(line.strip()) > 0 and line[0] == '>':
					self.fetchSubjctName(line)
				elif "Length=" in line:
#					print "self.queryCount: " + str(self.queryCount)
#					print "len(self.token.getNames()): " + str(len(self.token.getNames()))
					name = self.token.getNames()[self.queryCount]
#					print "name: " + str(name)
					self.fetchLength(name, line)
					self.setQueryState("subjctsLengthFetched")
			elif queryState == "subjctsLengthFetched":
				if "Identities" in line:
#					print "Identities"
					self.fetchIdentitiesAndGaps(line)
					self.setQueryState("identFetched")
			elif queryState == "identFetched":
				if "Query" in line:
					self.parseSequenceLine(line)
					self.setQueryState("queryStartEndFetched")
					
			elif queryState == "queryStartEndFetched":
				if "Sbjct" in line:
					self.parseSequenceLine(line)
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
