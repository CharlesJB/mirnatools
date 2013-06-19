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
		self.reverse_complement = False

	def setReverseComplement(self, value):
		self.reverse_complement = value

	def getReverseComplement(self):
		return self.reverse_complement

	def setID(self, ID):
		self.ID = ID

	def getID(self):
		return self.ID

	def setCount(self, count):
		if count == "":
			count = 0
		self.count = count
	
	def getCount(self):
		return self.count

	def getShortName(self, name):
		return str(name).split()[0]

	def addName(self, name):
		self.subjcts_names.append(name)

	def getNames(self):
		return self.subjcts_names


	def setData(self, name, datatype, value):
		if name == "Query":
			self.query.setData(datatype, value)
		else:
			shortName = self.getShortName(name)
			if shortName not in self.data:
				self.data[shortName] = Entry()
				self.addName(shortName)
			self.data[shortName].setData(datatype, value)

	def getData(self, name, datatype):
		if name == "Query":
			return self.query.getData(datatype)
		else:
			return self.data[name].getData(datatype)
	
	def removeName(self, name):
		if self.subjcts_names.count(name) == 1:
			self.subjcts_names.remove(name)
			if name in self.data:
				del self.data[name]
		else:
			isRemoved = False
			for i in range(len(self.subjcts_names)-1, -1, -1):
				if isRemoved == False:
					if self.subjcts_names[i] == name:
						self.subjcts_names.pop(i)
						isRemoved = True

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
		self.token.setID(tokens[1])
#		self.token.setID(' '.join(line.split()[1:len(line.split())]))
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
	
	def setReverseComplementTag(self):
		self.token.setReverseComplement(True)

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
			self.setReverseComplementTag()
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
#		self.token.addCount(count)
		self.token.setCount(self.token.getCount() + count)

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
					self.fetchScores(line)

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

import sys

# Tests
if __name__ == "__main__":
	tests_count = 0
	fail_count = 0
	error_messages = "Failed tests: \n"

	#***************************
	# Tests for class: Entry
	#***************************

	#### Set/Get
	# Preset Value
	setGetEntry = Entry()
	setGetEntry.setData("length", 100)
	if setGetEntry.getData("length") != 100:
		fail_count += 1
		error_messages += "\n\tEntry: Set/Get, preset values."
	tests_count += 1
	# Custom Value
	setGetEntry = Entry()
	setGetEntry.setData("asdf", 200)
	if setGetEntry.getData("asdf") != 200:
		fail_count += 1
		error_messages += "\n\tEntry: Set/Get, custom values."
	tests_count += 1
	# Invalid Value
	setGetEntry = Entry()
	try:
		a = setGetEntry.getData("qwert")
		fail_count += 1
		error_messages += "\n\tEntry: Set/Get, invalid values."
	except KeyError:
		a = 1
	tests_count += 1

	#***************************
	# Tests for class: Token"
	#***************************

	#### setReverseComplement/getReverseComplement 
	# Normal Usage
	setGetReverse= Token()
	setGetReverse.setReverseComplement("acgt")
	if setGetReverse.getReverseComplement() != "acgt":
		fail_count += 1
		error_messages += "\n\tToken: setReverseComplement/getReverseComplement, Normal Usage."
	tests_count += 1
	# Empty Values
	setGetReverse= Token()
	setGetReverse.setReverseComplement("")
	if setGetReverse.getReverseComplement() != "":
		fail_count += 1
		error_messages += "\n\tToken: setReverseComplement/getReverseComplement, Empty Values."
	tests_count += 1

	#### setID/getID
	# Normal Usage
	setGetID = Token()
	setGetID.setID("PP_01")
	if setGetID.getID() != "PP_01":
		fail_count += 1
		error_messages += "\n\tToken: setID/getID, Normal Usage."
	tests_count += 1
	# Empty Values
	setGetID = Token()
	setGetID.setID("")
	if setGetID.getID() != "":
		fail_count += 1
		error_messages += "\n\tToken: setID/getID, Empty Values."
	tests_count += 1

	#### setCount/getCount
	# Normal Usage
	setGetCount = Token()
	setGetCount.setCount(10)
	if setGetCount.getCount() != 10: 
		fail_count += 1
		error_messages += "\n\tToken: setCount/getCount, Normal Usage."
	tests_count += 1
	# Empty Values
	setGetCount = Token()
	setGetCount.setCount("")
	if setGetCount.getCount() != 0:
		fail_count += 1
		error_messages += "\n\tToken: setCount/getCount, Empty Values."
	tests_count += 1

	#### getShortName
	# Normal Usage
	getShortName = Token()
	if getShortName.getShortName("hsa-miR-200b-3p MIMAT0000318 Homo sapiens miR-200b-3p") != "hsa-miR-200b-3p":
		fail_count += 1
		error_messages += "\n\tToken: getShortName, Normal Usage."
	tests_count += 1
	# Empty Values
	getShortName = Token()
	try:
		getShortName.getShortName("")
		fail_count += 1
		error_messages += "\n\tToken: getShortName, Empty Values."
	except IndexError:
		a = 1
	tests_count += 1
		
	#### addName/getNames
	getName = Token()
	# No Value
	if len(getName.getNames()) != 0:
		fail_count += 1
		error_messages += "\n\tToken: addName/getNames, No value."
	tests_count += 1
	# One Value
	getName.addName("abc")
	if len(getName.getNames()) != 1 or getName.getNames()[0] != "abc":
		fail_count += 1
		error_messages += "\n\tToken: addName/getNames, One Value."
	tests_count += 1
	# Multiple Values
	getName.addName("def")
	if len(getName.getNames()) != 2 or getName.getNames()[0] != "abc" or getName.getNames()[1] != "def":
		fail_count += 1
		error_messages += "\n\tToken: addName/getNames, Multiple values."
	tests_count += 1
	#### setData/getData
	# Query
	setData = Token()
	setData.setData("Query", "query_sequence", "ggtt")
	if setData.getData("Query", "query_sequence") != "ggtt" or setData.query.getData("query_sequence") != "ggtt":
		fail_count += 1
		error_messages += "\n\tToken: setData/getData, Query."
	tests_count += 1
	# Other
	setData = Token()
	setData.setData("mir-123", "Length", 10)
	if setData.getData("mir-123", "Length") != 10 or setData.data["mir-123"].getData("Length") != 10:
		fail_count += 1
		error_messages += "\n\tToken: setData/getData, Other."
	tests_count += 1
	# Other Long Name
	setData = Token()
	setData.setData("hsa-miR-200b-3p MIMAT0000318 Homo sapiens miR-200b-3p", "Sequence", "acct")
	if setData.getData("hsa-miR-200b-3p", "Sequence") != "acct" or setData.data["hsa-miR-200b-3p"].getData("Sequence") != "acct":
		fail_count += 1
		error_messages += "\n\tToken: setData/getData, Other Long Name."
	tests_count += 1
	# Invalid Name
	setData = Token()
	try:
		setData.getData("hsa-miR-200b-3p", "Sequence")
		fail_count += 1
		error_messages += "\n\tToken: setData/getData, Invalid Name."
	except KeyError:
		a = 1
	tests_count += 1
	# Invalid Value
	setData = Token()
	setData.setData("mir-123", "Length", 10)
	try:
		setData.getData("mir-123", "Sequence")	
		fail_count += 1
		error_messages += "\n\tToken: setData/getData, Invalid Value."
	except KeyError:
		a = 1
	tests_count += 1
	#### removeName
	# One occurence First Name
	removeName = Token()
	removeName.addName("mir-456")
	removeName.removeName("mir-456")
	if removeName.subjcts_names.count("mir-456") != 0:
		fail_count += 1
		error_messages += "\n\tToken: removeName, One Occurence First Name."
	tests_count += 1
	# One occurence First Name With Data
	removeName = Token()
	removeName.setData("mir-789", "Sequence", "aaaT")
	removeName.removeName("mir-789")
	if removeName.subjcts_names.count("mir-789") != 0:
		fail_count += 1
		error_messages += "\n\tToken: removeName, One Occurence First Name With Data Test #1."
	else:
		try:
			removeName.getData("mir-789", "Sequence")
			fail_count += 1
			error_messages += "\n\tToken: removeName, One Occurence First Name With Data Test #2."
		except KeyError:
			a = 1
	tests_count += 1
	# One occurence Not First Name
	removeName = Token()
	removeName.addName("mir-123")
	removeName.addName("mir-456")
	removeName.removeName("mir-456")
	if removeName.subjcts_names.count("mir-456") != 0 or removeName.subjcts_names.count("mir-123") != 1:
		fail_count += 1
		error_messages += "\n\tToken: removeName, One occurence Not First Name."
	tests_count += 1
	# Multiple occurences First Name
	removeName = Token()
	removeName.addName("mir-456")
	removeName.addName("mir-456")
	removeName.removeName("mir-456")
	if removeName.subjcts_names.count("mir-456") != 1:
		fail_count += 1
		error_messages += "\n\tToken: removeName, Multiple Occurences First Name."
	tests_count += 1
	# Multiple occurences First Name With Data
	removeName = Token()
	removeName.setData("mir-789", "Sequence", "aaaT")
	removeName.addName("mir-789")
	removeName.addName("mir-345")
	removeName.removeName("mir-789")
	if removeName.subjcts_names.count("mir-789") != 1 or removeName.getData("mir-789", "Sequence") != "aaaT":
		fail_count += 1
		error_messages += "\n\tToken: removeName, Multiple Occurences First Name With Data Test #1."
	else:
		removeName.removeName("mir-789")
		try:
			removeName.getData("mir-789", "Sequence")
			fail_count += 1
			error_messages += "\n\tToken: removeName, Multiple Occurences First Name With Data Test #2."
		except KeyError:
			a = 1
	tests_count += 1
	# Invalid Name
	removeName = Token()
	try:
		removeName.removeName("mir-123")
	except KeyError:
		error_messages += "\n\tToken: removeName, Invalid Name."
		fail_count += 1
	tests_count += 1
	#### getNumberOfResult
	# No Result
	getNumberOfResult = Token()
	if getNumberOfResult.getNumberOfResult() != 0:
		fail_count += 1
		error_messages += "\n\tToken: getNumberOfResult, No Result."
	tests_count += 1
	# One Result
	getNumberOfResult = Token()
	getNumberOfResult.addName("mir-556")
	if getNumberOfResult.getNumberOfResult() != 1:
		fail_count += 1
		error_messages += "\n\tToken: getNumberOfResult, One Result."
	tests_count += 1
	# Two Results
	getNumberOfResult = Token()
	getNumberOfResult.addName("mir-566")
	getNumberOfResult.addName("mir-565")
	if getNumberOfResult.getNumberOfResult() != 2:
		fail_count += 1
		error_messages += "\n\tToken: getNumberOfResult, Two Results."
	tests_count += 1
	#### processScores
	# No Score
	processScores = Token()
	processScores.processScores()
	if processScores.getNumberOfResult() != 0:
		fail_count += 1
		error_messages += "\n\tToken: processScores, No Score."
	tests_count += 1
	# One Score
	processScores = Token()
	processScores.setData("mir-112", "score", 20)
	processScores.processScores()
	if processScores.getNumberOfResult() != 1 or processScores.getData("mir-112", "score") != 20:
		fail_count += 1
		error_messages += "\n\tToken: processScores, One Score."
	tests_count += 1
	# Multiple Scores Different Values
	processScores = Token()
	processScores.setData("mir-113", "score", 20)
	processScores.setData("mir-114", "score", 30)
	processScores.processScores()
	if processScores.getNumberOfResult() != 1 or processScores.getData("mir-114", "score") != 30:
		fail_count += 1
		error_messages += "\n\tToken: processScores, Multiples Scores Different Values Test #1."
	else:
		try:
			processScores.getData("mir-113", "score")
			fail_count += 1
			error_messages += "\n\tToken: processScores, Multiples Scores Different Values Test #2."
		except KeyError:
			a=1
	tests_count += 1
	# Multiple Scores Same Values
	processScores = Token()
	processScores.setData("mir-115", "score", 30)
	processScores.setData("mir-116", "score", 30)
	processScores.processScores()
	if processScores.getNumberOfResult() != 2 or processScores.getData("mir-115", "score") != 30 or processScores.getData("mir-116", "score") != 30:
		fail_count += 1
		error_messages += "\n\tToken: processScores, Multiples Scores Same Values."
	tests_count += 1


	#***************************
	# Tests for class: Parser"
	#***************************

	#### setQueryState/getQueryState
	# Normal Usage
	setQueryState = Parser("test.txt")
	setQueryState.setQueryState("tested")
	if setQueryState.getQueryState() != "tested" or setQueryState.queryState != "tested":
		fail_count += 1
		error_messages += "\n\tParser: setQueryState, Normal Usage."
	tests_count += 1
	# No Value
	setQueryState = Parser("test.txt")
	setQueryState.setQueryState("")
	if setQueryState.getQueryState() != "" or setQueryState.queryState != "":
		fail_count += 1
		error_messages += "\n\tParser: setQueryState, No Value."
	tests_count += 1
	#### setState/getState
	# Normal Usage
	setState = Parser("test.txt")
	setState.setState("tested")
	if setState.getState() != "tested" or setState.state != "tested":
		fail_count += 1
		error_messages += "\n\tParser: setState, Normal Usage."
	tests_count += 1
	# No Value
	setState = Parser("test.txt")
	setState.setState("")
	if setState.getState() != "" or setState.state != "":
		fail_count += 1
		error_messages += "\n\tParser: setState, No Value."
	tests_count += 1
	#### fetchQueryInfos
	# Valid Query Line
	query_line="Query= PP_7 Count: 58911"
	fetchQueryInfos = Parser("test.txt")
	fetchQueryInfos.fetchQueryInfos(query_line)
	if fetchQueryInfos.token.getID() != "PP_7" or fetchQueryInfos.token.getCount() != 58911:
		fail_count += 1
		error_messages += "\n\tParser: fetchQueryInfos, Valid Query Line."
	tests_count += 1
	#### fetchSubjctName
	# Valid Subjct Line
	subjct_line="> hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168"
	fetchSubjctName = Parser("test.txt")
	fetchSubjctName.fetchSubjctName(subjct_line)
	name = fetchSubjctName.token.getNames()[0]
	if name != "hsa-miR-3168" or fetchSubjctName.token.getData(name, "query_fullName") != subjct_line[2:]:
		fail_count += 1
		error_messages += "\n\tParser: fetchSubjctName, Valid Subjct Line."
	tests_count += 1
	#### fetchScores
	# Valid Score Line
	score_line="  hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168                     26.3    0.003"
	fetchScores= Parser("test.txt")
	fetchScores.fetchScores(score_line)
	name = fetchScores.token.getNames()[0]
	if name != "hsa-miR-3168" or fetchScores.token.getData(name, "score") != 26.3 or fetchScores.token.getData(name, "e_value") != 0.003:
		fail_count += 1
		error_messages += "\n\tParser: fetchScores, Valid Score Line."
	tests_count += 1
	#### fetchLength
	# Valid Length Line
	length_line="Length=19"
	name="mir-111"
	fetchLength = Parser("test.txt")
	fetchLength.fetchLength(name, length_line)
	if fetchLength.token.getData(name, "length") != 19:
		fail_count += 1
		error_messages += "\n\tParser: fetchLength, Valid Length Line."
	tests_count += 1
	#### fetchIdentitiesAndGaps
	# Valid Identity/Gap Line
	identity_line=" Identities = 12/13 (100%), Gaps = 1/13 (0%)"
	score_line="  hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168                     26.3    0.003"
	fetchIdentitiesAndGaps = Parser("test.txt")
	fetchIdentitiesAndGaps.fetchScores(score_line)
	fetchIdentitiesAndGaps.fetchIdentitiesAndGaps(identity_line)
	name = fetchIdentitiesAndGaps.token.getNames()[0]
	if fetchIdentitiesAndGaps.token.getData(name, "mismatches") != 1 or fetchIdentitiesAndGaps.token.getData(name, "gaps") != 1:
		fail_count += 1
		error_messages += "\n\tParser: fetchIdentitiesAndGaps, Valid Identities/Gaps Line."
	tests_count += 1
	#### fetchStartEnd
	# Valid StartEnd Query Line
	fetchStartEnd = Parser("test.txt")
	startend_query_line = "Query  2   GAGTTCTACAGTC  14"
	name = "mir-333"
	fetchStartEnd.fetchStartEnd(name, startend_query_line)
	if fetchStartEnd.token.getData(name, "query_start") != 2 or fetchStartEnd.token.getData(name, "query_end") != 14:
		fail_count += 1
		error_messages += "\n\tParser: fetchStartEnd, Valid StartEnd Query Line."
	tests_count += 1
	# Valid StartEnd Subjct Line
	startend_subjct_line = "Sbjct  1   GAGTTCTACAGTC  13"
	name = "mir-333"
	fetchStartEnd.fetchStartEnd(name, startend_subjct_line)
	if fetchStartEnd.token.getData(name, "start") != 1 or fetchStartEnd.token.getData(name, "end") != 13:
		fail_count += 1
		error_messages += "\n\tParser: fetchStartEnd, Valid StartEnd Subjct Line."
	tests_count += 1
	#### fetchSequence
	# Valid Sequence Query Line
	fetchSequence = Parser("test.txt")
	sequence_query_line = "Query  2   GAGTTCTACAGTC  14"
	name = "mir-333"
	fetchSequence.fetchSequence(name, sequence_query_line)
	if fetchSequence.token.getData(name, "query_sequence") != "GAGTTCTACAGTC":
		fail_count += 1
		error_messages += "\n\tParser: fetchSequence, Valid Sequence Query Line."
	tests_count += 1
	# Valid Sequence Subjct Line
	fetchSequence = Parser("test.txt")
	sequence_subjct_line = "Sbjct  1   GAGTTCTACAGTC  13"
	name = "mir-333"
	fetchSequence.fetchSequence(name, sequence_subjct_line)
	if fetchSequence.token.getData(name, "sequence") != "GAGTTCTACAGTC":
		fail_count += 1
		error_messages += "\n\tParser: fetchSequence, Valid Sequence Subjct Line."
	tests_count += 1
	#### parseSequenceLine
	# Valid Sequence Line
	sequence_line = "Query  2   GAGTTCTACAGTC  14"
	score_line="  hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168                     26.3    0.003"
	parseSequenceLine = Parser("test.txt")
	parseSequenceLine.fetchScores(score_line)
	name = parseSequenceLine.token.getNames()[0]
	parseSequenceLine.parseSequenceLine(sequence_line)
	if parseSequenceLine.token.getData(name, "query_start") != 2 or parseSequenceLine.token.getData(name, "query_end") != 14 or parseSequenceLine.token.getData(name, "query_sequence") != "GAGTTCTACAGTC":
		fail_count += 1
		error_messages += "\n\tParser: parseSequenceLine, Valid Sequence Line"
	tests_count += 1
	#### addCount
	# No Count
	addCount = Parser("test.txt")
	addCount.addCount(10)
	if addCount.token.getCount() != 10:
		fail_count += 1
		error_messages += "\n\tParser: addCount, No Count."
	tests_count += 1
	# One Count
	addCount.addCount(5)
	if addCount.token.getCount() != 15:
		fail_count += 1
		error_messages += "\n\tParser: addCount, One Count."
	tests_count += 1
	#### parseLine
	parseLine = Parser("test.txt")
	if parseLine.queryState != "newQuery" or parseLine.state != "noEntry":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, states are incorrectly initialized with constructor."
	tests_count += 1
	# Query Line 
	line_query = "Query= PP_7 Count: 58911\n"
	parseLine.parseLine(line_query)
	if parseLine.token.getID() != "PP_7" or parseLine.token.getCount() != 58911 or parseLine.state != "inEntry":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Query Line."
	tests_count += 1
	# Query Length Line
	line_query_length = "Length=19\n"
	parseLine.parseLine(line_query_length)
	name = "Query"
	if parseLine.token.getData(name, "length") != 19 or parseLine.state != "queryLengthFetched":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Query Length Line."
	tests_count += 1
	# Procucing Line
	line_producing= "Sequences producing significant alignments:                          (Bits)  Value\n"
	parseLine.parseLine(line_producing)
	if parseLine.state != "hasHit":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Producing Line."
	tests_count += 1
	# Score Line
	line_score = "  hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168                     26.3    0.003\n"
	parseLine.parseLine(line_score)
	name = parseLine.token.getNames()[0]
	if name != "hsa-miR-3168" or parseLine.token.getData(name, "score") != 26.3 or parseLine.token.getData(name, "e_value") != 0.003:
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Score Line."
	tests_count += 1
	# Subjct Name Line
	line_subjct_name = "> hsa-miR-3168 MIMAT0015043 Homo sapiens miR-3168\n"
	parseLine.parseLine(line_subjct_name)
	name = parseLine.token.getNames()[0]
	if name != "hsa-miR-3168" or parseLine.token.getData(name, "query_fullName") != line_subjct_name[2:].strip() or parseLine.state != "scoresFetched":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Valid Subjct Line."
	tests_count += 1
	# Subjct Length Line 
	line_length = "Length=17\n"
	parseLine.parseLine(line_length)
	name = parseLine.token.getNames()[0]
	if parseLine.token.getData(name, "length") != 17 or parseLine.state != "scoresFetched" or parseLine.queryState != "subjctsLengthFetched":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Subjct Length Line."
	tests_count += 1
	# Identity Line
	line_identity = " Identities = 12/13 (100%), Gaps = 1/13 (0%)\n"
	parseLine.parseLine(line_identity)
	name = parseLine.token.getNames()[0]
	if parseLine.token.getData(name, "mismatches") != 1 or parseLine.token.getData(name, "gaps") != 1 or parseLine.state != "scoresFetched" or parseLine.queryState != "identFetched":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Identity Line."
	tests_count += 1
	# Query StartEnd Line
	line_strand= " Strand=Plus/Plus\n"
	line_query_startend = "Query  2   GAGTTCTACAGTC  14\n"
	parseLine.parseLine(line_strand)
	parseLine.parseLine(line_query_startend)
	if parseLine.token.getData(name, "query_start") != 2 or parseLine.token.getData(name, "query_end") != 14 or parseLine.state != "scoresFetched" or parseLine.queryState != "queryStartEndFetched":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Query StartEnd Line."
	tests_count += 1
	# Subjct StartEnd Line
	line_junk= "           |||||||||||||\n"
	line_subjct_startend = "Sbjct  1   GAGTTCTACAGTC  13\n"
	parseLine.parseLine(line_junk)
	parseLine.parseLine(line_subjct_startend)
	if parseLine.token.getData(name, "start") != 1 or parseLine.token.getData(name, "end") != 13 or parseLine.state != "scoresFetched" or parseLine.queryState != "newQuery":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Subjct StartEnd Line."
	tests_count += 1
	# Lambda Line
	line_lambda = "Lambda      K        H\n"
	parseLine.parseLine(line_lambda)
	if parseLine.state != "noEntry":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, Lambda Line."
	tests_count += 1
	# No Hits Found Line
	parseLine = Parser("test.txt")
	line_query = "Query= PP_7 Count: 58911\n"
	line_query_length = "Length=19\n"
	line_noHit= "***** No hits found *****\n"
	parseLine.parseLine(line_query)
	parseLine.parseLine(line_query_length)
	parseLine.parseLine(line_noHit)
	if parseLine.queryState != "newQuery" or parseLine.state != "noEntry":
		fail_count += 1
		error_messages += "\n\tParser: parseLine, No Hits Found Line."
	tests_count += 1
	#### isEOF
	#### createNextToken

	#***************************
	# End of Tests!
	#***************************

	print ""
	if fail_count > 0:
		print " [" + str(fail_count) + " of " + str(tests_count) + "] test(s) failed."
		print ""
		print error_messages
	else:
		print " [" + str(tests_count) + "] tests worked correctly."
	print ""
