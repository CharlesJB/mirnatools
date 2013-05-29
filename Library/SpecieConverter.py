#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-01-21

"""
Convert binomial specie name to commun usage. i.e.: Homo sapiens to human (this is a library)

DataStore is a plain text file with 2 columns. The first is the binomial name and the second the commun usage.

The default DataStore is DataStore/speciesConversion.ds
"""

class SpecieConverter:
	def __init__(self):
		self.reset()
		self.loadDataStore("DataStore/speciesConversion.ds")

	def reset(self):
		self.m_speciesList = {}

	def loadDataStore(self, pathToDataStore):
		for line in open(pathToDataStore):
			binomialName = ' '.join(line.strip().split()[0:2]).lower()
			commonName = ' '.join(line.strip().split()[2:])
			if binomialName not in self.m_speciesList:
				self.m_speciesList[binomialName] = commonName

	def loadAlternativeDataStore(self, pathToDataStore):
		self.reset()
		self.loadDataStore(pathToDataStore)

	def convertSpecieName(self, binomialName):
		if binomialName.lower() in self.m_speciesList:
			return self.m_speciesList[binomialName.lower()]
