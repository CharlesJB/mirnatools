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

	def reset(self):
		m_speciesList = {}

	def loadDataStore(self, pathToDataStore):
		for line in open(pathToDataStore):
			binomialName = line.strip().split()[0].lower()
			commonName = line.strip().split()[1]
			if binomial not in m_speciesList:
				m_speciesList[binomialName] = commonName

	def loadAlternativeDataStore(self, pathToDataStore):
		self.reset()
		self.loadDataStore(pathToDataStore)

	def convertSpecieName(self, binomialName):
		return m_speciesList[binomialName.lower()]
