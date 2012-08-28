#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-08-23

"""
usage:
keep_hsapiens.py filename.fa
"""

class Printer:
	def __init__ (self):
		self.printing = False

	def parseLine(self, line):
		if '>' in line:
			if "Homo" in line:
				self.printing = True
				print line.strip()
			else:
				self.printing = False
		else:
			if self.printing == True:
				print line.strip()
				

import sys

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print __doc__
		sys.exit(1)

	printer = Printer()
	filename=sys.argv[1]
	for line in open(filename):		
		printer.parseLine(line)
