#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-09-25

import sys

class StatsHolder:
	def __init__(self):
		self.reset()

	def reset(self):
		self.name = ""
		self.mature_name = ""
		self.strand = ""
		self.uridylation_count = 0
		self.adenylation_count = 0
		self.other_count = 0
		self.none_count = 0

	def print_current(self):
		t = "\t"
		print self.name + t + self.mature_name + t + self.strand + t + str(self.uridylation_count) + t + str(self.adenylation_count) + t + str(self.other_count) + t + str(self.none_count)

	def get_name(self):
		return self.name

stats_holder = StatsHolder()

for line in sys.stdin:
	name = line.split()[0]
	mature_name = line.split()[1]
	if stats_holder.name != "" and (stats_holder.name != name or stats_holder.mature_name != mature_name):
		stats_holder.print_current()
		stats_holder.reset()
	
	stats_holder.name = line.split()[0]
	stats_holder.mature_name = line.split()[1]
	stats_holder.strand = line.split()[2]
	stats_holder.uridylation_count += int(line.split()[3])
	stats_holder.adenylation_count += int(line.split()[4])
	stats_holder.other_count += int(line.split()[5])
	stats_holder.none_count += int(line.split()[6])

stats_holder.print_current()
