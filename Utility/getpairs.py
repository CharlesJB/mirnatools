#!/usr/bin/python
# encoding: utf-8
# author: Sébastien Boisvert
# Modified: Charles Joly Beauparlant

import sys

file_mature=sys.argv[1]
file_hairpins=sys.argv[2]
file_premature=sys.argv[3]

matures={}
for line in open(file_mature):
	i=line.split("\t")
	matures[int(i[0])]=i[1].strip()

hairpins={}
for line in open(file_hairpins):
	i=line.split("\t")
	hairpins[int(i[0])]=i[2].strip()

for line in open(file_premature):
	i=line.split("\t")
	hairpin=hairpins[int(i[0])]
	mature=matures[int(i[1])]
	print hairpin+"\t"+mature
