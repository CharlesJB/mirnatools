#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-09-27

"""

Usage:
	fetchStats.py <file_pairs> <file_mature> <file_blast> <file_collapsed> <file_premature>

"""

# Classes definitions
class MatureAllocator:
	def __init__(self, file_pairs):
		self.reset()
		self._load_file(file_pairs)

	def reset(self):
		self.m_pair_list = {}

	def get_mature_miRNA(self, miRNA_name):
		if miRNA_name in self.m_pair_list:
			return self.m_pair_list[miRNA_name]
		else:
			return ""

	def _load_file(self, file_pairs):
		for line in open(file_pairs):
			tokens = line.split()
			genomic = tokens[0]
			mature = tokens[1]
			if genomic not in self.m_pair_list:
				self.m_pair_list[genomic] = []
			self.m_pair_list[genomic].append(mature)

class MatureSequence:
	def __init__(self, mature_sequence):
		self.m_sequence = mature_sequence
		self.m_count = 0

	def set_count(self, count):
		self.m_count = count
	
	def get_count(self):
		return self.m_count

	def get_sequence(self):
		return self.m_sequence

def get_strand(premature_seq, subject_start): # TODO: test this
	if subject_start < (len(premature_seq)/2):
		return "5p"
	return "3p"

def update_return(last, new):
	if last == "O":
		return "O"
	if last == "" and new == "U":
		return "U"
	if last == "U" and new == "U":
		return "U"
	if last == "U" and new == "A":
		return "O"
	if last == "" and new == "A":
		return "A"
	if last == "A" and new == "A":
		return "A"
	if last == "A" and new == "U":
		return "O"

def check_alternate_base(sequence, start, end):
	toReturn = ""
	if start != 1:
		for i in range(0,start-1):
			toReturn = update_return(toReturn, sequence[i])
	if end != len(sequence):
		for i in range(end, len(sequence)):
			toReturn = update_return(toReturn, sequence[i])
	return toReturn
	

import sys
if __name__=="__main__":
	if len(sys.argv) == 1:
		print __doc__
                sys.exit(1)

        if len(sys.argv)!=6:
                print __doc__
                sys.exit(1)

file_pairs = sys.argv[1]
file_mature = sys.argv[2]
file_blast= sys.argv[3]
file_collapsed = sys.argv[4]
file_premature= sys.argv[5]
mode = "stats"
#mode = sys.argv[6]

# 1. Load pairs
mature_allocator = MatureAllocator(file_pairs)

# 2. Load mature sequences in a {}
sequence_counts = {}
#sequence_id = {}
collapsed_sequences = {}
count = 0
identifier = ""
for line in open(file_collapsed):
	if '>' in line:
		count = int(line.split()[3])
		identifier = line.split()[1][3:]
		sequence_counts[identifier] = count
	else:
		sequence = line.strip()
		collapsed_sequences[identifier] = sequence
#	else:
#		sequence_counts[line.strip()] = count
#		sequence_id[identifier] = line.strip()

#for entry in sequence_counts:
#	print entry + ": " + str(sequence_counts[entry])

mature_sequences = {}
name = ""
sequence = ""
for line in open(file_mature):
	if '>' in line:
		name = line.split()[0][1:]
	else:
		if '-' not in line: # Some lines only have the value "--", we want to skip them
			sequence = line.strip()
			mature_sequences[name] = MatureSequence(sequence) # TODO: use a simple {} string string insead?
#			if sequence in sequence_counts:
#				mature_sequences[name].set_count(sequence_counts[sequence]) # TODO test this
#			else:
#				mature_sequences[name].set_count(0) # TODO test this

#for entry in mature_sequences:
#	print "--------------"
#	print entry + ": "
#	print "sequence: " + str(mature_sequences[entry].get_sequence())
#	print "count: " + str(mature_sequences[entry].get_count())

premature_sequences = {}
name = ""
sequence = ""
for line in open(file_premature):
	if '>' in line:
		name = line.strip().split()[0][1:]
	else:
		if '-' not in line: # Some lines only have the value "--", we want to skip them
			sequence = line.strip()
			premature_sequences[name] = sequence

#for entry in premature_sequences:
#	print entry + ": " + str(premature_sequences[entry])

# # 0	hsa-mir-103b-2	100.00	23	0	0	1	23	23	1	1e-08	46.1
# 3. Parse blast output
#3- Pour chaque lignes dans le fichier blast:
last_query = ""
toPrint = []
current_count = 0

# For counts.txt, check pre-normalize counts
#normalized_counts = {}
#total_entry = {}
#for line in open(file_blast):
#	entry = line.split()[0]
#	identifier = line.split()[1][3:]
#	if entry not in total_entry:
#		total_entry[entry] = 0
#	total_entry[entry] += 1

#for entry in total_entry:
#	normalized_counts[entry] = float(sequence_counts[identifier]) / float(total_entry[entry])

for line in open(file_blast):
	query = line.split()[0]
	pre_mature_name = line.split()[1]
#	mature_count = mature_sequences[mature_name].get_count() # TODO this is not ok, we need the query count!!!1
#	query_count = sequence_counts[sequence_id[query]]
	query_count = sequence_counts[query]
	percent_match = line.split()[2]
	align_length = int(line.split()[3])
	mismatch_count = int(line.split()[4])
	gap_count = int(line.split()[5])
	query_start =  int(line.split()[6])
	query_end=  int(line.split()[7])
	subject_start =  int(line.split()[8])
	subject_end=  int(line.split()[9])
	if query_start > query_end:
		tmp = query_start
		query_start = query_end
		query_end = tmp
	if subject_start > subject_end:
		tmp = subject_start
		subject_start = subject_end
		subject_end = tmp

	if last_query == "":
		last_query = query
	elif query != last_query:
		for line_to_print in toPrint:
			if line_to_print != "skip":
				if current_count != 0:
					if mode == "stats":
						print line_to_print + str(float(query_count)/float(current_count))
				else:
					if mode == "stats":
						print line_to_print + "0"
		last_query = query
		toPrint = []
		current_count = 0
	else:
		current_count += 1
	
#	toPrint = pre_mature_name + "\t" + mature_name + "\t" + strand + "\t"
		
#  3-1 Déterminer sur quel brin est le miRNA (tmp fix: $9 > 30; best:  checker la position dans la sequence)
	strand = get_strand(premature_sequences[pre_mature_name], subject_start) # TODO: fix this
#	print "---------"
#	print "strand: " + strand
#	print "pre_mature_name: " + pre_mature_name
	mature_names = mature_allocator.get_mature_miRNA(pre_mature_name)
#	print "len(mature_names): " + str(len(mature_names))
	mature_name = ""
	for name in mature_names:
		if strand in name:
			mature_name = name
	if mature_name == "":
		strand = "-"
		if len(mature_names) > 0:
			mature_name = mature_names[0]
		else:
			mature_name = "-"
#	print "mature_name: " + mature_name
	toPrint.append(pre_mature_name + "\t" + mature_name + "\t" + strand + "\t")
#  3-2 S'assurer que c'est max +-3 de chaque côté
	if mature_name in mature_sequences and mature_name != "-" and mature_name not in mature_sequences:
		mature_name = mature_name[:len(mature_name)-3]
	if mature_name in mature_sequences and mature_name != "-" and query_start < 4 and query_end > (len(mature_sequences[mature_name].get_sequence()) - 3):
#  3-3 Si pleine longueur et pas de mismatch -> none
		if query_start == 1 and query_end == align_length and percent_match == "100.00":
			toPrint[len(toPrint)-1] += "none" + "\t"
#			toPrint += "none" + "\t" + mature_count + "\n"
#  3-4 Si mismatch -> other
		elif query_start == 1 and query_end == align_length:
			toPrint[len(toPrint)-1] += "other" + "\t"
#			toPrint += "other" + "\t" + mature_count + "\n"
#  3-5 Sinon:
		else:
#    3-5-1 Checker les bases modifiées:
			alternate_base = check_alternate_base(mature_sequences[mature_name].get_sequence(), query_start, query_end)
#      3-5-1-1 Si U uniquement -> uridylation
			if alternate_base == "U":
				toPrint[len(toPrint)-1] += "uridylation" + "\t"
#				toPrint += "uridylation" + "\t" + mature_count + "\n"
#      3-5-1-2 Si A uniquement -> adenylation
			elif alternate_base == "A":
				toPrint[len(toPrint)-1] += "adenylation" + "\t"
#				toPrint += "adenylation" + "\t" + mature_count + "\n"
#      3-5-1-3 Si A et U -> other
#			elif alternate_base == "O":
			else:
				toPrint[len(toPrint)-1] += "other" + "\t"
#				toPrint[len(toPrint)-1] += "both" + "\t"
#				toPrint += "other" + "\t" + mature_count + "\n"
	else:
		toPrint[len(toPrint)-1] = "skip"
