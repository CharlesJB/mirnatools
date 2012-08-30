#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2012-07-23

"""
This script calculates the distribution of sequence's length from a fasta file
Usage:
cat joe.fasta | ./LengthDistribution.py <maxLength> <NbSeq> > distribution.txt
maxLength: Maximum sequence length to compute
NbSeq: Optional parameter. Will compare to that value when calculating percentage.
"""

class DistributionCalculator:
        def __init__(self, maxLength, total):
                self.maxLength = maxLength
		self.total = total
                self.clear()

        def clear(self):
                self.distribution = [0] * self.maxLength

        def processSequence(self, sequence, count):
                length = len(sequence)
                if length <= self.maxLength:
                        self.distribution[length] += count

        def printResults(self):
		# Calculate total number of sequences (if not specified by user)
		if len(self.total) == 0:
			total = 0
			for i in range(self.maxLength):
				total += self.distribution[i]
		else:
			total = self.total
		# Print percentage of total for each lenght
                for i in range(self.maxLength):
			percent = self.distribution[i] / float(total) * 100
                        print str(i) + "\t" + str(percent)
	

import sys

if __name__=="__main__":
        if len(sys.argv)!=2 or len(sys.argv)!=3:
                print __doc__
                sys.exit(1)

        maxLength=int(sys.argv[1])
	total=int(sys.argv[2])
        distributionCalculator = DistributionCalculator(maxLength, total)

        i =0
        for line in sys.stdin:
		if i==0:
			count = int(line.split()[3])
                if i==1:
                        sequence = line.strip()
                        distributionCalculator.processSequence(sequence, count)

                i+=1

                if i==2:
                        i=0

        distributionCalculator.printResults()

