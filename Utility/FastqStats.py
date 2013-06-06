#!/usr/bin/python
# encoding: utf-8
# author: Charles Joly Beauparlant
# 2013-06-01

"""

This will count the number of sequences in a fastq file, append the results to a file 
and redirect the sequences to stdout.

Usage:

cat filename.fastq | FastqStats.py <filename> <message>
	message: A line that describe the current status of the fastq stream.
	filename: The name of the output file.

"""

VERSION="1.1"

import sys

if __name__=="__main__":
	if len(sys.argv) == 1:
		print __doc__
                sys.exit(1)

	if sys.argv[1] == "version":
		print "FastqStats.py v." + str(VERSION)
		sys.exit()

	if sys.argv[1] == "tests":
		Tests()
		sys.exit()

        if len(sys.argv)!=3:
                print __doc__
                sys.exit(1)

	filename = sys.argv[1]
	message = sys.argv[2]

	# Count the number of lines in the file
	count = 0
	for line in sys.stdin:
		sys.stdout.write(line)
		count += 1

	# Save the results to file
	try:
		with open(filename, "a") as output:
			output.write(message + "\t")
			if count > 0:
				output.write(str(count/4) + "\n")
			else:
				output.write("0\n")
	except IOError:
		sys.stderr.write("FastqStats.py: Invalid filename: " + filename)
		sys.stderr.write(message + "\t")
		sys.stderr.write(str(count/4) + "\n")
