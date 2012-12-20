#!/bin/bash

# Note: the Usage function must be defined in the script calling ValidateFile,
#       before calling ValidateFile!
ValidateFile() {
	fileToCheck=$1
	if [ ! -e "$fileToCheck" ]
	then
		echo ""
		echo "Invalid file: $fileToCheck"
		Usage
		exit
	fi
}

