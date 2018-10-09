#!/usr/bin/python
#
# @author Sergey Menis
# @email sergey.menis@gmail.com
#
# Identify all common nmers in two FASTA format files
#

import sys

# Hamming Distance Comparer
# Default is distance of ZERO
# Modified from: http://en.wikipedia.org/wiki/Hamming_distance
def hammingDistanceComparer(s1, s2, maxMismatchScore=None):
    if maxMismatchScore is None:
	maxMismatchScore = 0
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    hammingDistance =sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
    if hammingDistance <= maxMismatchScore:
	return 1
    else:
	return 0
	
# Blossum Match - The score must not be 
#def BlossumComparer(StringA, StringB, minimumMatchScore=None):
#	if minimumMatchScore is None:
#		minimumMatchScore = 0

# Read FASTA 
def readFastaEntry( fp ):
	seq = []
	#Drop the first line
	title = fp.readline()
    	while 1:
		line = fp.readline().rstrip()
		if line == "":
            		break
        	seq.append(line)
	sequence = "".join(seq)
	return sequence

def generateNmers(sequence, size):
	nmers = []
	beg = 0
	end = 0 + size
	while end <= len(sequence):
		nmers.append((sequence[beg:end], beg+1, end+1))
		beg = beg + 1
		end = beg + size	
	return nmers

def removeSubranges(records):
	observedRanges = []
	cleanRecords = []
	for record in records:
		found = 0
		for observed in observedRanges:
			if ((record[1] >= observed[0]) and (record[2] <= observed[1]) and (record[4] >= observed[2]) and (record[5] <= observed[3])):
				found = 1
		if found == 0:
			cleanRecords.append((record[0], record[1], record[2], record[3], record[4], record[5]))
			observedRanges.append((record[1], record[2], record[4], record[5]))	

	return cleanRecords
	

def main():
	
	if len(sys.argv) < 2:
		print "Usage: findCommondNmers.py -f1 firstFASTA -f2 secondFasta"
		print "       (opt)-from 10 (opt)-to 20 (opt)-minMatchPercentage 100"
		sys.exit()

	if '-f1' in sys.argv:
		fileOne = sys.argv[sys.argv.index('-f1') + 1]
	else:
		sys.exit("ERROR: You must provide a FASTA file for -f1")

        if '-f2' in sys.argv:
                fileTwo = sys.argv[sys.argv.index('-f2') + 1]
        else:
                sys.exit("ERROR: You must provide a FASTA file for -f2")	
	
	# Load the files from FASTA
	fp = open(fileOne)
	sequenceOne = readFastaEntry(fp)
	fp.close()
	
	fp = open(fileTwo)
	sequenceTwo = readFastaEntry(fp)
	fp.close()

        # Defaults
        fromNmer = 10
        toNmer = max(len(sequenceOne), len(sequenceTwo))
        minMatch = 100 

        if '-from' in sys.argv:
                fromNmer = int(sys.argv[sys.argv.index('-from') + 1]) 
        if '-to' in sys.argv:
                toNmer = int(sys.argv[sys.argv.index('-to') + 1])    

        if '-minMatchPercentage' in sys.argv:
                minMatch = int(sys.argv[sys.argv.index('-minMatchPercentage') + 1]) 




	# Generate all possible nmers and compare them 
	nmer = toNmer
	commonNmers = []
	while nmer >= fromNmer:
		nmersOne = generateNmers(sequenceOne, nmer)
		nmersTwo = generateNmers(sequenceTwo, nmer)
				
		for r1 in nmersOne:
			for r2 in nmersTwo:
				 if hammingDistanceComparer(r1[0], r2[0], int((nmer*(100-minMatch))/100)) == 1:
					commonNmers.append((r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]))
		nmer = nmer - 1

	# Clean up duplicates: subranges of larger strings
	rangesToPrint = removeSubranges(commonNmers)

	# Print the result
	for r in rangesToPrint:
		print r[0], r[1], r[2], r[3], r[4], r[5]


if __name__ == "__main__":
    main()
