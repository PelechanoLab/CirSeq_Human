from ConsensusModule import *
import sys
import os
import gzip
import numpy

#check args
if len(sys.argv) < 2:
	sys.stderr.write("usage: %s output_directory file_1.fastq file_2.fastq ... file_n.fastq\n" % (sys.argv[0]))
	sys.exit(1)

#validate input files
infiles = sys.argv[2:]
workdir = sys.argv[1]

if not os.path.isdir(workdir):
	sys.stderr.write("Workdir '%s' doesn't exist, aborting\n" % (workdir))
	sys.exit(1)

for f in infiles:
	if not os.path.isfile(f):
		sys.stderr.write("Input files '%s' does not exist, aborting\n" % (f))
		sys.exit(1)

outfile = gzip.open(workdir + "/1_consensus.fastq.gz","wb")
parameterfile1 = open(workdir + "/ProcessingStats.txt","w")
parameterfile2 = open(workdir + "/RepeatLengthDistribution.txt","w")
parameterfile3 = open(workdir + "/RepeatCopyDistribution.txt","w")
parameterfile4 = open(workdir + "/Length_Copy_Distribution.txt","w")
PoorQuality = 0
NoRepeats = 0
ShorterRepeatLength = 0
LongerRepeatLength = 0
LowIdentity = 0
ConsensusSequences = 0
TotalReads = 0

RepeatLengths = [0]*110
RepeatCopies = [0]*15
Length_Copy_Matrix = [[0]*110]*15

#read through fastq.gz files, generate consensus sequences using ConsensusModule and report read parameters
for f in infiles:
	sys.stderr.write("Processing file '%s'\n" % f)
	infile = gzip.open(f, 'rb')
	counter_PoorQuality, counter_NoRepeats, counter_ShorterRepeatLength, counter_LongerRepeatLength, counter_LowIdentity, counter_ConsensusSequences, counter_TotalReads, counter_RepeatLengths,counter_Copies, matrix_Len_Copies = Consensus(infile, outfile)
	infile.close()
	
	PoorQuality += counter_PoorQuality
	NoRepeats += counter_NoRepeats
	ShorterRepeatLength += counter_ShorterRepeatLength
	LongerRepeatLength += counter_LongerRepeatLength
	LowIdentity += counter_LowIdentity
	ConsensusSequences += counter_ConsensusSequences
	TotalReads += counter_TotalReads

	RepeatLengths = numpy.add(RepeatLengths, counter_RepeatLengths)
	RepeatCopies = numpy.add(RepeatCopies, counter_Copies)
	Length_Copy_Matrix = numpy.add(Length_Copy_Matrix, matrix_Len_Copies)

#write read parameters
parameterfile1.write("Consensus Generation Stats\n")
parameterfile1.write(str(PoorQuality) + "\tPoor quality\n")
parameterfile1.write(str(NoRepeats) + "\tNo repeats detected\n")
parameterfile1.write(str(ShorterRepeatLength) + "\tShort repeat length (<25nt)\n")
parameterfile1.write(str(LongerRepeatLength) + "\tAbnormal repeat length\n")
parameterfile1.write(str(LowIdentity) + "\tLow sequence identity between repeats\n")
parameterfile1.write(str(ConsensusSequences) + "\tConsensus sequences\n")
parameterfile1.write(str(TotalReads) + "\tTotal reads\n" + "\n")
parameterfile1.write(str(round((ConsensusSequences*100.0)/TotalReads, 2)) + "\t% Consensus sequence generated\n")

#write lengths of consensus reads
parameterfile2.write("Length\tNumberOfReads\n")
i = 0
for Length in RepeatLengths:
	parameterfile2.write(str(i) + "\t" + str(Length) + "\n")
	i += 1

parameterfile3.write("Copy\tNumberOfReads\n")
i = 0
for Copy in RepeatCopies:
        parameterfile3.write(str(i) + "\t" + str(Copy) + "\n")
        i += 1
	
parameterfile4.write("Copy\tLength\tNumberOfReads\n")
i = 0
for Copy in Length_Copy_Matrix:
	j = 0
	for Length in Copy:
		parameterfile4.write(str(i) + "\t" + str(j) + "\t" + str(Length) + "\n")
		j += 1
	i += 1

outfile.close()
parameterfile1.close()
parameterfile2.close()
parameterfile3.close()
parameterfile4.close()
