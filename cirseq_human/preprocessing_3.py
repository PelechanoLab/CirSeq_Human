import sys,gzip,pysam

workdir = sys.argv[1]

infile = pysam.AlignmentFile(workdir + "/9_alignment.bam","rb")
outfile = pysam.AlignmentFile(workdir + "/10_alignment.bam","wb",template=infile)

def CompareAS(infile,outfile):	
	#running lists
	Alignment = []
	AlignmentScore = []
	
	#current read
	SequenceID = ""
	
	for read in infile.fetch(until_eof=True):
		
		#add alignments of the same read to the running list
		if read.query_name == SequenceID:
			#add only alignments lacking gaps and clipped bases
			cigar=read.cigarstring
			if cigar.count("D") == 0 and cigar.count("I") == 0 and cigar.count("S") == 0:
				AS = read.get_tag("AS") 
				AlignmentScore.append(AS)
				Alignment.append(read)
				
		else:
			if len(Alignment) > 0:
				#write alignment with the best alignment score
				outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])
			
			SequenceID = read.query_name
			
			Alignment = []
			AlignmentScore = []
			cigar=read.cigarstring
			if cigar.count("D") == 0 and cigar.count("I") == 0 and cigar.count("S") == 0:
				AS = read.get_tag("AS")
				AlignmentScore.append(AS)
				Alignment.append(read)



	if len(Alignment) > 0:
		#write alignment with the best alignment score
		outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])

CompareAS(infile,outfile)
infile.close()
outfile.close()
