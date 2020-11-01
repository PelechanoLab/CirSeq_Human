import sys,gzip,pysam
from pysam import qualities_to_qualitystring as decode

workdir = sys.argv[1]

infile = pysam.AlignmentFile(workdir + "/2_alignment.bam","rb")
outfile1 = pysam.AlignmentFile(workdir + "/3_alignment.bam","wb",template=infile)
outfile2 = gzip.open(workdir + "/4_rearranged.fastq.gz","wb") #clipped reads with order of sequence blocks swapped
outfile3 = gzip.open(workdir + "/5_rotated.fastq.gz","wb") # every possible rotation of gapped and multi-clipped reads

def Rotate(read):
	i = 0
	while i < read.query_length:
		outfile3.write("@" + read.query_name + "\n")
		outfile3.write(read.query_sequence[-i:] + read.query_sequence[:-i] + "\n")
		outfile3.write("+" + "\n")
		outfile3.write(decode(read.query_qualities[-i:]) + decode(read.query_qualities[:-i]) + "\n")
		i += 1

for read in infile.fetch(until_eof=True):
		
	#make every possible rotation of reads with gaps and more than one clipped sequence
	cigar=read.cigarstring
	if read.is_unmapped:
		Rotate(read)
	elif (cigar.count("D")> 0 or cigar.count("I") > 0 or cigar.count("S") > 1) and read.is_secondary == False:
		Rotate(read)

	elif cigar.count("S") == 0:
		#fetch number of mismatches, save perfectly matched and ungapped alignments
		if read.get_tag("NM") == 0:
			outfile1.write(read)
		#make every possible rotation of reads with no gaps but have mismatches. Coincidentally matched bases near the
		#ends may suppress sequence clipping
		elif read.is_secondary == False:
			Rotate(read)
		
	#rearrange sequences with a single block of clipped reads
	elif cigar.count("S") == 1 and read.is_secondary == False:
		
		#parse CIGAR
		numbers = ""
		letters = []
		for character in cigar:
			if character == "M" or character == "S":
				letters.append(character)
				numbers += "X"
			else:
				numbers += character
		numbers = numbers.split("X")
		
		#rearrange seuences and quality scores
		Sequence = read.query_sequence[int(numbers[0]):] + read.query_sequence[0:int(numbers[0])]
		QualityScores = decode(read.query_qualities[int(numbers[0]):]) + decode(read.query_qualities[0:int(numbers[0])])
		outfile2.write("@" + read.query_name + "\n")
		outfile2.write(Sequence + "\n")
		outfile2.write("+" + "\n")
		outfile2.write(QualityScores + "\n")

infile.close()
outfile1.close()
outfile2.close()
outfile3.close()
