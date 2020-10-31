CirSeq
======
Description
---------------------------------------------------------------------------
CirSeq is an application that generates consensus sequences from tandem 
repeat sequence reads. This application accepts sequence data as gzipped 
FASTQ files (extension .fastq.gz) and requires a user-supplied reference 
sequence in FASTA format (extension .fasta). The output generated by CirSeq
includes the following files:

1. consensus.fastq.gz

	Contains all consensus sequences generated in FASTQ format. Per 
	base quality scores can be converted to estimated error 
	probabilities using the formula: 10**(3*Quality score/-10).

2. data.sam.gz

	Contains all consensus sequences that map to the user-supplied 
	reference sequence. Data is in SAM format.

3. Q[positive integer <= 41]threshold.txt

	Summarizes the counts of bases mapped to each position of the 
	reference sequence at and above a quality threshold. The default 
	quality threshold is 20 (estimated error probability of 10**-6). 
	The columns in this file are as follows: Numeric position 
	(according to reference sequence provided), the reference base 
	(wild type sequence according to reference sequence provided) and 
	four columns with the counts of each base, A, C, G and T, 
	respectively, observed at each reference position.

4. RepeatLengthDistribution.txt

	Contains the counts of reads with repeats from 25 - 99 bases long. 
	This distribution can aid in diagnosis of issues with sequencing 
	library preparation. The columns in this file are as follows: 
	Repeat length (in bases), counts of reads with the indicated 
	length.

5. QualityMetrics.txt

	Summarizes the observed mutation frequency and 
	transition:transversion rates with respect to quality scores. These 
	metrics can be used to define an appropriate quality threshold. The 
	columns in this file are as follows: Quality score, mismatches, 
	aligned bases, transition mutations and transversion mutations.

6. ProcessingStats.txt

	Summarizes important statistics of the sequencing data processing 
	that may aid in evaluating the quality of sequencing libraries and 
	diagnosing problems with library preparation and sequencing.


System requirements
---------------------------------------------------------------------------
The following packages are prerequisites for running CirSeq:

1. Python (version 2.7.5)
2. Cython (version 0.19.1)
3. NumPy (version 1.7.1)
4. SciPy (version 0.13.3)
5. Bowtie2 (version 2.1.0)
6. R (version 3.0.1) OPTIONAL for generation of plotted outputs

NOTE 1: Cython requires a compiler. For OSX this may require installation 
of Xcode.

NOTE 2: Bowtie2 and R must be in the PATH.


Setup
---------------------------------------------------------------------------
Execute the following command in the script directory to compile the 
ConsensusGeneration module:

python setup.py build_ext —inplace


Usage
---------------------------------------------------------------------------
.[Script directory]/run.sh [output directory] [reference file (FASTA format)] [Script directory] [gzipped FASTQ file(s)]

