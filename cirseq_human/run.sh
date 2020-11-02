#!/bin/bash

#$ -l h_rt=336:0:0
#$ -cwd
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -N cirseq

set -e

WORKDIR=$1
REFDIR_IDX=$2
THREAD=$3
SCRIPTDIR=$4
SPLICESITE=$5
FASTQS="${@:6}"

#check if working directory exists
mkdir -p $WORKDIR

#check if script directory exists
if [ ! -d "$SCRIPTDIR" ];then 
	echo "Directory for scripts does not exist!"
	exit 0
fi

# generate consensus
python ${SCRIPTDIR}/ConsensusGeneration.py $WORKDIR $FASTQS

# align consensus
echo "Aligning all consensus sequences..."
hisat2 -p $THREAD -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/1_consensus.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/2_alignment.bam

# preprocess 1
echo "Extracting perfectly mapped reads and rearranging/rotating the rest..."
python ${SCRIPTDIR}/preprocessing_1.py ${WORKDIR}

# align again
echo "Aligning the rearranged sequences..."
hisat2 -p $THREAD --no-unal -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/4_rearranged.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/6_alignment.bam

# preprocess 2
echo "Extracting perfectly mapped reads and rotating the rest..."
python ${SCRIPTDIR}/preprocessing_2.py ${WORKDIR} 

echo "Aligning the rotated sequences..."
hisat2 -p $THREAD --no-unal -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/5_rotated.fastq.gz,${WORKDIR}/8_rotated.fastq.gz | samtools view -@ $THREAD -b -o ${WORKDIR}/9_alignment.bam

# preprocess 3
echo "Extracting the best hit of each perfectly mapped sequence..."
python ${SCRIPTDIR}/preprocessing_3.py ${WORKDIR} 

echo "Combining alignment files and sorting..."
samtools merge -f -@ $THREAD ${WORKDIR}/consensus_alignment.bam ${WORKDIR}/3_alignment.bam ${WORKDIR}/7_alignment.bam ${WORKDIR}/10_alignment.bam
samtools sort -@ $THREAD -o ${WORKDIR}/consensus_align_sorted.bam ${WORKDIR}/consensus_alignment.bam
rm -f ${WORKDIR}/4_rearranged.fastq.gz ${WORKDIR}/*_rotated.fastq.gz ${WORKDIR}/*_alignment.bam
samtools index -@ $THREAD ${WORKDIR}/consensus_align_sorted.bam
