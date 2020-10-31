#!/bin/bash

#$ -l h_rt=336:0:0
#$ -cwd
#$ -l arch=linux-x64
#$ -S /bin/bash
#$ -N cirseq

set -e

WORKDIR=$1
REFFILE=$2
THREAD=$3
SCRIPTDIR=$4
SPLICESITE=$5
FASTQS="${@:6}"

#check if working directory exists
mkdir -p $WORKDIR

# generate consensus
python ${SCRIPTDIR}/ConsensusGeneration.py $WORKDIR $FASTQS

# align consensus
sat2 -p $THREAD --no-hd --no-sq -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/1_consensus.fastq.gz | gzip -c > ${WORKDIR}/2_alignment.sam.gz

# preprocess 1
python ${SCRIPTDIR}/preprocessing_1.py ${WORKDIR}

# align again
hisat2 -p $THREAD --no-unal --no-hd --no-sq -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/4_rearranged.fastq.gz | gzip -c > ${WORKDIR}/6_alignment.sam.gz

# preprocess 2
python ${SCRIPTDIR}/preprocessing_2.py ${WORKDIR} 

hisat2 -p $THREAD --no-unal --no-hd --no-sq -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/5_rotated.fastq.gz | gzip -c > ${WORKDIR}/9_alignment.sam.gz

hisat2 -p $THREAD --no-unal --no-hd --no-sq -x ${REFDIR_IDX} --known-splicesite-infile ${SPLICESITE} -U ${WORKDIR}/8_rotated.fastq.gz | gzip -c > ${WORKDIR}/10_alignment.sam.gz

# preprocess 3
python ${SCRIPTDIR}/preprocessing_3.py ${WORKDIR} 

cat ${WORKDIR}/3_alignment.sam.gz ${WORKDIR}/7_alignment.sam.gz ${WORKDIR}/11_alignment.sam.gz > ${WORKDIR}/data.sam.gz
rm -f ${WORKDIR}/3_alignment.sam.gz ${WORKDIR}/7_alignment.sam.gz ${WORKDIR}/11_alignment.sam.gz ${WORKDIR}/2_alignment.sam.gz ${WORKDIR}/4_rearranged.fastq.gz ${WORKDIR}/5_rotated.fastq.gz ${WORKDIR}/6_alignment.sam.gz ${WORKDIR}/8_rotated.fastq.gz ${WORKDIR}/9_alignment.sam.gz ${WORKDIR}/10_alignment.sam.gz

