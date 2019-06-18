#!/usr/bin/bash

# Script to assembly YFV genomes with SPADES
# Needs to be in a folder along with reference fasta file and raw fastq files from a sample

#Create folder structure
mkdir -p fastq reference mapped passedQC variants assembly human_mapped

# Take sample file name
SAMPLE1=$(basename fastq/*R1* .fastq.gz)
SAMPLE2=$(basename fastq/*R2* .fastq.gz)
# Generalize this
SAMPLE=$(basename fastq/*R1* _R1_001.fastq.gz)
# Reference variable
# Very important!!!!
# Specify the name of your reference of choice inside the quotes!!!
REFERENCE_HUMAN="/stgisilon_reference/basespace_hg19/Homo_sapiens/UCSC/hg19/FASTA_NC_012920_1/genome.fa"
REFERENCE_YFV="reference/YFV_RJ104.fasta"

# Align to human genome
bwa mem -t 12 -M $REFERENCE_HUMAN passedQC/${SAMPLE1}_cutadapt.fastq passedQC/${SAMPLE2}_cutadapt.fastq | samtools view -b > human_mapped/${SAMPLE}_mapped_human.bam

# Save unmapped reads
samtools view -f 12 human_mapped/${SAMPLE}_mapped_human.bam > human_mapped/${SAMPLE}_mapped_human_unmapped.bam

# Convert to FASTQ
bamToFastq -i human_mapped/${SAMPLE}_mapped_human_unmapped.bam -fq human_mapped/${SAMPLE}_unmapped_human_R1.fastq -fq2 human_mapped/${SAMPLE}_unmapped_human_R2.fastq

# Assembly
spades.py -1 human_mapped/${SAMPLE}_unmapped_human_R1.fastq -2 human_mapped/${SAMPLE}_unmapped_human_R2.fastq --trusted-contigs $REFERENCE_YFV -t 12 -m 20 -o assembly/








