#!/usr/bin/bash

# Script to execute Pipeline to analyze NGS data from YFV whole genome sequencing
# Needs to be in a folder along with reference fasta file and raw fastq files from a sample

#Create folder structure
mkdir fastq reference mapped passedQC variants assembly human_mapped

# Take sample file name
FILE1=$(basename *R1* .fastq.gz)
FILE2=$(basename *R2* .fastq.gz)
# Generalize this
FILE=$(basename *R1* _R1_001.fastq.gz)
# Reference variable
# Very important!!!!
# Specify the name of your reference of choice inside the quotes!!!
REFERENCE="YFV_genome_referenceRJ97.fasta"

# Move Files to suitable place
mv *.fastq.gz fastq/
mv *.fasta reference/

# Quality filtering
cutadapt -q 30,30 -m 50 --max-n 0 -u 9 -u -5 -U 9 -U -5 -o passedQC/${FILE1}_cutadapt.fastq -p passedQC/${FILE2}_cutadapt.fastq fastq/${FILE1}.fastq.gz fastq/${FILE2}.fastq.gz > passedQC/summary_cutadapt.log


# Indexing reference
bwa index reference/${REFERENCE}

# Alignment
bwa mem -B 2 -O 3 -t 12 -M reference/${REFERENCE} passedQC/${FILE1}_cutadapt.fastq passedQC/${FILE2}_cutadapt.fastq | samtools view -b > mapped/${FILE}_mapped_${REFERENCE}.bam

# Cleaning alignments
samtools view -F 2304 -f 2 -q 30 -b mapped/${FILE}_mapped_${REFERENCE}.bam > mapped/${FILE}_mapped_${REFERENCE}_clean.bam

# Ordering
samtools sort mapped/${FILE}_mapped_${REFERENCE}_clean.bam > mapped/${FILE}_mapped_${REFERENCE}_clean_sorted.bam
samtools sort mapped/${FILE}_mapped_${REFERENCE}.bam > mapped/${FILE}_mapped_${REFERENCE}_sorted.bam

# Indexing BAM
samtools index mapped/${FILE}_mapped_${REFERENCE}_clean_sorted.bam
samtools index mapped/${FILE}_mapped_${REFERENCE}_sorted.bam

# Calling variants 
freebayes -f reference/${REFERENCE} -b mapped/${FILE}_mapped_${REFERENCE}_clean_sorted.bam -v variants/variants_${FILE}_${REFERENCE}.vcf -p 4 -F 0.1








