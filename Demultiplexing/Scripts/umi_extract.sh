#!/bin/bash

#$ -pe parallel 2
#$ -l h_vmem=12G
#$ -N umi_extract
#$ -o ./umi_extract.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V

#Script for extracting the UMI from RNAseq fastqs
R1=$1;
R2=$2;

/fml/chones/local/bin/umi_tools extract --bc-pattern=NNNNNNNNNN --stdin $R2 --stdout $R2.UMI.fastq.gz --read2-in $R1 --read2-out $R1.UMI.fastq.gz -L extract.log
