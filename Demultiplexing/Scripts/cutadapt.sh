#!/bin/bash

#$ -pe parallel 20
#$ -l h_vmem=8G
#$ -N cutadapt
#$ -o ./cutadapt.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/vsoltys" ];then mkdir /tmp/vsoltys;fi

file=$1;
file2=$2
echo $file
echo $file2

/fml/chones/local/bin/cutadapt -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG --nextseq-trim 20 --cores=20 -O 5 -m 30 -o $file.cutadapt -p $file2.cutadapt --too-short-output=$file.tooshort.fastq.gz --too-short-paired-output=$file2.tooshort.fastq.gz $file $file2
