#!/bin/bash

#$ -pe parallel 2
#$ -l h_vmem=8G
#$ -N Subsample
#$ -o ./SubSam.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V

file=$1;
echo $file

int=$2
echo $int

/fml/chones/local/bin/samtools view -@ 10 -bh -s 4225.$int $file > ./$int.subs.bam
samtools index -@ 10 $int.subs.bam
