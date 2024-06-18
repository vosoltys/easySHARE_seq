#!/bin/bash

#$ -pe parallel 5
#$ -l h_vmem=4G
#$ -N split
#$ -o ./split.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd

input=$1
output=$2
echo "Input file is: "$1
echo "Output file is: "$2


zcat $1 | split -a 2 -d -l 200000000 - $2
