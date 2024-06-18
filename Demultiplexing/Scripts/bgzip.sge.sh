#!/bin/bash

#$ -pe parallel 25
#$ -l h_vmem=4G
#$ -N bgzip
#$ -o ./bgzip.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/vsoltys" ];then mkdir /tmp/vsoltys;fi

file=$1;
echo $file

bgzip -@25 $file
