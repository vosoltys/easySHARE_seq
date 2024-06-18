#!/bin/bash

#$ -pe parallel 2
#$ -l h_vmem=4G
#$ -N Index1Split
#$ -o ./Read1.split.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/vsoltys" ];then mkdir /tmp/vsoltys;fi

main_dir=$1
var=$2
fbname=$3

for i in $var; do pi=`printf %02d $i`; for fastq in $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R1_001.fastq.gz.cutadapt"; do stem=`basename $fastq _001.fastq.cutadapt.gz`; cat $fastq | awk '/BX/ {match($0,"BX:Z"); d=substr($0,RSTART+13,4)};{print $0 | " gzip -c > '$main_dir'/cutadapt/'$stem'_split_"d"_001.fastq.gz"}' ;done;done;
