#!/bin/bash

#$ -pe parallel 2
#$ -l h_vmem=8G
#$ -N Filter
#$ -o ./Filter.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V

file=$1;
echo $file

if [ ! -e "/tmp/vsoltys" ]
                         then mkdir /tmp/vsoltys
fi

samtools view -H $file > $file.header; samtools view -F 4 -q 30 $file | grep -v chrM | grep -v A000 | grep -v C000 > $file.sam; cat $file.header $file.sam | samtools view -bh - > $file.filtered.bam; samtools index $file.filtered.bam; samtools view $file.filtered.bam | wc -l > $file.count; echo $file | cut -c3-11 | sed 's/Pool-//g'| awk '{print "5"$0}' | paste - $file.count >> FilteredBams.Count; rm $file.count $file.sam $file.header
