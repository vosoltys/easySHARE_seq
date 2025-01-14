#!/bin/bash

#$ -pe parallel 20
#$ -l h_vmem=8G
#$ -N placement
#$ -o ./mm10_placement
#$ -j y
#$ -S /bin/bash
#$ -cwd



file=$1;
echo $file

fbname=$(basename $file _R1_001.fastq.gz)
echo $fbname

/local/bin/bwa mem -C -t 10 /genome/gbdb/mm10/mm10.fa $file ${file/R1_001/R2_001} -R "@RG\tID:$fbname\tSM:$fbname\tLB:SHARE_libraries\tPL:Illumina.NovaSeq.2x150" |
/local/bin/samtools view -bh - > ./$fbname.mm10.bam
/local/bin/samtools sort -@ 20 -l 9 -T ./$fbname.tmpsort -o ./$fbname.sorted.mm10.bam ./$fbname.mm10.bam
/local/bin/samtools index ./$fbname.sorted.mm10.bam
