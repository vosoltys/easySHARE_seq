#!/bin/bash
#$ -pe parallel 2
#$ -l h_vmem=12G
#$ -N ValReads_extract
#$ -o ./ValidReads_extract.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V
#Script for extracting the valid reads from RNAseq fastqs and 
# pi= N5 index (e.g 002) and the reads are the full read names (or absolute path)
pi=$1
R1=$2
R2=$3
Pool=$4

zcat $R2 | /fml/chones/local/bin/bioawk -c fastx seq | awk '{if (($2 ~ /^TTT[AGCT]T/) || ($2 ~ /^[AGCT]TTTT/) || ($2 ~ /^T[AGCT]TTT/) || ($2 ~ /^TT[AGCT]TT/) || ($2 ~ /^TTTT[AGCT]/) || ($2 ~ /^[AGCT]TTTT/)) print $0}' | awk '{print $1}' > "$Pool"_"$pi"_Valid_Reads; cat "$Pool"_"$pi"_Valid_Reads | split -a 2 -d -l 5000000 - "$Pool"_Valid_Reads_"$pi"_

#x is the number of files each ValidReads file gets split up into 
for i in `seq 0 30`; do x=`printf %02d $i`; zcat $R1 | grep -a -A3 -F -f "$Pool"_Valid_Reads_"$pi"_"$x" | awk '{if ($0 != "--") print $0}' | gzip -c - >> "$R1".filtered.fastq.gz; done

#As a last step, move the barcode of Read1 into the readname
zcat "$R1".filtered.fastq.gz | /fml/chones/local/bin/bioawk -c fastx seq | awk '{print "@"$1"_"substr($4,RSTART+6,12) "\n" $2 "\n" "+" "\n" $3}' | gzip -c > "$pi"_"$Pool"_R1_UMI_001.fastq.gz
