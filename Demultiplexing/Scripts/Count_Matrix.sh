#!/bin/bash

#$ -pe parallel 10
#$ -l h_vmem=32G
#$ -N umi_count
#$ -o ./umi_dedup.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V

#Script for deduplicating bamfile and counting UMIs, generating a matrix

file=$1
num=$2

/fml/chones/local/bin/umi_tools count --extract-umi-method=tag --cell-tag=CR --umi-tag=BX --per-gene --gene-tag=XT --per-cell --assigned-status-tag=XS --stdin=$file --log=$num.log  --stdout=$num.csv
