#!/bin/bash

#$ -pe parallel 5
#$ -l h_vmem=8G
#$ -N Demux
#$ -o ./SHARE_cluster_demult.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/vsoltys" ];then mkdir /tmp/vsoltys;fi

input_file=$1;
output_file=$2;
main_dir=$3;
echo $input_file
echo $output_file
cp $main_dir/Scripts/BC_A.txt ./
cp $main_dir/Scripts/BC_C.txt ./
cp $main_dir/Scripts/BC_D.txt ./
cp $main_dir/Scripts/SHARE_demult_fastq.o ./

./SHARE_demult_fastq.o $input_file $output_file
