#!/bin/bash

#$ -pe parallel 5
#$ -l h_vmem=8G
#$ -N Full_Demux
#$ -o ./Full_Demux_SHARE.out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd



#Requirements: in main directory folder /fastqs with the undetermined fastqs (4 of them!)
#usage: ./Full... BASENAME_INPUT_FILE MAIN_DIRECTORY
#example: ./Full_demultiplex_easySHARE.sh Test.Undetermined_S0_L002 my_favourite_dir/

input_file=$1;
echo "Input files are : "$input_file
main_dir=$2;
echo "Main directory is: "$main_dir
fbname=$(basename $input_file Undetermined_S0_L001)
echo "Current user is "$USER
#First is splitting input Undetermined fastqs into chunks of 2M

echo "Splitting files..."

mkdir $main_dir/split

cp $main_dir/Scripts/split.sh $main_dir
qsub $main_dir/split.sh $main_dir/fastqs/$fbname"_R1_001.fastq.gz" $main_dir/split/$fbname"_R1_001.fastq.gz_"
qsub $main_dir/split.sh $main_dir/fastqs/$fbname"_R2_001.fastq.gz" $main_dir/split/$fbname"_R2_001.fastq.gz_"
qsub $main_dir/split.sh $main_dir/fastqs/$fbname"_I1_001.fastq.gz" $main_dir/split/$fbname"_I1_001.fastq.gz_"
qsub $main_dir/split.sh $main_dir/fastqs/$fbname"_I2_001.fastq.gz" $main_dir/split/$fbname"_I2_001.fastq.gz_"

#these 2 lines make the script wait until all scripts have run on the cluster. Not necessary when running scripts locally.
cluster_act=$(qstat -u $USER);
until [[ "$cluster_act" == "" ]]; do echo "Waiting for cluster job to finish.."; sleep 50; cluster_act=$(qstat -u $USER); done


echo "Splitting done, starting to zip them...";
mv $main_dir/split.* $main_dir/split

#Rename the split files
ls -ltrh $main_dir/split | grep $fbname"_I2_001.fastq.gz_" | wc -l | awk '{print $1-1}' > $main_dir/number;
var=$(cat $main_dir/number);

for i in `seq 0 $var`; do pi=`printf %02d $i`; mv $main_dir/split/$fbname"_R1_001.fastq.gz_"${pi} $main_dir/split/"Split_"$pi"_"$fbname"_R1_001.fastq"; done;
for i in `seq 0 $var`; do pi=`printf %02d $i`; mv $main_dir/split/$fbname"_R2_001.fastq.gz_"${pi} $main_dir/split/"Split_"$pi"_"$fbname"_R2_001.fastq"; done;
for i in `seq 0 $var`; do pi=`printf %02d $i`; mv $main_dir/split/$fbname"_I1_001.fastq.gz_"${pi} $main_dir/split/"Split_"$pi"_"$fbname"_I1_001.fastq"; done;
for i in `seq 0 $var`; do pi=`printf %02d $i`; mv $main_dir/split/$fbname"_I2_001.fastq.gz_"${pi} $main_dir/split/"Split_"$pi"_"$fbname"_I2_001.fastq"; done;

#Zip them
echo "Generated "$var+1" files per fastq";
echo "Splitting done, starting to zip the split files...";

cp $main_dir/Scripts/bgzip.sge.sh $main_dir;
for i in `seq 0 $var`; do pi=`printf %02d $i`;qsub $main_dir/bgzip.sge.sh $main_dir/split/"Split_"$pi"_"$fbname"_R1_001.fastq";done;
for i in `seq 0 $var`; do pi=`printf %02d $i`;qsub $main_dir/bgzip.sge.sh $main_dir/split/"Split_"$pi"_"$fbname"_R2_001.fastq";done;
for i in `seq 0 $var`; do pi=`printf %02d $i`;qsub $main_dir/bgzip.sge.sh $main_dir/split/"Split_"$pi"_"$fbname"_I1_001.fastq";done;
for i in `seq 0 $var`; do pi=`printf %02d $i`;qsub $main_dir/bgzip.sge.sh $main_dir/split/"Split_"$pi"_"$fbname"_I2_001.fastq";done;

#Loop until all cluster jobs are finished
cluster_act=$(qstat -u $USER);
until [[ "$cluster_act" == "" ]]; do echo "Waiting for cluster job to finish.."; sleep 50; cluster_act=$(qstat -u $USER); done;
echo "Zipping done, starting to place Barcodes...";
mv $main_dir/bgzip* $main_dir/split

#Place BCs with Script
cp $main_dir/Scripts/SHARE_Cluster_demult.sh $main_dir;
mkdir $main_dir/Demux;

#Caution: One of the files (BC_D.txt) does not exist in /fml/chones/projects/PC060_scCREs/Scripts. Needs to be done beforehand!

for i in `seq 0 $var`; do pi=`printf %02d $i`; qsub $main_dir/SHARE_Cluster_demult.sh $main_dir/split/Split_${pi}"_"$fbname"_" $main_dir/Demux/"Demux_Split_"${pi}"_"$fbname $main_dir;done;

cluster_act=$(qstat -u $USER);
until [[ "$cluster_act" == "" ]]; do echo "Waiting for cluster job to finish.."; sleep 50; cluster_act=$(qstat -u $USER); done

mv $main_dir/SHARE* $main_dir/Demux/
mv $main_dir/BC_* $main_dir/Demux/
echo "Placing Barcodes done, starting to cutadapt...";

#Run cutadapt on resulting files
cp $main_dir/Scripts/cutadapt.sh $main_dir;
mkdir $main_dir/cutadapt;

for i in `seq 0 $var`; do pi=`printf %02d $i`; qsub $main_dir/cutadapt.sh $main_dir/Demux/"Demux_Split_"${pi}"_"$fbname"_R1_001.fastq.gz" $main_dir/Demux/"Demux_Split_"${pi}"_"$fbname"_R2_001.fastq.gz";done;
#output files will have .cutadapt at the end

cluster_act=$(qstat -u $USER);
until [[ "$cluster_act" == "" ]]; do echo "Waiting for cluster job to finish.."; sleep 50; cluster_act=$(qstat -u $USER); done

mv $main_dir/Demux/*.cutadapt $main_dir/cutadapt;
mv $main_dir/cutadapt* $main_dir/cutadapt;
mv $main_dir/Demux/*tooshort* $main_dir/cutadapt;
echo "Done with cutadapt... starting to split by i5 indices";

#Split by i5 index
cp $main_dir/Scripts/Index1Split.sh $main_dir;
cp $main_dir/Scripts/Index2Split.sh $main_dir;


for i in `seq 0 65`; do pi=`printf %02d $i`; qsub ./Index1Split.sh /fml/obi/data/PC060_scCREs/Raw_Reads/Pool_1 $pi $fbname; done
for i in `seq 0 65`; do pi=`printf %02d $i`; qsub ./Index2Split.sh /fml/obi/data/PC060_scCREs/Raw_Reads/Pool_1 $pi $fbname; done

cluster_act=$(qstat -u $USER);
until [[ "$cluster_act" == "" ]]; do echo "Waiting for cluster job to finish.."; sleep 50; cluster_act=$(qstat -u $USER); done

echo "Done splitting by indices... Combining all fastq's with the same index...";

#resulting files are named Demux_Split_${pi}_$fbname_R2_001.fastq.gz.cutadapt_split_D${pu}_001.fastq.gz" (pi = original split number, pu=index number (0-XX)!
#Combine the same indices into one file (so maximum of 96 files)

mkdir $main_dir/index_split;
#R1
for i in `seq 0 496`; do pu=`printf %03d $i`; for x in `seq 0 $var`; do pi=`printf %02d $x`; cat $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R1_001.fastq.gz.cutadapt_split_D"${pu}"_001.fastq.gz"; done >> $main_dir/index_split/"D"${pu}"_"$fbname"_R1_001.fastq.gz";done;

#R2
for i in `seq 0 496`; do pu=`printf %03d $i`; for x in `seq 0 $var`; do pi=`printf %02d $x`; cat $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R2_001.fastq.gz.cutadapt_split_D"${pu}"_001.fastq.gz"; done >> $main_dir/index_split/"D"${pu}"_"$fbname"_R2_001.fastq.gz";done;

echo "All files are split according to i5 index...Deleting intermediate files"
mv $main_dir/Index* $main_dir/cutadapt

#Need to delete all files in between
#also possible to delete in between

#split files from before placing BCs " $main_dir/split/"Split_"$pi"_"$fbname"_R1_001.fastq"
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/split/"Split_"$pi"_"$fbname"_R1_001.fastq.gz"; done;
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/split/"Split_"$pi"_"$fbname"_R2_001.fastq.gz"; done;
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/split/"Split_"$pi"_"$fbname"_I1_001.fastq.gz"; done;
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/split/"Split_"$pi"_"$fbname"_I2_001.fastq.gz"; done;

#cutadapted files AFTER placing BCs
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R1_001.fastq.gz.cutadapt"; done;
#for i in `seq 0 $var`; do pi=`printf %02d $i`; rm $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R2_001.fastq.gz.cutadapt"; done;

#all single i5-split files (main part)
#for i in `seq 0 496`; do pu=`printf %03d $i`; for x in `seq 0 $var`; do pi=`printf %02d $x`; rm $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R1_001.fastq.gz.cutadapt_split_D"${pu}"_001.fastq.gz"; done;done;
#for i in `seq 0 496`; do pu=`printf %03d $i`; for x in `seq 0 $var`; do pi=`printf %02d $x`; rm $main_dir/cutadapt/"Demux_Split_"${pi}"_"$fbname"_R2_001.fastq.gz.cutadapt_split_D"${pu}"_001.fastq.gz"; done;done;

rm $main_dir/number
#leftover files: original FQs, FQs with placed BCs (no cutadapt) and FQs ready to map, one for each index!
