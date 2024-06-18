#!/bin/bash

#$ -pe parallel 20
#$ -l h_vmem=45G
#$ -N STAR.Placement
#$ -o ./STAR_Cluster_out.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -V

if [ ! -e "/tmp/vsoltys" ]
                         then mkdir /tmp/vsoltys
fi

file=$1
echo $file

/fml/chones/local/bin/STAR --outFileNamePrefix ./$file"_" --genomeDir /fml/chones/genome/gbdb/mm10/STAR/ --readFilesIn $file --readFilesCommand zcat --chimOutType WithinBAM --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within
/fml/chones/local/bin/STAR --outFileNamePrefix ./$file"_Aligned.sortedByCoord.out.bam_" --genomeDir /fml/chones/genome/gbdb/mm10/STAR/STAR_consensus_GTF --readFilesIn $file"_Aligned.sortedByCoord.out.bam" --readFilesCommand samtools view --chimOutType WithinBAM --limitSjdbInsertNsj=2000000 --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate --readFilesType SAM SE --sjdbFileChrStartEnd ./$file"_SJ.out.tab" --outSAMunmapped Within
