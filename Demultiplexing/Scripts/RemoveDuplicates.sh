#!/bin/bash

#$ -pe parallel 20
#$ -l h_vmem=16G
#$ -N RemoveDuplicates
#$ -o ./RemDup.txt
#$ -j y
#$ -S /bin/bash
#$ -cwd


file=$1;
echo $file

if [ ! -e "/tmp/vsoltys" ]
                         then mkdir /tmp/vsoltys
fi


/usr/bin/java -Xmx4g -jar /fml/chones/local/picard-2.18.25/picard.jar MarkDuplicates INPUT=$file OUTPUT=$file.remDupl.bam METRICS_FILE=$file.dup_metrics READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=TRUE

