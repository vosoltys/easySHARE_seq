#Processing demultiplexed easySHARE-seq ATAC-seq data
#replace INDEX with your sample (i5 index) number, e.g. 001
#all scripts in /Scripts

#mapping to mm10
for i in /fml/obi/data/PC060_scCREs/Raw_Reads/ATACseq_fastqs/*R1*.gz; do qsub ./mm10.bwa.placement.sh $i; done

#Removing duplicates
for i in *.sorted.mm10.bam; do qsub ./RemoveDuplicates.sh $i; done

#filtering
for i in *.sorted.mm10.bam.remDupl.bam; do qsub ./Filter_Bam.sh $i; done

#Calling Peaks
#for calling peaks, subsample all available bamfiles to a simmilar depth, then use MACS2
for i in INDEX_LANE ; do x=$(echo "$i" | sed 's/_/_Pool-/g'); samtools flagstat -@ 30./D"$x".sorted.mm10.bam.remDupl.bam.filtered.bam > "$i".flagstat; done
for i in *.flagstat; do head -1 $i | awk -v x="$i" '{print $1"\t"x}' >> temp; done
echo -e Reads"\t"Sample | cat - temp | sed 's/.flagstat//g' > Read.Count; rm temp
cat /fml/obi/data/PC060_scCREs/ATACseq_Analysis/ATAC_Samples.txt  | sort -k2 > temp.all
tail -n+2 Read.Count | sort -k2 | join -1 2 -2 2 - temp.all | awk '{print $1"\t"$2"\t"$3}' | sed 's/_/\t/g' | awk '{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > temp; echo -e ID"\t"Reads"\t"Genotype"\t"Sex"\t"Tissue | cat - temp > Read.Count 
rm temp temp.all
tail -n+2 Read.Count | sed 's/_/\t/g' | awk '{print $1"_"$2"\t"$3"\t"$4"_Pool-"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > temp; head -1 Read.Count > temp2; cat temp2 temp > Read.Count ; rm temp temp

#change number in awk command to smallest available read number
tail -n+2 Read.Count | sort -k2 -n | awk '{print $0"\t"19325129/$2}' > temp; head -1 Read.Count > temp2; cat temp2 temp | sed 's/_/_Pool-/g' | sed 's/0\.//g' > Read.Count ; rm temp temp2

for i in `seq 1 8`; do file=$(tail -n+2 Read.Count | awk -v x="$i" 'NR==x' | awk '{print $1}'); int=$(tail -n+2 Read.Count | awk -v x="$i" 'NR==x' | awk '{print $5}'); qsub ./Subsample_Bam.sh /fml/obi/data/PC060_scCREs/ATACseq_Analysis/bamfiles/final_bamfiles/D"$file".sorted.mm10.bam.remDupl.bam.filtered.bam.final.bam $int; done
ls -ltrh | grep bam$ | awk '{print $9}' > merge.list
samtools merge -@ 25 -b merge.list B6_PWD.pooled.bam
rm merge.list *subs.bam*

/fml/chones/local/bin/macs2 callpeak -t B6_PWD.pooled.bam -g mm -n B6_PWD_Liver.peaks -B --nomodel -f BAMPE --extsize 150 --nomodel --min-length 100 --cutoff-analysis --keep-dup all --trackline


#making fragment files (not on subsampled bams)
for i in *bam; do echo $i; sinto fragments -b $i -f $i.fragments.txt -p 3 -t "BX"; done
for i in *txt; do echo $i | cut -c1-11 | sed 's/Pool-//g' > temp; name=$(cat temp); mv $i $name.txt; done
for i in *txt; do sort -k 1,1 -k2,2n $i > $i.frags.sort.bed; bgzip $i.frags.sort.bed; tabix -p bed $i.frags.sort.bed.gz; done

#getting the top 5000 barcodes as not to generate a gigantic matrix
for i in ./*filtered.bam; do bxtools stats $i > $i.stats; done
for i in *stats; do echo $i | cut -c1-11 | sed 's/Pool-//g' > temp; name=$(cat temp); mv $i $name.stats; rm temp; done
for i in *.stats; do cat $i | sort -k2,2 -n -r | grep -v A000 | grep -v C000 | head -5000 | cut -f1 > $i.Cells; done
