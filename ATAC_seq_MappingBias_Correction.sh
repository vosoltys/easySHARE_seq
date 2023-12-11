
## General workflow fro calculating per-peak mapping bias correction factor

#download Sanger PWK genome and gDNA sequencing data
#Importantly, scATAC-seq data and this gDNA data should the experience the same filtering requirements
#call peaks


#First, get fasta for peaks
bedtools getfasta -fi mm10.fa -bed Liver_peaks.bed -fo Liver_peaks.fa
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%400==0){file=sprintf("myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < Liver_peaks.fa

#BLAT the mm10 fastas for the peaks onto the PWK genome 
blat -noHead -minScore=88 PWK.fa $file $file.out.psl
for i in *.out.psl; do cat $i | awk '{if ($1/$11 > 0.90) print $0}' | pslFilter stdin stdout -noHead >> All.BLAT.psl; done
#filter out aligned segments lower than 90% aligned and segments with more than 80 mismatches
cat All.BLAT.psl | awk '{if (($1+$2)/$1 >= 0.90) print $0}' | awk '{if ($2 <= 80) print $0}' > All.filter.BLAT.psl 
#if multiple hits, extract the one with the longest aligned length. T~he filter out hits across chromosomes
cat Liver_peaks.bed | awk '{print $1":"$2"-"$3}' > Liver_peaks.names
for i in `seq 1 $num_peaks`; do peak=$(head -"$i" Liver_peaks.names | tail -1); cat All.filter.BLAT.psl | grep $peak | sort -k1 -n -r | head -1 >> All.filter.BLAT.final.psl; done
All.filter.BLAT.final.psl | sed 's/:/\t/g' | awk 'BEGIN {OFS="\t"}{if ($10 == $15) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10":"$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22}' > temp; mv temp All.filter.BLAT.final.psl

#extract the coordinates
cat All.filter.BLAT.final.psl | awk '{print $10"\t"$11"\t"$14"\t"$16"\t"$17"\t"$17-$16}' | sed 's/-/\t/g' | sed 's/:/\t/g' | bedtools sort -i - > Liver_peaks.PWD
for i in `seq 1 $num_hits`; do echo $i >> temp; done; paste Liver_peaks.PWD temp > temp2; mv temp2 Liver_peaks.PWD
rm temp 

#Make sure new coordinates are not overlapinng
cat Liver_peaks.PWD | cut -f5,6,7,9 | bedtools sort -i - > all
cat Liver_peaks.PWD | cut -f5,6,7 | bedtools sort -i - > all.test
bedtools intersect -a all -b all.test -c | awk '{if ($5 == "1") print $0}' > 1s
bedtools intersect -a all -b all.test -c | awk '{if ($5 == "2") print $0}' | sed -n 1~2p | cut -f1-4 | cat 1s - | cut -f1-4 | bedtools sort -i - > temp
cat temp | cut -f4 | sort > ids
cat Liver_peaks.PWD | sort -k9 | join -1 9 -2 1 - ids | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$1}' | bedtools sort -i - > temp2
echo -e BL6_chr"\t"BL6_Start"\t"BL6_End"\t"BL6_length"\t"PWD_chr"\t"PWD_Start"\t"PWD_end"\t"PWD_length"\t"ID | cat - temp2 > Liver_peaks.final


#make files for BL6 and PWD only
tail -n+2 Liver_peaks.final | cut -f1,2,3,4,9 > Liver_peaks.final.BL6 
tail -n+2 Liver_peaks.final | cut -f5,6,7,8,9 | bedtools sort -i - > Liver_peaks.final.PWD

#Now, we extract the exact same amount of reads per matched peak in both chromosomes
for i in Liver_peaks.final.BL6; do cat $i | awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $5, $1, $2+1, $3,"."}' > Liver_peaks.final.BL6.saf; done
for i in Liver_peaks.final.PWD; do cat $i | awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $5, $1, $2+1, $3, "."}' > Liver_peaks.final.PWD.saf; done

#counting with the processed Sanger gDNA files
featureCounts --fracOverlap 0.9 -T 35 --primary -a Liver_peaks.final.BL6.saf -F SAF -o BL6.AllCounts C57BL_6NJ.filtered.bam
featureCounts --fracOverlap 0.9 -T 25 --primary -a Liver_peaks.final.PWD.saf -F SAF -o PWD.AllCounts PWK_PhJ.filtered.bam

#How maany reads can we use?
tail -n+3 .AllCounts | sort -k1 -n > temp; tail -n+3 BL6.AllCounts | paste - temp | awk '{print $2,$3,$4,$6,$9,$10,$11,$13,$1,$7,$14}' | sed 's/\s/\t/g' | awk 'BEGIN{FS=OFS="\t"; print "Chr_BL6\tStart_BL6\tEnd_BL6\tLength_BL6\tChr_PWD\tStart_PWD\tEnd_PWD\tLength_PWD\tWindow_ID\tFragments_BL6\tFragments_PWD"}{print $0}' > BL6_PWD.Fragment_Counts_All
tail -n+2 BL6_PWD.Fragment_Counts_All | awk '{if (($10 <8) || ($11<8)) skip; else print $0}' | awk '{if ($10 > $11) print $0"\t"$11; else print $0"\t"$10}' | awk 'BEGIN{FS=OFS="\t"; print "Chr_BL6\tStart_BL6\tEnd_BL6\tLength_BL6\tChr_PWD\tStart_PWD\tEnd_PWD\tLength_PWD\tWindow_ID\tFragments_BL6\tFragments_PWD\tFrag_to_extract"}{print $0}' > temp; mv temp BL6_PWD.Fragment_Counts_All

#Extract
for i in `seq 2 65133`; do 
region=$(awk -v var="$i" 'NR==var' BL6_PWK.Fragment_Counts_All | awk '{print $1":"$2"-"$3}')
number=$(awk -v var="$i" 'NR==var' BL6_PWK.Fragment_Counts_All | awk '{print $12}'); 
samtools view -@ 20 C57BL_6NJ.filtered.bam "$region" | shuf -n $number >> BL6.R1_overlap_windows.sam; done
for i in `seq 2 65133`; do
region=$(awk -v var="$i" 'NR==var' BL6_PWK.Fragment_Counts_All | awk '{print $5":"$6"-"$7}')
number=$(awk -v var="$i" 'NR==var' BL6_PWK.Fragment_Counts_All | awk '{print $12}')
samtools view -@ 20 PWK_PhJ.filtered.bam "$region" | shuf -n $number >> PWK.R1_overlap_windows.sam; done

#Convert to fastq and map to mm10
bedtools bamtofastq -i BL6.R1_overlap_windows.bam -fq ./BL6.R1.fastq; pigz -p 20 BL6_R1_001.fastq
bedtools bamtofastq -i PWK.R1_overlap_windows.bam -fq ./PWK.R1.fastq; pigz -p 20 PWK_R1_001.fastq

bwa mem -t 10 /mm10.fa $file -R "@RG\tID:$fbname\tSM:$fbname\tLB:SHARE_AB_BL6\tPL:Illumina.NovaSeq.1x150" | samtools view -bh - > ./$fbname.mm10.bam
samtools sort -@ 20 -l 9 -T ./$fbname.tmpsort -o ./$fbname.sorted.mm10.bam ./$fbname.mm10.bam; samtools index ./$fbname.sorted.mm10.bam
samtools view -@ 20 -bh -F 4 -q 30 $file.sorted.mm10.bam > $file.filtered.bam; samtools index -@ 20 $file.filtered.bam

#Finally, calculate how much mapping bias per peak there is
featureCounts --fracOverlap 0.9 -T 35 --primary -a Liver_peaks.final.BL6.saf -F SAF -o BL6.mm10.AllCounts BL6.filtered.bam
featureCounts --fracOverlap 0.9 -T 35 --primary -a Liver_peaks.final.BL6.saf -F SAF -o PWD.mm10.AllCounts PWD.filtered.bam
paste BL6.mm10.AllCounts PWD.mm10.AllCounts | tail -n+3 |awk '{print $2"\t"$3"\t"$4"\t"$6"\t"$9"\t"$10"\t"$11"\t"$13"\t"$1"\t"$7"\t"$14}' | awk '{if (($10==0) || ($11==0)) skip; else print $0}' | awk '{if (($10==0) && ($11==0)) skip; else print $0}' | awk '{print $0"\t"$10/($10+$11)}' > BL6_PWD.mm10.Bias
cat BL6_PWD.mm10.Bias | awk '{if (($10+$11) < 12) skip; else print $1"\t"$2"\t"$3"\t"$12}' > BL6_PWD.mm10.peakbias
#add peaks that have been filtered out back in and set them to 0.5
bedtools intersect -v -a Liver_PeakSet_peaks.narrowPeak -b BL6_PWD.mm10.peakbias | awk '{print $1"\t"$2"\t"$3"\t0.5"}' | cat - BL6_PWD.mm10.peakbias | bedtools sort -i - > BL6vPWD_PeakBias.final
cat BL6vPWD_PeakBias.final | awk '{print $0"\t"$4/0.5}' | awk '{print $0"\t"2-$5}' > temp; echo -e chr"\t"start"\t"end"\t"Bias_BL6_to_PWD"\t"PWD_multiplier"\t"BL6_multiplier | cat - temp > BL6vPWD_PeakBias.final




