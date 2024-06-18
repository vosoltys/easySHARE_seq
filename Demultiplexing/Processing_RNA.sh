#Processing demultiplexed easySHARE-seq RNA-seq data
#replace INDEX with your sample (i5 index) number, e.g. 001
#all scripts in /Scripts

#extract UMIs from Read 1 and add to the sequence name
for i in INDEX; do pi=`printf %02d $i`; qsub ./umi_extract.sh ./D0${pi}_Pool-X_R1_001.fastq.gz ./D0${pi}_Pool-X_R2_001.fastq.gz; done

#only keep reads with TTTTT at beginning of polyT tail (> 96%). Also appends cell barcode to readname
for i in INDEX; do pi=`printf %03d $i`; qsub ./Calc_Valid_Reads.sh $pi D"${pi}"_Pool-X_R1_001.fastq.gz.UMI.fastq.gz D"${pi}"_Pool-X_R2_001.fastq.gz.UMI.fastq.gz 1; done

#Mapping with STAR (twopass). This also outputs unmapped reads!
for i in *_R1_UMI_001.fastq.gz; do qsub ./STAR.sh $i; done

#Moving Barcode and UMI from sequencing name into easy to work with bam flags
for i in INDEX_LANE; do samtools view -H ./"$i"_R1_UMI_001.fastq.gz_Aligned.sortedByCoord.out.bam_Aligned.sortedByCoord.out.bam > "$i".prelim.header; 
samtools view ./"$i"_R1_UMI_001.fastq.gz_Aligned.sortedByCoord.out.bam_Aligned.sortedByCoord.out.bam | awk '/^@/ {print;next} {N=split($1,n,"_");print $0 "\tCR:Z:" n[N]}' | awk '/^@/ {print;next} {print substr($1,1,length($1)-13) "\t" $0}' | cut -f1,3- | awk '/^@/ {print;next} {N=split($1,n,"_");print $0 "\tBX:Z:" n[N]}' | awk '/^@/ {print;next} {print substr($1,1,length($1)-11) "\t" $0}' | cut -f1,3- > "$i".prelim.sam; cat "$i".prelim.header "$i".prelim.sam | samtools view -@ 10 -bh - > "$i".bam; rm "$i".prelim.header "$i".prelim.sam; done

#sort and index 
for i in INDEX_LANE; ; do samtools sort -@ 25 $i.bam > $i.sorted.bam; samtools index -@ 25 $i.sorted.bam; done

#Counting overlapping transcripts using featureCounts. 
for i in *.sorted.bam; do featureCounts -a ./GTFs/Mus_musculus.GRCm38.102.ConsensusFiltered.gtf.gz -o "$i".exononly.featureCounts -g gene_name -t gene -R BAM -T 12 -M -s 1 --primary $i; done 

#filter
for i in *featureCounts.bam; do samtools view -@ 20 -h -F 256 -q 5 -d HI:1 $i | grep -v A000 | grep -v C000 | samtools view -@ 25 -bh - | samtools sort -@ 30 - > "$i".filtered.sorted.bam; samtools index -@ 20 "$i".filtered.sorted.bam; done

#generate count matrix
for i in *filtered.sorted.bam; do num=$(echo $i | cut -c1-5); qsub ./Count_Matrix.sh $i $num; done

#Done
