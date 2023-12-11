#!/bin/bash

## Describes the process of artificial genome mapping starting with demultiplexed files

#Build artificial genome, starting with VCF filtered for indicative positions
#Positions need to be homozygous
zcat PWK.ALTHET.SNP.vcf.gz | grep -v "#" | awk '{$7="PASS"; print $0}' | sed 's/\s/\t/g' > tmp.PWD.SNPs
cat PWK.Indels.vcf | grep -v "#" | awk '{$7="PASS"; print $0}' | sed 's/\s/\t/g' > tmp.PWD.InDels
zcat PWK.ALTHET.SNP.vcf.gz | grep "#" > tmp.SNPs.head
cat PWK.Indels.vcf | grep "#" > tmp.InDels.head
cat tmp.SNPs.head tmp.PWD.SNPs | gzip -c - > SNPs.AllHom.PWD.vcf.gz
zcat SNPs.AllHom.PWD.vcf.gz | awk 'BEGIN{FS=OFS="\t"} {gsub(/^[0-3].[0-3]/,"1/1",$10)} 1' | gzip -c - > temp2; mv temp2 SNPs.AllHom.PWD.vcf.gz 
cat tmp.InDels.head tmp.PWD.InDels | gzip -c - > InDels.AllHom.PWD.vcf.gz
zcat InDels.AllHom.PWD.vcf.gz | awk 'BEGIN{FS=OFS="\t"} {gsub(/^[0-3].[0-3]/,"1/1",$10)} 1' | gzip -c - > temp; mv temp InDels.AllHom.PWD.vcf.gz 

#Construct genome
java -Xmx40000m -jar ../vcf2diploid_v0.2.6a/vcf2diploid.jar -id PWK_PhJ -chr mm10.fa -vcf SNPs.AllHom.PWD.vcf.gz InDels.AllHom.PWD.vcf.gz -outDir ./ > pwd.log.txt
for i in chr*paternal.fa; do echo $i; cat $i >> PWD_Artif_Genome.fa;done
mv paternal.chain mm10_TO_PWD_Art.chain

#Lift over GTF + VCFs
#
liftOver -gff mm10.gtf mm10_TO_PWD_Art.chain PWD_Art.gtf unmapped.gtf
java -Xmx16g -jar /fml/chones/local/picard-2.18.25/picard.jar CreateSequenceDictionary R=PWD_Artif_Genome.fa O=PWD_Artif_Genome.dict
java -Xmx32g -jar picard.jar LiftoverVcf I=PWK.ALTHET.SNP.vcf.gz  O=PWK.SNP.liftover.vcf.gz  CHAIN=mm10_TO_PWD_Art.chain REJECT=out.vcf.gz R=PWD_Artif_Genome.fa RECOVER_SWAPPED_REF_ALT=true MAX_RECORDS_IN_RAM=400000
java -Xmx32g -jar picard.jar LiftoverVcf I=PWK.Indels.vcf O=PWK.InDel.liftover.vcf.gz  CHAIN=mm10_TO_PWD_Art.chain REJECT=out.InDel.vcf.gz R=PWD_Artif_Genome.fa RECOVER_SWAPPED_REF_ALT=true MAX_RECORDS_IN_RAM=400000



## Process scRNA-seq fastqs
#extract UMIs, add to BC name
umi_tools extract --bc-pattern=NNNNNNNNNN --stdin $R2 --stdout $R2.UMI.fastq.gz --read2-in $R1 --read2-out $R1.UMI.fastq.gz -L extract.log

#Only retain reads with polyT-tail
pi=$1
R1=$2
R2=$3
Pool=$4

zcat $R2 | bioawk -c fastx seq | awk '{if (($2 ~ /^TTT[AGCT]T/) || ($2 ~ /^[AGCT]TTTT/) || ($2 ~ /^T[AGCT]TTT/) || ($2 ~ /^TT[AGCT]TT/) || ($2 ~ /^TTTT[AGCT]/) || ($2 ~ /^[AGCT]TTTT/)) print $0}' | awk '{print $1}' > "$Pool"_"$pi"_Valid_Reads; cat "$Pool"_"$pi"_Valid_Reads | split -a 2 -d -l 5000000 - "$Pool"_Valid_Reads_"$pi"_
for i in `seq 0 30`; do x=`printf %02d $i`; zcat $R1 | grep -a -A3 -F -f "$Pool"_Valid_Reads_"$pi"_"$x" | awk '{if ($0 != "--") print $0}' | gzip -c - >> "$R1".filtered.fastq.gz; done

#As a last step, move the barcode of Read1 into the readname
zcat "$R1".filtered.fastq.gz | /fml/chones/local/bin/bioawk -c fastx seq | awk '{print "@"$1"_"substr($4,RSTART+6,12) "\n" $2 "\n" "+" "\n" $3}' | gzip -c > "$pi"_"$Pool"_R1_UMI_001.fastq.gz

#Firstpass mapping with STAR to both mm10 and AG
STAR --outFileNamePrefix /firstpass_bamfiles_mm10/$file"_" --genomeDir mm10/STAR/ --readFilesIn $file --readFilesCommand zcat --chimOutType WithinBAM --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate
STAR --outFileNamePrefix /firstpass_bamfiles_Art/$file"_" --genomeDir PWD_ArtificialGenome/STAR/ --readFilesIn $file --readFilesCommand zcat --chimOutType WithinBAM --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate

#Twopass mapping with STAR with newly found SJs
STAR --outFileNamePrefix twopass_bamfiles_mm10/$file"_" --genomeDir mm10/STAR --readFilesIn /firstpass_bamfiles_mm10/$file --readFilesCommand samtools view --chimOutType WithinBAM --limitSjdbInsertNsj=2000000 --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate --readFilesType SAM SE --sjdbFileChrStartEnd /firstpass_bamfiles_mm10/$file_SJ.out.tab
STAR --outFileNamePrefix twopass_bamfiles_Art/$file"_" --genomeDir /PWD_ArtificialGenome/STAR/ --readFilesIn /firstpass_bamfiles_Art/$file --readFilesCommand samtools view --chimOutType WithinBAM --limitSjdbInsertNsj=2000000 --limitOutSJcollapsed 5000000 --outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --outSAMtype BAM SortedByCoordinate --readFilesType SAM SE --sjdbFileChrStartEnd $file_SJ.out.tab

#For each twopass output, move UMI and cell barcode from read name to bam tag
samtools view -H $bam > "$i".prelim.header; samtools view $bam | awk '/^@/ {print;next} {N=split($1,n,"_");print $0 "\tCR:Z:" n[N]}' | awk '/^@/ {print;next} {print substr($1,1,length($1)-13) "\t" $0}' | cut -f1,3- | awk '/^@/ {print;next} {N=split($1,n,"_");print $0 "\tBX:Z:" n[N]}' | awk '/^@/ {print;next} {print substr($1,1,length($1)-11) "\t" $0}' | cut -f1,3- > "$i".prelim.sam; cat "$i".prelim.header "$i".prelim.sam | samtools view -@ 10 -bh - > "$i".bam; rm "$i".prelim.header "$i".prelim.sam; done

#After sorting, count reads overlapping cDNA and add tag for mm10 and AG
featureCounts -a mm10.gtf -o "$i".exononly.featureCounts -g gene_name -t gene -R BAM -T 12 -M -s 1 --primary $i.mm10.sorted.bam; done
featureCounts -a PWD_ArtificialGenome/PWX_Art.gtf -o "$i".exononly.featureCounts -g gene_name -t gene -R BAM -T 12 -M -s 1 --primary $i.art.sorted.bam; done 

#Compare alignment quality across both alignments and select better mapping read
i=$1
samtools view -@ 30 -F 260 -d HI:1 "$i".mm10.sorted.bam.featureCounts.bam | cut -f1,5 | sort -k1 --parallel=15 -S 50% > "$i".mm10.mapQ; 
samtools view -@ 30 -F 260 -d HI:1 "$i".art.sorted.bam.featureCounts.bam | cut -f1,5 | sort -k1 --parallel=15 -S 50% > "$i".pwd.mapQ;

#get reads that are only in one or the other!
cat "$i".mm10.mapQ | sed 's/\s/\t/g' | awk '{print $1}' | uniq > "$i".mm10.names;
cat "$i".pwd.mapQ | sed 's/\s/\t/g' | awk '{print $1}' | uniq > "$i".pwd.names;
grep -v -F -f "$i".pwd.names "$i".mm10.names  >> "$i".mm10.reads;
grep -v -F -f "$i".mm10.names "$i".pwd.names >> "$i".pwd.reads;

#what to do with both?
join -j1 "$i".mm10.mapQ "$i".pwd.mapQ | sed 's/\s/\t/g' | awk '{if (($2==$3) || ($2 > $3)) print $1}' | uniq >> "$i".mm10.reads;
join -j1 "$i".mm10.mapQ "$i".psc.mapQ | sed 's/\s/\t/g' | awk '{if ($2 < $3) print $1}' | uniq >> "$i".pwd.reads;

#filter the bamfiles
samtools view -@ 30 -bh -N "$i".mm10.reads "$i".mm10.sorted.bam.featureCounts.bam | samtools sort -@ 15 - > "$i".mm10.bam; samtools index -@ 30 "$i".mm10.bam;
samtools view -@ 30 -bh -N "$i".pwd.reads "$i".art.sorted.bam.featureCounts.bam | samtools sort -@ 15 - > "$i".pwd.bam; samtools index -@ 30 "$i".pwd.bam;

#De-duplication for each bam
umi_tools dedup --extract-umi-method=tag -I $file -S $file.dedup.new.bam --umi-tag=BX --per-gene --gene-tag=XT --per-cell --cell-tag=CR --multimapping-detection-method=NH
samtools index -@ 5 $file.dedup.new.bam

#Merge & de-duplicate a second time
samtools view -h -@ 30 "$i".art.dedup.new.bam | sed 's/_paternal//g' | samtools view -bh -@ 20 - > "$i".art.bam; 
samtools index -@ 20 "$i".art.bam; samtools merge -o "$i".combined.bam "$i".mm10.dedup.new.bam "$i".art.dedup.new.bam; samtools sort -@ 20 "$i".combined.bam  > "$i".temp; mv "$i".temp "$i".combined.bam; samtools index -@ 20 "$i".combined.bam; 
umi_tools dedup --extract-umi-method=tag -I "$i".combined.bam -S $file.final.bam --umi-tag=BX --per-gene --gene-tag=XT --per-cell --cell-tag=CR --multimapping-detection-method=NH

#Make single-cell matrices & Sample matrices
umi_tools count --umi-tag=BX --extract-umi-method=tag --per-gene --gene-tag=XT --assigned-status-tag=XS --stdin=$file.finnal.bam --log=$file.count.log --wide-format-cell-counts --stdout=$file.csv
umi_tools count --extract-umi-method=tag --cell-tag=CR --umi-tag=BX --per-gene --gene-tag=XT --per-cell --assigned-status-tag=XS --stdin=$file --log=$num.log  --stdout=$file.cells.csv







