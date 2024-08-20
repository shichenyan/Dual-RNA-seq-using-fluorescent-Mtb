#!/bin/bash

wdir="/home/scratch/dualseq"
readsdir="/home/Data/Dualseq"

cd $wdir
mkdir QC
mkdir trimmed
mkdir align
mkdir count

######## fastqc #########
cd $wdir/QC
fq1=$readsdir/${1}_1.fq.gz
fq2=$readsdir/${1}_2.fq.gz
fastqc -t 10 $fq1 -o ./
fastqc -t 10 $fq2 -o ./

######## trimming ###########
cd $wdir/trimmed
trim_galore -q 25 --phred33 --illumina --length 35 --stringency 2 --paired -o ./ $fq1 $fq2
ls ./${1}*.gz | xargs fastqc -t 10  -o ./

######## Bowtie2 – to separate rRNA ###########
### build rRNA index ###
grep ‘rRNA’ Homo_sapiens.GRCh38.99.gtf  
gffread hs.GRCh38.rrna.gtf -g Homo_sapiens.GRCh38.dna.primary_assembly.fa -w hs.GRCh38.rrna.fa
grep 'gene_biotype "rRNA"' GCF_000016145.1_ASM1614v1_genomic.gtf > GCF_000016145.1_Ra.rrna.gtf
gffread GCF_000016145.1_Ra.rrna.gtf -g GCF_000016145.1_ASM1614v1_genomic.fna -w GCF_000016145.1_Ra.rrna.fa
cat GCF_000016145.1_Ra.rrna.fa Mycobacterium_tuberculosis_2021-05-05_silva.fasta hs.GRCh38.rrna.fa > Mtb_Hs_rrna.fasta
bowtie2-build -f Mtb_Hs_rrna.fasta Mtb_Hs_rrna

### mapping to remove rRNA ###
cd $wdir/trimmed
cat ../trimmed/${1}_1_val_1.fq.gz ../trimmed/${1}_2_val_2.fq.gz > ${1}.fq.gz
pigz -d -p 20 -c ${1}.fq.gz > ${1}.fq
bowtie2 -p 20 -x /home/DB/rRNA_SILVA/Mtb_Hs_rrna -U ${1}.fq \
--sensitive --al ${1}_aligned_rrna.fastq --un ${1}_unaligned_rrna.fastq \
-S ${1}.rrna.sam
rm ${1}.fq

samtools view -Su ${1}.rrna.sam | samtools sort - -o ${1}.rrna.sort.bam
samtools index ${1}.rrna.sort.bam
samtools flagstat ${1}.rrna.sort.bam > $(basename ${1}.rrna ".sort.bam").flagstat
rm ${1}.rrna.sam

######## Bowtie2 – align TB reads and split non-matching ones ###########
cd $wdir/align
bowtie2 -x /home/DB/index/genome_gtf/TB/H37Ra/Mtb -U ${1}_unaligned_rrna.fastq \
--very-sensitive --al ${1}_Mtb.fastq --un ${1}_hs.fastq -S ${1}_Mtb.bowtie2.sam
	

######## Hisat2 – align Mtb reads ###########
cd $wdir/align
hisat2 -p 20 --rg-id=TB_${1} --rg SM:TB_${1} --rg LB:TB_${1} --dta \
-x /home/DB/index/genome_gtf/TB/H37Ra/GCF_000016145.1_hisat2_index \
--rna-strandness R -U ${1}_Mtb.fastq -S ${1}_Mtb.hisat2.sam
samtools sort -@ 20 -o ${1}_Mtb.hisat2.bam ${1}_Mtb.hisat2.sam
find ${1}_Mtb.hisat2.bam -exec echo samtools index {} \; | sh
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse \
--minaqual 1 --type gene --idattr gene_id ${1}_Mtb.hisat2.bam \
/home/DB/index/genome_gtf/TB/H37Ra/GCF_000016145.1_ASM1614v1_genomic.gtf > ../count/${1}_Mtb.tsv

samtools flagstat ${1}_Mtb.hisat2.bam > $(basename ${1}_Mtb ".hisat2.bam").flagstat	
	
	
######## Hisat2 – align Hs reads ###########
cd $wdir/align
hisat2 -p 20 --rg-id=hs_${1} --rg SM:hs_${1} --rg LB:hs_${1} --dta \
-x /home/DB/genome_gtf/Hs_GRCh38/Homo_sapiens.GRCh38 \
--rna-strandness R -U ${1}_hs.fastq -S ${1}_hs.hisat2.sam
samtools sort -@ 20 -o ${1}_hs.hisat2.bam ${1}_hs.hisat2.sam
find ${1}_hs.hisat2.bam -exec echo samtools index {} \; | sh
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse \
--minaqual 1 --type exon --idattr gene_id ${1}_hs.hisat2.bam \
/home/DB/index/genome_gtf/Hs_GRCh38/Homo_sapiens.GRCh38.99.gtf > ../count/${1}_hs.tsv

samtools flagstat ${1}_hs.hisat2.bam > $(basename ${1}_hs ".hisat2.bam").flagstat
