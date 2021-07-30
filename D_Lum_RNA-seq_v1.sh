# Lum_RNA-seq_v1.sh
# Enrico Barrozo, BCM, Aagaard Lab, 2021

# Data from Lum, et al., Clinical & Translational Immunology (2019)	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6831931/pdf/CTI2-8-e01082.pdf
# GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139181

## 33 Illumina sequencing files (PE 2x125 bp reads); see SraRunTable.txt-5 for Metadata
## 5 patients: 2 healthy controls, 3 patients with ZIKV infections at different trimesters including 1, 2, and 3
## 3 samples included: 1) fetal membrane, 2a-b) placental disc separated for (a) CD45+ or (b) CD45- fractions
## 3 technical replicates for all except patient 3 (trimester 3)


# login to aagaardlab3
ssh <>

# make lum directory and tree
mkdir lum
cd lum
mkdir scripts
mkdir data
mkdir docs
mkdir results
cd results
mkdir qc
cd ..

docker pull inutano/sra-toolkit
docker run -it -v $PWD:$PWD --name "Lum-download" inutano/sra-toolkit bash
vdb-config --interactive
# set to default, save, then exit
cd /home/ebarrozo/lum/data
fasterq-dump --help
fasterq-dump --progress --split-files --skip-technical SRR10318591 SRR10318592 SRR10318596 SRR10318597 SRR10318598 SRR10318599 SRR10318600 SRR10318601 SRR10318605 SRR10318606 SRR10318607 SRR10318608 SRR10318609 SRR10318610 SRR10318614 SRR10318615 SRR10318616 SRR10318617 SRR10318618 SRR10318619 SRR10318621 SRR10318622
fasterq-dump --progress --split-files --skip-technical SRR10318590 SRR10318593 SRR10318594 SRR10318595 SRR10318602 SRR10318603 SRR10318604 SRR10318611 SRR10318612 SRR10318613 SRR10318620 
## Once again errors with the second group of data (foetal membranes). Going to leave those out for now. 

cd /home/ebarrozo/lum/data
docker logs Lum-download
docker start -i Lum-download

## if it is running- check logs as they come out
docker logs -f Lum-download
#ctrl + c to exit
docker stop Lum-download
docker rm Lum-download


gunzip SRR10318591_1.fastq
fastqc --threads 32 --outdir /home/ebarrozo/lum/results/qc SRR10318591_1.fastq.gz
fastqc --threads 32 --outdir /home/ebarrozo/lum/results/qc SRR10318590_1.fastq.gz

cd /home/ebarrozo/lum
docker run -it -v $PWD:$PWD --name "lum-docker" continuumio/anaconda3 bash 

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install trimmomatic fastqc star picard samtools
conda install -c bioconda bedtools ucsc-bedgraphtobigwig htseq
## htseq takes a long time
conda install -c biobuilds fastx-toolkit
docker logs -f lum-docker
##		ctrl + c to exit [dont do ctrl+z twice, which exits conda]

# docker stop lum-docker
# docker rm lum-docker

# docker start -i lum-docker

cd /home/ebarrozo/lum/data

# Perform QC on all files and determine length, adapter/index retention.
for f in *.fastq; do
fastqc --threads 32 --outdir /home/ebarrozo/lum/results/qc $f
done

## Are adapters retained?
## What is the read length? Is it 125?

## Trim adapter and barcodes using trimmomatic, installed in a conda environment
#for f in *.fastq.gz;
#do
#trimmomatic SE -threads 32 $f ${f%.fastq.gz}.trimmed.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20
#done

## Trimmomatic PE loop
for file in *_1.fastq; do
withpath="${file}"
filename=${withpath##*/}
base="${filename%*_*.fastq}"
echo "${base}"
trimmomatic PE -threads 32 "${base}"*_1.fastq "${base}"*_2.fastq "${base}"_1.trimmed_PE.fastq "${base}"_1.trimmed_SE.fastq "${base}"_2.trimmed_PE.fastq "${base}"_2.trimmed_SE.fastq SLIDINGWINDOW:4:20 MINLEN:20
done

## Remove the files that won't be used
for f in *SE.fastq; do
rm $f 
done
for f in *2.trimmed_PE.fastq; do
rm $f 
done

## see if adding the adapters found in the first fastqc run fixed the problem
for f in *_1.trimmed_PE.fastq; do
fastqc --threads 32 --outdir /home/ebarrozo/lum/results/qc $f
done

# Unzip fastq files for alignment
#for f in *trimmed.fastq.gz; do
#gunzip $f
#done

## Index STAR genome at appropriate lengths (nt from qc minus 1)
## 	Generate a STAR genome build
cd /home/ebarrozo/lum/docs
mkdir ref
cd ref
mkdir STAR
cd /home/ebarrozo/lum/docs/ref/STAR
cp /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa .
cp /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf .

sed 's/gene_id "gene-/gene_id "/' /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf > /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.2.gtf
sed 's/transcript_name "rna-/transcript_name "/' /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.2.gtf > /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.3.gtf
cd /home/ebarrozo/lum/docs/ref/STAR
cp /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.3.gtf .

## What is the read length after trimming? --sjdbOverhang 124?
STAR \
--runThreadN 32 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/lum/docs/ref/STAR \
--genomeFastaFiles /home/ebarrozo/lum/docs/ref/STAR/genome.fa \
--sjdbGTFfile /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.gtf \
--sjdbOverhang 149 \
--genomeSAindexNbases 10

cd /home/ebarrozo/lum/data
## Align to custom ZIKV+hg38 genome (see zikv-hg38_eb.sh)
for f in *_1.trimmed_PE.fastq; do
STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/lum/docs/ref/STAR \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 3 \
--outFilterMismatchNmax 2 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignEndsType Extend5pOfRead1 \
--sjdbGTFfile /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.gtf \
--outFilterType BySJout \
--sjdbOverhang 149 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## 89% reads too short. Let's change filters. 

## From Aagaard HITS-CLIP: 
## 60-66% of reads too short. Let's change the filter outFilterScoreMinOverLread 0.66 to
#(0.2=0%) 0.5=62%, 0.3=14%. let's go with that filter.
## --outFilterMatchNminOverLread  0.1 \ 0.66 default

# Convert logs to txt files and save to local
for f in *Log.final.out; do
cat $f > ${f%.out}.txt
done

## Remove PCR duplicates
for f in *Aligned.sortedByCoord.out.bam; do
picard MarkDuplicates \
INPUT=$f \
OUTPUT=${f%Aligned.sortedByCoord.out.bam}.duprm.bam \
METRICS_FILE=${f%Aligned.sortedByCoord.out.bam}.met \
REMOVE_DUPLICATES=true
done

## Make counts files for biological replicates
for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/lum/results
done

bedtools bamtobed -i SRR10318622_1.trimmed_PE.fastq.duprm.sorted.bam > SRR10318622_1.trimmed_PE.fastq.duprm.sorted.bed

samtools sort -@ 32 -n SRR10318622_1.trimmed_PE.fastq.duprm.bam -o SRR10318622.sorted.bam
samtools index SRR10318622.sorted.bam
htseq-count -f bam -s no SRR10318622.sorted.bam /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.gtf > SRR10318622.htseq.counts.txt
head SRR10318622.htseq.counts.txt
tail SRR10318622.htseq.counts.txt
head SRR10318622_1.trimmed_PE.fastq.duprm.sorted.bed
tail SRR10318622_1.trimmed_PE.fastq.duprm.sorted.bed

cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs
mkdir dang
cd /home/ebarrozo/dang/docs/ref/STAR
cp mirbase.genes2.filtered.gtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/dang

cp mirbase.genes2.filtered.gtf /home/ebarrozo/lum/docs/ref

cd /home/ebarrozo/lum
docker start -i lum-docker

cd /home/ebarrozo/lum/data
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 SRR10318622.sorted.bam /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf --additional-attr=gene_name > SRR10318622.htseq.counts.2.txt
tail -n 50 SRR10318622.htseq.counts.2.txt
head -n 50 SRR10318622.htseq.counts.2.txt

tail /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf

head -n 3000 /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf
gene_id ENSG gene_name LINC01409
tail -n 50 /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf

samtools sort -@ 32 -n SRR10318622_1.trimmed_PE.fastq.duprm.bam -o SRR10318622.sorted.bam
samtools index SRR10318622.sorted.bam
rm SRR10318622.sorted.bam.bai

htseq-count -f bam -s no --max-reads-in-buffer=160000000000 SRR10318622.sorted.bam /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.3.gtf --additional-attr=gene_name > SRR10318622.htseq.counts.2.txt
__no_feature		3820823
__ambiguous		632
__too_low_aQual		0
__not_aligned		0
__alignment_not_unique		160568
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 SRR10318622.sorted.bam /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf --additional-attr=gene_name > SRR10318622.htseq.counts.2.txt
__no_feature		2486875
__ambiguous		83445
__too_low_aQual		0
__not_aligned		0
__alignment_not_unique		160568
htseq-count -f bam -s yes --max-reads-in-buffer=160000000000 SRR10318622.sorted.bam /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf --additional-attr=gene_name > SRR10318622.htseq.counts.2.txt
tail -n 50 SRR10318622.htseq.counts.2.txt
__no_feature		3631606
__ambiguous		2275
__too_low_aQual		0
__not_aligned		0
__alignment_not_unique		160568
head -n 50 SRR10318622.htseq.counts.2.txt

head -n 50 /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.gtf
tail -n 50 /home/ebarrozo/lum/docs/ref/STAR/fixed.genes.gtf

################################# Generate RPM-scaled bigwig files for UCSC

cd /home/ebarrozo/lum/data

ls *.duprm.bam > BAM.list
SRR10318591_1.trimmed_PE.fastq.duprm.bam
SRR10318592_1.trimmed_PE.fastq.duprm.bam
SRR10318596_1.trimmed_PE.fastq.duprm.bam
SRR10318597_1.trimmed_PE.fastq.duprm.bam
SRR10318598_1.trimmed_PE.fastq.duprm.bam
SRR10318599_1.trimmed_PE.fastq.duprm.bam
SRR10318600_1.trimmed_PE.fastq.duprm.bam
SRR10318601_1.trimmed_PE.fastq.duprm.bam
SRR10318605_1.trimmed_PE.fastq.duprm.bam
SRR10318606_1.trimmed_PE.fastq.duprm.bam
SRR10318607_1.trimmed_PE.fastq.duprm.bam
SRR10318608_1.trimmed_PE.fastq.duprm.bam
SRR10318609_1.trimmed_PE.fastq.duprm.bam
SRR10318610_1.trimmed_PE.fastq.duprm.bam
SRR10318614_1.trimmed_PE.fastq.duprm.bam
SRR10318615_1.trimmed_PE.fastq.duprm.bam
SRR10318616_1.trimmed_PE.fastq.duprm.bam
SRR10318617_1.trimmed_PE.fastq.duprm.bam
SRR10318618_1.trimmed_PE.fastq.duprm.bam
SRR10318619_1.trimmed_PE.fastq.duprm.bam
SRR10318621_1.trimmed_PE.fastq.duprm.bam
SRR10318622_1.trimmed_PE.fastq.duprm.bam

cd45(-)
nano healthy.list
SRR10318591_1.trimmed_PE.fastq.duprm.bam
SRR10318596_1.trimmed_PE.fastq.duprm.bam
SRR10318597_1.trimmed_PE.fastq.duprm.bam
SRR10318598_1.trimmed_PE.fastq.duprm.bam

nano patient1.list
SRR10318605_1.trimmed_PE.fastq.duprm.bam
SRR10318606_1.trimmed_PE.fastq.duprm.bam
SRR10318607_1.trimmed_PE.fastq.duprm.bam

nano patient2.list
SRR10318614_1.trimmed_PE.fastq.duprm.bam
SRR10318615_1.trimmed_PE.fastq.duprm.bam
SRR10318616_1.trimmed_PE.fastq.duprm.bam

nano patient3.list
SRR10318621_1.trimmed_PE.fastq.duprm.bam

samtools merge -f ctrl.merged.bam -b healthy.list
samtools index ctrl.merged.bam
samtools merge -f p1.merged.bam -b patient1.list
samtools index p1.merged.bam
samtools merge -f p2.merged.bam -b patient2.list
samtools index p2.merged.bam
samtools merge -f p3.merged.bam -b patient3.list
samtools index p3.merged.bam


cd45(+)
nano healthy.list
SRR10318592_1.trimmed_PE.fastq.duprm.bam
SRR10318599_1.trimmed_PE.fastq.duprm.bam
SRR10318600_1.trimmed_PE.fastq.duprm.bam
SRR10318601_1.trimmed_PE.fastq.duprm.bam

nano patient1.list
SRR10318617_1.trimmed_PE.fastq.duprm.bam
SRR10318618_1.trimmed_PE.fastq.duprm.bam
SRR10318619_1.trimmed_PE.fastq.duprm.bam

nano patient2.list
SRR10318614_1.trimmed_PE.fastq.duprm.bam
SRR10318615_1.trimmed_PE.fastq.duprm.bam
SRR10318616_1.trimmed_PE.fastq.duprm.bam

nano patient3.list
SRR10318622_1.trimmed_PE.fastq.duprm.bam

samtools merge -f ctrl.cd45p.merged.bam -b healthy.list
samtools index ctrl.cd45p.merged.bam
samtools merge -f p1.cd45p.merged.bam -b patient1.list
samtools index p1.cd45p.merged.bam
samtools merge -f p2.cd45p.merged.bam -b patient2.list
samtools index p2.cd45p.merged.bam
samtools merge -f p3.cd45p.merged.bam -b patient3.list
samtools index p3.cd45p.merged.bam

## For RPM calculations see star.logs and lum.meta.rpm.xlsx 

docker run -it -v $PWD:$PWD --name "barrozo1" continuumio/anaconda3 bash 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda bedtools 
conda install -c bioconda ucsc-bedgraphtobigwig 

cd /home/ebarrozo/lum/data
#docker start -i barrozo1
## For RPM calculations see star.logs and lum.meta.rpm.xlsx 


## Make a bigWig file from biological replicates merged
bedtools bamtobed -i ctrl.merged.bam > ctrl.merged.bed
bedtools genomecov -bg -strand + -5 -i ctrl.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.merged.plus.bedgraph > ctrl.merged.plus.sorted.bedgraph
bedGraphToBigWig ctrl.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i ctrl.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.merged.minus.bedgraph > ctrl.merged.minus.sorted.bedgraph
bedGraphToBigWig ctrl.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.merged.minus.sorted.bw
bedtools genomecov -scale 57.400923 -bg -strand + -5 -i ctrl.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.merged.rpm.plus.bedgraph > ctrl.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig ctrl.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 57.400923 -bg -strand + -5 -i ctrl.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.merged.rpm.minus.bedgraph > ctrl.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig ctrl.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p1.merged.bam > p1.merged.bed
bedtools genomecov -bg -strand + -5 -i p1.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.merged.plus.bedgraph > p1.merged.plus.sorted.bedgraph
bedGraphToBigWig p1.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p1.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.merged.minus.bedgraph > p1.merged.minus.sorted.bedgraph
bedGraphToBigWig p1.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.merged.minus.sorted.bw
bedtools genomecov -scale 6.862442 -bg -strand + -5 -i p1.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.merged.rpm.plus.bedgraph > p1.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p1.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 6.862442 -bg -strand + -5 -i p1.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.merged.rpm.minus.bedgraph > p1.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p1.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p2.merged.bam > p2.merged.bed
bedtools genomecov -bg -strand + -5 -i p2.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.merged.plus.bedgraph > p2.merged.plus.sorted.bedgraph
bedGraphToBigWig p2.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p2.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.merged.minus.bedgraph > p2.merged.minus.sorted.bedgraph
bedGraphToBigWig p2.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.merged.minus.sorted.bw
bedtools genomecov -scale 72.554865 -bg -strand + -5 -i p2.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.merged.rpm.plus.bedgraph > p2.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p2.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 72.554865 -bg -strand + -5 -i p2.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.merged.rpm.minus.bedgraph > p2.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p2.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p3.merged.bam > p3.merged.bed
bedtools genomecov -bg -strand + -5 -i p3.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.merged.plus.bedgraph > p3.merged.plus.sorted.bedgraph
bedGraphToBigWig p3.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p3.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.merged.minus.bedgraph > p3.merged.minus.sorted.bedgraph
bedGraphToBigWig p3.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.merged.minus.sorted.bw
bedtools genomecov -scale 24.571453 -bg -strand + -5 -i p3.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.merged.rpm.plus.bedgraph > p3.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p3.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 24.571453 -bg -strand + -5 -i p3.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.merged.rpm.minus.bedgraph > p3.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p3.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.merged.rpm.minus.sorted.bw

## Make a bigWig file from biological replicates merged
bedtools bamtobed -i ctrl.cd45p.merged.bam > ctrl.cd45p.merged.bed
bedtools genomecov -bg -strand + -5 -i ctrl.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.cd45p.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.cd45p.merged.plus.bedgraph > ctrl.cd45p.merged.plus.sorted.bedgraph
bedGraphToBigWig ctrl.cd45p.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.cd45p.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i ctrl.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.cd45p.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.cd45p.merged.minus.bedgraph > ctrl.cd45p.merged.minus.sorted.bedgraph
bedGraphToBigWig ctrl.cd45p.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.cd45p.merged.minus.sorted.bw
bedtools genomecov -scale 91.132685 -bg -strand + -5 -i ctrl.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.cd45p.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.cd45p.merged.rpm.plus.bedgraph > ctrl.cd45p.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig ctrl.cd45p.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.cd45p.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 91.132685 -bg -strand + -5 -i ctrl.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > ctrl.cd45p.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n ctrl.cd45p.merged.rpm.minus.bedgraph > ctrl.cd45p.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig ctrl.cd45p.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai ctrl.cd45p.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p1.cd45p.merged.bam > p1.cd45p.merged.bed
bedtools genomecov -bg -strand + -5 -i p1.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.cd45p.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.cd45p.merged.plus.bedgraph > p1.cd45p.merged.plus.sorted.bedgraph
bedGraphToBigWig p1.cd45p.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.cd45p.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p1.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.cd45p.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.cd45p.merged.minus.bedgraph > p1.cd45p.merged.minus.sorted.bedgraph
bedGraphToBigWig p1.cd45p.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.cd45p.merged.minus.sorted.bw
bedtools genomecov -scale 6.862442 -bg -strand + -5 -i p1.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.cd45p.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.cd45p.merged.rpm.plus.bedgraph > p1.cd45p.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p1.cd45p.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.cd45p.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 6.862442 -bg -strand + -5 -i p1.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p1.cd45p.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p1.cd45p.merged.rpm.minus.bedgraph > p1.cd45p.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p1.cd45p.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p1.cd45p.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p2.cd45p.merged.bam > p2.cd45p.merged.bed
bedtools genomecov -bg -strand + -5 -i p2.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.cd45p.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.cd45p.merged.plus.bedgraph > p2.cd45p.merged.plus.sorted.bedgraph
bedGraphToBigWig p2.cd45p.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.cd45p.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p2.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.cd45p.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.cd45p.merged.minus.bedgraph > p2.cd45p.merged.minus.sorted.bedgraph
bedGraphToBigWig p2.cd45p.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.cd45p.merged.minus.sorted.bw
bedtools genomecov -scale 23.003814 -bg -strand + -5 -i p2.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.cd45p.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.cd45p.merged.rpm.plus.bedgraph > p2.cd45p.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p2.cd45p.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.cd45p.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 23.003814 -bg -strand + -5 -i p2.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p2.cd45p.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p2.cd45p.merged.rpm.minus.bedgraph > p2.cd45p.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p2.cd45p.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p2.cd45p.merged.rpm.minus.sorted.bw

bedtools bamtobed -i p3.cd45p.merged.bam > p3.cd45p.merged.bed
bedtools genomecov -bg -strand + -5 -i p3.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.cd45p.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.cd45p.merged.plus.bedgraph > p3.cd45p.merged.plus.sorted.bedgraph
bedGraphToBigWig p3.cd45p.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.cd45p.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i p3.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.cd45p.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.cd45p.merged.minus.bedgraph > p3.cd45p.merged.minus.sorted.bedgraph
bedGraphToBigWig p3.cd45p.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.cd45p.merged.minus.sorted.bw
bedtools genomecov -scale 28.124362 -bg -strand + -5 -i p3.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.cd45p.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.cd45p.merged.rpm.plus.bedgraph > p3.cd45p.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig p3.cd45p.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.cd45p.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 28.124362 -bg -strand + -5 -i p3.cd45p.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > p3.cd45p.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n p3.cd45p.merged.rpm.minus.bedgraph > p3.cd45p.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig p3.cd45p.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai p3.cd45p.merged.rpm.minus.sorted.bw


for f in *.rpm.minus.sorted.bw; do
cp $f /home/ebarrozo/lum/results
done
for f in *.rpm.plus.sorted.bw; do
cp $f /home/ebarrozo/lum/results
done

docker stop barrozo1
docker rm barrozo1
## aagaard hits-clip tracks
## plus		https://genome.ucsc.edu/s/ebarrozo/aagaard%2Dzikv%2Dplus
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/D803E0E4-0795-4641-A61D-1A59C2AA3B81/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/87943640-A21D-43A9-94D8-5D22B922E849/mock.merged.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/4303E901-B7B1-4D1E-9BAF-044DF4E70C91/inf.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/BB623BB5-563D-429B-A445-FF7F5E947A1B/mock.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E733EFFE-C208-4A78-B739-E00921BD6619/inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/6C810DBE-95C8-4D02-A69D-008BAB596369/mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
## minus	https://genome.ucsc.edu/s/ebarrozo/aagaard%2Dzikv%2Dminus
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/A0DA8C90-7506-4F17-A75D-03255D96AB51/inf.merged.rpm.minus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/C218656F-3E67-4A4D-8E2E-F9448144C949/mock.merged.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/1FE9CF60-2FFE-41DD-A462-20B634826595/inf.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/648CA5F2-8E26-4780-86FC-A8CD879B174D/mock.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/F4026399-C1C6-42E8-AC7B-B7962DB618D5/inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E0062A69-E48E-47E1-B071-B9BD1E12F9EA/mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw

#### CD45-
## lum plus	https://genome.ucsc.edu/s/ebarrozo/lum%2Dcd45neg%2Dplus
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45+-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
## lum minus	https://genome.ucsc.edu/s/ebarrozo/lum%2Dcd45neg%2Dneg
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/459FC839-A779-4FB3-8B52-06B57BD1DB36/p1.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/9547C431-A235-46A4-94F3-AF898F213C1A/p2.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/00B08BC1-C9D9-4485-AF37-05A068ADEA60/p3.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/837720E3-12B2-45FD-AEB8-424448481E91/ctrl.merged.rpm.minus.sorted.bw

#### CD45+
## lum plus		https://genome.ucsc.edu/s/ebarrozo/lum%2Dcd45pos%2Dplus
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
## lum minus		https://genome.ucsc.edu/s/ebarrozo/lum%2Dcd45pos%2Dneg
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2A1248B5-EDA4-4074-8492-B2DCE233BDF0/p1.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/682B092C-D931-4C26-8D3E-9222C06AF35D/p2.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/7BC71040-18B4-422D-B41D-4A55CF25F3BE/p3.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/B5149479-00F6-4CEF-81E5-5A1ECE193128/ctrl.cd45p.merged.rpm.minus.sorted.bw



## final plus placenta https://genome.ucsc.edu/s/ebarrozo/zikv%2Dplacenta%2Dplus
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45+-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/D803E0E4-0795-4641-A61D-1A59C2AA3B81/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/87943640-A21D-43A9-94D8-5D22B922E849/mock.merged.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/4303E901-B7B1-4D1E-9BAF-044DF4E70C91/inf.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/BB623BB5-563D-429B-A445-FF7F5E947A1B/mock.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E733EFFE-C208-4A78-B739-E00921BD6619/inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/6C810DBE-95C8-4D02-A69D-008BAB596369/mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw

## final neg placenta https://genome.ucsc.edu/s/ebarrozo/zikv%2Dplacenta%2Dminus
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/459FC839-A779-4FB3-8B52-06B57BD1DB36/p1.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/9547C431-A235-46A4-94F3-AF898F213C1A/p2.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/00B08BC1-C9D9-4485-AF37-05A068ADEA60/p3.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/837720E3-12B2-45FD-AEB8-424448481E91/ctrl.merged.rpm.minus.sorted.bw
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2A1248B5-EDA4-4074-8492-B2DCE233BDF0/p1.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/682B092C-D931-4C26-8D3E-9222C06AF35D/p2.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/7BC71040-18B4-422D-B41D-4A55CF25F3BE/p3.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/B5149479-00F6-4CEF-81E5-5A1ECE193128/ctrl.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/A0DA8C90-7506-4F17-A75D-03255D96AB51/inf.merged.rpm.minus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/C218656F-3E67-4A4D-8E2E-F9448144C949/mock.merged.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/1FE9CF60-2FFE-41DD-A462-20B634826595/inf.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/648CA5F2-8E26-4780-86FC-A8CD879B174D/mock.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/F4026399-C1C6-42E8-AC7B-B7962DB618D5/inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E0062A69-E48E-47E1-B071-B9BD1E12F9EA/mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw

## Final placenta and hNPC plus (Aagaard + Lum + Dang)
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45+-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2651F398-F329-4B96-AFB3-1F688CD9DA26/p1.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/6FBF5F56-7889-4A6A-979E-F62610287E72/p2.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/636AF749-1767-4D49-B250-4AD829286156/p3.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/47521AAE-FF19-4DFE-A7D7-BC37A3B2255F/ctrl.cd45p.merged.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/D803E0E4-0795-4641-A61D-1A59C2AA3B81/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/87943640-A21D-43A9-94D8-5D22B922E849/mock.merged.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/4303E901-B7B1-4D1E-9BAF-044DF4E70C91/inf.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/BB623BB5-563D-429B-A445-FF7F5E947A1B/mock.merged.piranha.rpm.plus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E733EFFE-C208-4A78-B739-E00921BD6619/inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/6C810DBE-95C8-4D02-A69D-008BAB596369/mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw


track type=bigWig name="MR776-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/EE31FCF5-7D85-4E37-A8DD-CA58457330ED/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="hNPC-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="MR776-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/EE31FCF5-7D85-4E37-A8DD-CA58457330ED/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="hNPC-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="MR776-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/EE31FCF5-7D85-4E37-A8DD-CA58457330ED/inf.merged.rpm.plus.sorted.bw
track type=bigWig name="hNPC-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="MR776-miRseq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="Paraiba-miRseq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="hNPC-miRseq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="MR776-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="Paraiba-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw
track type=bigWig name="hNPC-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/ABA37365-B159-415D-8A7F-1E9CCA2EE347/uninf.merged.rpm.plus.sorted.bw



## Final placenta and hNPC minus (Aagaard + Lum + Dang)
track type=bigWig name="p1-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/459FC839-A779-4FB3-8B52-06B57BD1DB36/p1.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/9547C431-A235-46A4-94F3-AF898F213C1A/p2.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/00B08BC1-C9D9-4485-AF37-05A068ADEA60/p3.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45+RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/837720E3-12B2-45FD-AEB8-424448481E91/ctrl.merged.rpm.minus.sorted.bw
track type=bigWig name="p1-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/2A1248B5-EDA4-4074-8492-B2DCE233BDF0/p1.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p2-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/682B092C-D931-4C26-8D3E-9222C06AF35D/p2.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="p3-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/7BC71040-18B4-422D-B41D-4A55CF25F3BE/p3.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="ctrl-CD45-RNA-seq" bigDataUrl=https://de.cyverse.org/dl/d/B5149479-00F6-4CEF-81E5-5A1ECE193128/ctrl.cd45p.merged.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/A0DA8C90-7506-4F17-A75D-03255D96AB51/inf.merged.rpm.minus.sorted.bw
track type=bigWig name="Uninf-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/C218656F-3E67-4A4D-8E2E-F9448144C949/mock.merged.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/1FE9CF60-2FFE-41DD-A462-20B634826595/inf.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="Uninf-Piranha" bigDataUrl=https://de.cyverse.org/dl/d/648CA5F2-8E26-4780-86FC-A8CD879B174D/mock.merged.piranha.rpm.minus.sorted.bw
track type=bigWig name="ZIKV-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/F4026399-C1C6-42E8-AC7B-B7962DB618D5/inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
track type=bigWig name="Uninf-PureClip" bigDataUrl=https://de.cyverse.org/dl/d/E0062A69-E48E-47E1-B071-B9BD1E12F9EA/mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw




