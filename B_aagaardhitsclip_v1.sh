## aagaardhitsclip_v1.sh
# Enrico Barrozo, BCM, Aagaard Lab, 2021

## Unpublished AGO-HITS-CLIP data from Max Seferovic and Mark Hamilton (2013)
## Uninfected or ZIKV-infected primary human trophoblast cultures

# login to aagaardlab2
aagaard2

# make aagaardhitsclip directory tree
mkdir aagaardhitsclip
cd aagaardhitsclip
mkdir scripts
mkdir data
mkdir docs
mkdir results
cd results
mkdir qc
cd ..

## Copy data to working directory
cd /home/ebarrozo/aagaardhitsclip/data
cd /media/jochum00/Aagaard_Raid2/seferovi/Placental_Clip_Files/A23_unprocessed_sequence_files
zcat A23_huTropho_Mock_1-2.txt.gz | head
zcat A23_huTropho_Mock_1.fastq.gz | head
zcat A23_huTropho_Mock_1.fastq.gz | tail

for f in *fastq.gz; do
cp $f /home/ebarrozo/aagaardhitsclip/data
done

cd /media/jochum00/Aagaard_Raid2/seferovi/Placental_Clip_Files/A23_unprocessed_sequence_files
for f in *txt.gz; do
cp $f /home/ebarrozo/aagaardhitsclip/data/txt
done
cd /home/ebarrozo/aagaardhitsclip/data/txt
for f in *.txt.gz; do
mv $f ${f%.txt.gz}.fastq.gz 
done

cd /media/jochum00/Aagaard_Raid2/seferovi/Placental_Clip_Files
for f in *fastq; do
cp $f /home/ebarrozo/aagaardhitsclip/data
done

## Make a docker container with all packages installed
cd /home/ebarrozo/aagaardhitsclip
docker pull continuumio/anaconda3
docker run -it -v $PWD:$PWD --name "aagaardhitsclipdocker" continuumio/anaconda3 bash 

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install trimmomatic
conda install fastqc
conda install star
conda install picard
conda install samtools
conda install -c bioconda htseq
conda install -c bioconda piranha 
conda install -c bioconda bedtools 
conda install -c biobuilds fastx-toolkit
conda install -c bioconda ucsc-bedgraphtobigwig 

#docker start -i aagaardhitsclipdocker
cd /home/ebarrozo/aagaardhitsclip/data
####################################################################################
#					Make custom reference genome 
####################################################################################
cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
gunzip refdata-gex-GRCh38-and-mm10-2020-A.tar.gz
tar refdata-gex-GRCh38-and-mm10-2020-A.tar
docker run -it -v $PWD:$PWD --name "cellranger" continuumio/anaconda3 bash 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c hcc cellranger 
cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs

docker pull cumulusprod/cellranger:3.0.2
docker run -it -v $PWD:$PWD --name "cellranger-sc" cumulusprod/cellranger:3.0.2 bash

cd /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs
cd mirbase.v22
# convert to gtf
gffread hsa.gff3 -T -o hsa.mirbase.gtf
grep "hsa-*" -A 1 mature.fa > hsa.mature.fa
sed -i 's/U/T/g' hsa.mature.fa
grep "exon" mirbase.fixed.genes.gtf > mirbase.genes.gtf
cut --complement -f1 mirbase.genes.gtf > mirbase.genes2.gtf
head mirbase.genes2.gtf
sed -i 's/Derives_from=/" mature_id "/g' mirbase.genes2.gtf

grep "exon" hsa.mirbase.2.gtf > hsa.mirbase.3.gtf
sed -i 's/;"/";/g' hsa.mirbase.3.gtf

## Filter gtf files to make sure they are formatted for cellranger. Extracts transcripts with "exon" and "protein_coding".
cellranger mkgtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mirbase.v22/hsa.mirbase.3.gtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mirbase.v22/hsa.mirbase.filtered.gtf --attribute=gene_biotype:protein_coding
cellranger mkgtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_latest_genomic.gtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_latest_genomic.filtered.gtf --attribute=gene_biotype:protein_coding
cellranger mkgtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/zikv/NC_012532.1.v2.gtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/zikv/NC_012532.1.v2.filtered.gtf --attribute=gene_biotype:protein_coding

## Add mirbase to the end of GRCh38 
cat /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_latest_genomic.filtered.gtf /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/mirbase.v22/hsa.mirbase.filtered.gtf > /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_mirbase.v22.gtf
## Make custom ref
cellranger mkref --nthreads=48 --genome=GRCh38 \
--fasta=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
--genes=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_mirbase.v22.gtf \
--genome=NC_012532.1.v2 \
--fasta=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/zikv/NC_012532.1.fa \
--genes=/media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/zikv/NC_012532.1.v2.filtered.gtf

## A script to edit the GRCh38_and_NC_012532.1.v2 reference from CellRanger for STAR
cd /home/ebarrozo/aagaardhitsclip/docs/ref
mkdir STAR
cd STAR
## Copy and paste the reference files to this directory 
cp /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_012532.1.v2/fasta/genome.fa .
cp /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/GRCh38_and_NC_012532.1.v2/genes/genes.gtf .

## Edit the gtf fil for dropseq by adding "transcript_name" in addition to "transcipt_id" using mike's script
sed 's/\r//g' genes.gtf \
|  awk -F $'\t' '{if(index($9,"transcript_id")!=0)\
{split($9,a,";");if(index(a[1],"transcript_id")!=0){changeVar=a[1];}\
else if(index(a[3],"transcript_id")!=0){changeVar=a[3]}else{changeVar=a[4]};};\
gsub("transcript_id","transcript_name",changeVar);\
print $0"; "changeVar;}' > fixed.genes.gtf

## Remove the prefixes "GRCh38_________" and "NC_012532.1.v2_" added by CellRanger
sed -i 's/GRCh38_________//g' genome.fa
sed -i 's/NC_012532.1.v2_//g' genome.fa
sed -i 's/GRCh38_________//g' fixed.genes.gtf
sed -i 's/NC_012532.1.v2_//g' fixed.genes.gtf

## 	Generate a STAR genome build
cd /home/ebarrozo/aagaardhitsclip/docs/ref/STAR
STAR \
--runThreadN 48 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/aagaardhitsclip/docs/ref/STAR \
--genomeFastaFiles /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa \
--sjdbGTFfile /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf \
--sjdbOverhang 50 

docker stop cellranger
docker rm cellranger
docker stop cellranger-sc
docker rm cellranger-sc
####################################################################################
###### Docker commands: 
# escape sequence and it will still run
ctrl + pq
## after it's done running
docker start -i aagaardhitsclipdocker
aagaardhitsclipdocker --help
cd /home/ebarrozo/aagaardhitsclip/data

## if it is running- check logs as they come out
docker logs -f aagaardhitsclipdocker
#ctrl + z or c to exit
#or 
docker logs aagaardhitsclipdocker
ctrl + z or c to exit
docker stop aagaardhitsclipdocker
docker rm aagaardhitsclipdocker

cd /home/ebarrozo/aagaardhitsclip/data
## Run QC on the fastq files
#fastqc --extract --nogroup --outdir /home/ebarrozo/aagaardhitsclip/results/qc SRR7062770.1.fastq
for f in *.fastq.gz; do
fastqc --threads 48 --outdir /home/ebarrozo/aagaardhitsclip/data $f
done

## Are adapters retained? I think so?
## What is the read length? 51
## Are reads PE or SE? SE (I think)


for f in *.fastq; do
fastqc --threads 48 --outdir /home/ebarrozo/aagaardhitsclip/data $f
done

## Trim adapter and barcodes using trimmomatic
cd /home/ebarrozo/aagaardhitsclip/data
for f in *.fastq.gz; do
trimmomatic SE -threads 48 $f ${f%.fastq.gz}.trimmed.fastq.gz SLIDINGWINDOW:4:20 MINLEN:20 ILLUMINACLIP:/home/ebarrozo/aagaardhitsclip/docs/adapters.fa:2:30:10
done

## see if adding the adapters found in the first fastqc run fixed the problem
for f in *.trimmed.fastq.gz; do
fastqc --threads 48 --outdir /home/ebarrozo/aagaardhitsclip/data $f
done
## What is the read length after trimming? 
## Are adapters retained? 

#	docker start -i aagaardhitsclipdocker
cd /home/ebarrozo/aagaardhitsclip/data

# Convert fastq to fasta
#for f in *.trimmed.fastq.gz; do
#zcat $f | fastq_to_fasta -n -r | gzip > ${f%.trimmed.fastq.gz}.fasta.gz
#done

# Unzip fastq files for alignment
for f in *trimmed.fastq.gz; do
gunzip $f
done

## Align to custom ZIKV+hg38 genome (see zikv-hg38_eb.sh)
for f in *trimmed.fastq; do
STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/aagaardhitsclip/docs/ref/STAR \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 3 \
--outFilterMismatchNmax 2 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignEndsType Extend5pOfRead1 \
--sjdbGTFfile /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf \
--outFilterType BySJout \
--sjdbOverhang 50 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## Using Lum_RNA-seq_v1 alignment settings
mkdir v2
for f in *trimmed.fastq; do
cp $f /home/ebarrozo/aagaardhitsclip/data/v2
done

for f in *trimmed.fastq; do
STAR --runMode alignReads --runThreadN 32 \
--genomeDir /home/ebarrozo/aagaardhitsclip/docs/ref/STAR \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 3 \
--outFilterMismatchNmax 2 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignEndsType Extend5pOfRead1 \
--sjdbGTFfile /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf \
--outFilterType BySJout \
--sjdbOverhang 50 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--outFileNamePrefix ${f%.trimmed.fastq}
done
### not working...
BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file A20_PTPE_12442_STARtmp//BAMsort/19/48
SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.


#inspect alignment logs
cat A19_SPT_12016Log.final.out
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
samtools sort -@ 8 -n $f -o ${f%.bam}.sorted.bam
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no $f /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf > ${f%.bam}.htseq.counts.txt
done

for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/aagaardhitsclip/results
done


## Make counts files for biological replicates
cd /home/ebarrozo/aagaardhitsclip/data
for f in *.sorted.bam; do
rm $f
done
for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done
for f in *.htseq.counts.txt; do
rm $f
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/aagaardhitsclip/results
done

cd /home/ebarrozo/aagaardhitsclip/data/txt
for f in *.sorted.bam; do
rm $f
done
for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done
for f in *.htseq.counts.txt; do
rm $f
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/lum/docs/ref/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/aagaardhitsclip/results/txt
done

####################################################################################
## bedGraphToBigWig not working in docker but working w/o docker.

## Merge BAM files and index for pureclip peak calling
ls *Zika_*.duprm.bam > BAM.list
samtools merge -f inf.merged.bam -b BAM.list
samtools index inf.merged.bam
#samtools fqidx genome.fa
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i inf.merged.bam > inf.merged.bed
bedtools genomecov -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.plus.bedgraph > inf.merged.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.minus.bedgraph > inf.merged.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.plus.bedgraph > inf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.minus.bedgraph > inf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.rpm.minus.sorted.bw
#############################################################################
## Merge BAM files and index for pureclip peak calling
ls *Mock_*.duprm.bam > BAM.list
samtools merge -f mock.merged.bam -b BAM.list
samtools index mock.merged.bam
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i mock.merged.bam > mock.merged.bed
bedtools genomecov -bg -strand + -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.plus.bedgraph > mock.merged.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.minus.bedgraph > mock.merged.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.rpm.plus.bedgraph > mock.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.rpm.minus.bedgraph > mock.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.rpm.minus.sorted.bw
#############################################################################
#############################################################################
## Merge BAM files and index for pureclip peak calling
ls *TPE_*.duprm.bam > BAM.list
samtools merge -f TPE.merged.bam -b BAM.list
samtools index TPE.merged.bam
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i TPE.merged.bam > TPE.merged.bed
bedtools genomecov -bg -strand + -5 -i TPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.plus.bedgraph > TPE.merged.plus.sorted.bedgraph
bedGraphToBigWig TPE.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i TPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.minus.bedgraph > TPE.merged.minus.sorted.bedgraph
bedGraphToBigWig TPE.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.rpm.plus.bedgraph > TPE.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig TPE.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.rpm.minus.bedgraph > TPE.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig TPE.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.rpm.minus.sorted.bw
#############################################################################
#############################################################################
## Merge BAM files and index for pureclip peak calling
ls *PTPE_*.duprm.bam > BAM.list
samtools merge -f PTPE.merged.bam -b BAM.list
samtools index PTPE.merged.bam
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i PTPE.merged.bam > PTPE.merged.bed
bedtools genomecov -bg -strand + -5 -i PTPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.plus.bedgraph > PTPE.merged.plus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i PTPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.minus.bedgraph > PTPE.merged.minus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.rpm.plus.bedgraph > PTPE.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.rpm.minus.bedgraph > PTPE.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.rpm.minus.sorted.bw
#############################################################################
#############################################################################
## Merge BAM files and index for pureclip peak calling
ls *TC_*.duprm.bam > BAM.list
samtools merge -f TC.merged.bam -b BAM.list
samtools index TC.merged.bam
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i TC.merged.bam > TC.merged.bed
bedtools genomecov -bg -strand + -5 -i TC.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.plus.bedgraph > TC.merged.plus.sorted.bedgraph
bedGraphToBigWig TC.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i TC.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.minus.bedgraph > TC.merged.minus.sorted.bedgraph
bedGraphToBigWig TC.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.rpm.plus.bedgraph > TC.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig TC.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.rpm.minus.bedgraph > TC.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig TC.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.rpm.minus.sorted.bw
#############################################################################
#############################################################################
## Merge BAM files and index for pureclip peak calling
ls *TC_*.duprm.bam > BAM.list
samtools merge -f SPT.merged.bam -b BAM.list
samtools index SPT.merged.bam
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i SPT.merged.bam > SPT.merged.bed
bedtools genomecov -bg -strand + -5 -i SPT.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.plus.bedgraph > SPT.merged.plus.sorted.bedgraph
bedGraphToBigWig SPT.merged.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.plus.sorted.bw
bedtools genomecov -bg -strand - -5 -i SPT.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.minus.bedgraph > SPT.merged.minus.sorted.bedgraph
bedGraphToBigWig SPT.merged.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.minus.sorted.bw
############## Examine STAR final out logs to determine scaling factor used below
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.rpm.plus.bedgraph > SPT.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig SPT.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.rpm.minus.bedgraph > SPT.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig SPT.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.rpm.minus.sorted.bw
#############################################################################


####################### Call peaks using piranha on combined bed
docker run -it -v $PWD:$PWD --name "aagaardhitsclip-piranha" continuumio/anaconda3 bash 

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda piranha 
conda install -c bioconda bedtools 
conda install -c bioconda ucsc-bedgraphtobigwig 

#docker start -i aagaardhitsclip-piranha
cd /home/ebarrozo/aagaardhitsclip/data

for f in *.merged.bed; do
Piranha $f -v -s -b 50 -p 0.05 -o ${f%.bed}.piranha.bed
done

docker stop aagaardhitsclip-piranha
docker rm aagaardhitsclip-piranha

## Make UCSC tracks of Piranha peaks using RPM scaling.
#############################################################################
#############################################################################
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.plus.bedgraph > inf.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.minus.bedgraph > inf.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.piranha.rpm.plus.bedgraph > mock.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.piranha.rpm.minus.bedgraph > mock.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.piranha.rpm.plus.bedgraph > TPE.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig TPE.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.piranha.rpm.minus.bedgraph > TPE.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig TPE.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.piranha.rpm.plus.bedgraph > PTPE.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.piranha.rpm.minus.bedgraph > PTPE.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.piranha.rpm.plus.bedgraph > TC.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig TC.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.piranha.rpm.minus.bedgraph > TC.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig TC.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.piranha.rpm.plus.bedgraph > SPT.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig SPT.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.piranha.rpm.minus.bedgraph > SPT.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig SPT.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.piranha.rpm.minus.sorted.bw
#############################################################################


## PureCLIP requires an older version of python. I have it installed on another conda "snakes".
## It also takes a reallly long time, so I run it in the background in a docker container

cd /home/ebarrozo/aagaardhitsclip
docker stop aagaardhitsclip-pureclip
docker rm aagaardhitsclip-pureclip
docker pull zavolab/pureclip:1.0.5
docker run -it -v $PWD:$PWD --name "aagaardhitsclip-pureclip"  zavolab/pureclip:1.0.5 bash
cd /home/ebarrozo/aagaardhitsclip/data

for f in *.merged.bam; do
pureclip -i $f -bai $f.bai -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa -ld -nt 24 -o ${f%.bed}.PureCLIP.crosslink_sites.bed -or ${f%.bed}.PureCLIP.crosslink_regions.bed
done

#pureclip -i /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam -bai /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam.bai -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fasta/genome.fa -ld -nt 24 -o /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam.PureCLIP.crosslink_sites.bed -or /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam.PureCLIP.crosslink_regions.bed
# docker start -i aagaardhitsclipdocker
# ctrl p+q
## if it is running- check logs as they come out
docker logs -f aagaardhitsclip-pureclip
#ctrl + z or c to exit
#or 
docker logs aagaardhitsclip-pureclip

#### Remove 7th coaagaardhitsclipn of PureCLIP output file
cd /home/ebarrozo/aagaardhitsclip/data
for f in *.merged.bam.PureCLIP.crosslink_sites.bed; do
#cat inf.merged.bam.PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6 > inf.merged.bam.PureCLIP.crosslink_sites_short.bed
cat $f | cut -f 1,2,3,4,5,6 > ${f%.merged.bam.PureCLIP.crosslink_sites.bed}.short.bed
done

## Make UCSC tracks of pureCLIP peaks using RPM scaling. 
#############################################################################
#############################################################################
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 0.504287 -bg -strand + -5 -i mock.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > TPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig TPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 1.865082 -bg -strand + -5 -i TPE.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > TPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig TPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 1.355666 -bg -strand + -5 -i PTPE.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai PTPE.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > TC.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig TC.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i TC.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > TC.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n TC.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > TC.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig TC.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai TC.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > SPT.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig SPT.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 3.031766 -bg -strand + -5 -i SPT.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > SPT.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n SPT.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > SPT.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig SPT.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai SPT.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################

## change viral chromosome to chromosome 1
# Make the viral chromosome chr1, so it can be displayed on hg38
cp /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fasta/genome.fa .
sed -i 's/NC_012532.1/chr1/g' genome.fa
samtools faidx genome.fa
head genome.fa

cd /home/ebarrozo/aagaardhitsclip/results
mkdir viral
cd viral
cp /home/ebarrozo/aagaardhitsclip/data/inf.merged.bed .
cp /home/ebarrozo/aagaardhitsclip/data/inf.merged.piranha.bed .
cp /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam.PureCLIP.crosslink_sites.bed .

# Make the viral chromosome chr1, so it can be displayed on hg38??
head inf.merged.rpm.plus.sorted.bedgraph
sed -i 's/NC_012532.1/chr1/g' inf.merged.bed

## ## 2.10.20_rpm-calcs.xlsx
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/results/viral/genome.fa.fai > inf.merged.rpm.plus.bedgraph 
head inf.merged.rpm.plus.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.plus.bedgraph > inf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/results/viral/genome.fa.fai inf.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 0.227706 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/results/viral/genome.fa.fai > inf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.minus.bedgraph > inf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/results/viral/genome.fa.fai inf.merged.rpm.minus.sorted.bw


###### Merge all files and make big wigs 
cd /home/ebarrozo/aagaardhitsclip/data
mkdir merged
cd merged
cp /home/ebarrozo/aagaardhitsclip/data/inf.merged.bam .
cp /home/ebarrozo/aagaardhitsclip/data/mock.merged.bam .
## Manually rename
cp /home/ebarrozo/aagaardhitsclip/data/txt/inf.merged.bam .
cp /home/ebarrozo/aagaardhitsclip/data/txt/mock.merged.bam .
ls *inf* > BAM.list
samtools merge -f inf.merged.bam -b BAM.list
samtools index inf.merged.bam
#samtools fqidx genome.fa
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i inf.merged.bam > inf.merged.bed
ls *mock* > BAM.list
samtools merge -f mock.merged.bam -b BAM.list
samtools index mock.merged.bam
#samtools fqidx genome.fa
## Make a bigWig file from biological replicates merged
bedtools bamtobed -i mock.merged.bam > mock.merged.bed


bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.plus.bedgraph > inf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.minus.bedgraph > inf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.rpm.minus.sorted.bw
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.rpm.plus.bedgraph > mock.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.rpm.plus.sorted.bw
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.rpm.minus.bedgraph > mock.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.rpm.minus.sorted.bw

for f in *.merged.bam; do
pureclip -i $f -bai $f.bai -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa -ld -nt 24 -o ${f%.bed}.PureCLIP.crosslink_sites.bed -or ${f%.bed}.PureCLIP.crosslink_regions.bed
done
for f in *.merged.bam.PureCLIP.crosslink_sites.bed; do
#cat inf.merged.bam.PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6 > inf.merged.bam.PureCLIP.crosslink_sites_short.bed
cat $f | cut -f 1,2,3,4,5,6 > ${f%.merged.bam.PureCLIP.crosslink_sites.bed}.short.bed
done
#############################################################################
bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.bedgraph > mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.bam.PureCLIP.crosslink_sites.rpm.plus.sorted.bw
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.bam.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.bedgraph > mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.bam.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
#############################################################################



for f in *.merged.bed; do
Piranha $f -v -s -b 50 -p 0.05 -o ${f%.bed}.piranha.bed
done
bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.plus.bedgraph > inf.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 0.341672 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > inf.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.minus.bedgraph > inf.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai inf.merged.piranha.rpm.minus.sorted.bw
#############################################################################
#############################################################################
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.piranha.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.piranha.rpm.plus.bedgraph > mock.merged.piranha.rpm.plus.sorted.bedgraph
bedGraphToBigWig mock.merged.piranha.rpm.plus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.piranha.rpm.plus.sorted.bw
bedtools genomecov -scale 0.827166 -bg -strand + -5 -i mock.merged.piranha.bed -g /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai > mock.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n mock.merged.piranha.rpm.minus.bedgraph > mock.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig mock.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/genome.fa.fai mock.merged.piranha.rpm.minus.sorted.bw


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
