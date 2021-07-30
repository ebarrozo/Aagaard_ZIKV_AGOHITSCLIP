###### zikv-hg38_eb_v2.sh
## 1.15.2021, Enrico Barrozo, Aagaard Lab, BCM

### v2: redoing the STAR indexing with shorter --genomeSAindexNbases 7 instead of the standard 14
### Also changing --sjdbOverhang from 100 to 74, matching the avgspotlength (75 minus 1) from the SRArunselecter

## V3 added miRbase miRNAs to gtf # https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgRna


## in order to make a custom human + zikv transcriptome, I utilized Cell Ranger
## I downloaded the ZIKV genome and annotation from NCBI NC_012532.1.gff3
## I manually modified the annotation (NC_012532.1.v2) such that each protein was separated into its own transcript
## IE: 5pUTR-anchC-M-E-NS1-NS2A-NS2B-NS3-NS4A-2k-NS4B-NS5-3pUTR
## I then combined this custom ZIKV annotation with the latest human transcriptome hg38
## using CellRanger mkgtf and mkref commands
## There are a few troubleshooting steps that included removing extraneous labels added by CellRanger
## Lastly, the genome is indexed using STAR 
## The output is a custom genome (genome.fa) and annotation (genes.gtf) prepared for STAR alignments
## This is referred to as GRCh38_and_NC_012532.1.v2

# convert gff3 to gtf
#If you need to convert .gff3 from ncbi to .gtf use the following commands
#module load gffread
#gffread NC_001806.gff3 -T -o NC_001806.gtf

cd /home/ebarrozo/dang/docs/ref
gffread GRCh38_latest_genomic.gff -T -o GRCh38_latest_genomic.gtf
gffread NC_012532.1.gff3 -T -o NC_012532.1.gtf

export PATH=/home/ebarrozo/cellranger-5.0.1:$PATH
cellranger mkgtf NC_012532.1.v2.gtf NC_012532.1.v2.filtered.gtf --attribute=gene_biotype:protein_coding
cellranger mkref \
--genome=NC_012532.1_genome \
--fasta=NC_012532.1.fasta \
--genes=NC_012532.1.v2.filtered.gtf

cellranger mkref --genome=GRCh38 --fasta=/home/ebarrozo/refdata-gex-GRCh38-2020-A/fasta/genome.fa --genes=/home/ebarrozo/refdata-gex-GRCh38-2020-A/genes/genes.gtf --nthreads=32 --genome=NC_012532.1.v2 --fasta=/home/ebarrozo/dang/docs/ref/NC_012532.1.fasta --genes=/home/ebarrozo/dang/docs/ref/NC_012532.1.v2.filtered.gtf

## A script to edit the GRCh38_and_NC_012532.1.v2 reference from CellRanger for STAR
cd /home/ebarrozo/dang/docs/ref
mkdir STAR
cd STAR
## Copy and paste the reference files to this directory 
cp /home/ebarrozo/dang/docs/ref/GRCh38_and_NC_012532.1.v2/fasta/genome.fa .
cp /home/ebarrozo/dang/docs/ref/GRCh38_and_NC_012532.1.v2/genes/genes.gtf .

track type=bigWig name="ZIKV-AGO-CLIP" bigDataUrl=https://de.cyverse.org/dl/d/12E73924-9429-4920-A5D5-64EF16C939B6/inf.merged.piranha.rpm.plus.sorted.bw


## Edit the gtf fil for dropseq by adding "transcript_name" in addition to "transcipt_id" using mike's script
sed 's/\r//g' genes.gtf \
|  awk -F $'\t' '{if(index($9,"transcript_id")!=0)\
{split($9,a,";");if(index(a[1],"transcript_id")!=0){changeVar=a[1];}\
else if(index(a[3],"transcript_id")!=0){changeVar=a[3]}else{changeVar=a[4]};};\
gsub("transcript_id","transcript_name",changeVar);\
print $0"; "changeVar;}' > fixed.genes.gtf

## Remove the prefixes "GRCh38_________" and "NC_012532.1.v2_" added by CellRanger
sed -i 's/GRCh38_________//g' genome.fa
sed -i 's/JN555585_1_10x_//g' genome.fa
sed -i 's/GRCh38_________//g' fixed.genes.gtf
sed -i 's/NC_012532.1.v2_//g' fixed.genes.gtf

## 	Generate a STAR genome build
cd /home/ebarrozo/dang/docs/ref/STAR
STAR \
--runThreadN 36 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/dang/docs/ref/STAR \
--genomeFastaFiles /home/ebarrozo/dang/docs/ref/STAR/genome.fa \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/STAR/fixed.genes.gtf \
--sjdbOverhang 75 \
--genomeSAindexNbases 7

### v2: redoing the STAR indexing with shorter --genomeSAindexNbases 7 instead of the standard 14
### Also changing --sjdbOverhang from 100 to 74, matching the avgspotlength (75 minus 1) from the SRArunselecter

## download hsa.gff3 from mirbase link here: # https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgRna
cd /home/ebarrozo/dang/docs/ref/mirbase
# convert to gtf
gffread hsa.gff3 -T -o hsa.mirbase.gtf
## tried going through gtf filter in cellranger but didn't work. try adding it manually. 

## IN R combine the fixed.genes.gtf and  hsa.miRbase.gtf
setwd("/Users/enricobarrozo/Box/AagaardLab/dang/docs")
library(tidyverse)
## for hsa.mirbase.2.txt; convert hsa.gff3 to gtf manually and include gene_id
mirbase<-as_tibble(read.table("hsa.mirbase.2.txt",sep = "\t",header = F))
fixed<-as_tibble(read.table("fixed.genes.txt",sep = "\t",header = F))
## remove viral genes
tail(fixed)
fixed <- fixed[-c(2765970:2765982),]
## remove extra chromosomes
fixed <- fixed[-c(2765347:2765969),]
res <- union_all(mirbase, fixed)
viral<-as_tibble(read.table("NC_012532.1.v2.filtered.txt",sep = "\t",header = F))
res <- union_all(res, viral)
#res$V1 = NULL
write.table(res, "mirbase.fixed.genes.txt", sep="\t")

cd /home/ebarrozo/dang/docs/ref/STAR
export PATH=/home/ebarrozo/cellranger-5.0.1:$PATH
cellranger mkgtf mirbase.genes2.gtf mirbase.genes2.filtered.gtf --attribute=gene_biotype:protein_coding

cellranger mkref \
--genome=NC_012532.1_genome \
--fasta=NC_012532.1.fasta \
--genes=NC_012532.1.v2.filtered.gtf






--sjdbGTFfeatureExon \

## Generate fasta file for hsa miRbase
grep "hsa-*" -A 1 mature.fa > hsa.mature.fa
sed -i 's/U/T/g' hsa.mature.fa

grep "exon" mirbase.fixed.genes.gtf > mirbase.genes.gtf
cut --complement -f1 mirbase.genes.gtf > mirbase.genes2.gtf
head mirbase.genes2.gtf
sed -i 's/Derives_from=/" mature_id "/g' mirbase.genes2.gtf


export PATH=/home/ebarrozo/cellranger-5.0.1:$PATH
cellranger mkgtf mirbase.genes2.gtf mirbase.genes2.filtered.gtf --attribute=gene_biotype:protein_coding


## run star index using mirbase.fixed.genes.gtf
## 	Generate a STAR genome build
cd /home/ebarrozo/dang/docs/ref/STAR
STAR \
--runThreadN 36 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/dang/docs/ref/STAR \
--genomeFastaFiles /home/ebarrozo/dang/docs/ref/STAR/genome.fa \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf \
--sjdbOverhang 74 \
--genomeSAindexNbases 7
