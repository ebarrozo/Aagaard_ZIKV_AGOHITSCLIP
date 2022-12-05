## recount_script.sh

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

# redoing bulk RNA-seq analysis of Aagaard Lab agohitsclip data using optimized tropho_bulkRNA-seq_preprocessing_v1.sh script

## analyze on aagaard2 Rstudio

## make a directory tree for this project
cd /home/ebarrozo/aagaardhitsclip
mkdir recount
cd /home/ebarrozo/aagaardhitsclip/recount
mkdir scripts
mkdir data
mkdir docs
mkdir results
cd results
mkdir qc
cd ..


cd /media/jochum00/Aagaard_Raid2/seferovi/Placental_Clip_Files/A23_unprocessed_sequence_files
for f in *txt.gz; do
cp $f /home/ebarrozo/aagaardhitsclip/recount/data
done
for f in *fastq.gz; do
cp $f /home/ebarrozo/aagaardhitsclip/recount/data
done
cd /home/ebarrozo/aagaardhitsclip/recount/data
ls


for f in *.txt.gz; do
mv $f ${f%.txt.gz}.fastq.gz 
done


############################################################################################################
####################################  Run FASTQC on 1 sample to determine if there are adapters and what the read length is  #############################################
############################################################################################################

cd /home/ebarrozo/aagaardhitsclip/recount/data

## open and view the file and write notes below (Safari doesn't work- so use Chrome or Firefox)
	## 51 read length
	## Illumina RNA PCR Primer, Index 1 detected

		## Make sure STAR generate genome with sjdb 50 and trimmomatic



############################################################################################################
####################################  custom human + ZIKV + mirTAR reference made in aagaardhitsclip_v1.sh  #############################################
############################################################################################################

## see STAR manual 2.7.0a https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
############################################################################################################
####################################  Run trimmomatic to remove adapters and filter reads #############################################
############################################################################################################
cd /home/ebarrozo/aagaardhitsclip/recount/data

## Added universaladapter to the adapter_v2.fa below
	## Sequences from Illumina documentation https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2019/03/illumina-adapter-sequences-2019-1000000002694-10.pdf

## Trim adapter and barcodes using trimmomatic	

# conda install trimmomatic
cd /home/ebarrozo/aagaardhitsclip/data
for f in *.fastq.gz; do
trimmomatic SE -threads 32 $f ${f%.fastq.gz}_1.trimmed_PE.fastq.gz SLIDINGWINDOW:4:15 MINLEN:20 ILLUMINACLIP:/media/jochum00/Aagaard_Raid2/ebarrozo/adapters_v2.fa:2:30:10
done
## These are not PE, but PE included for continuity in downstream scripts


## Check QC to see if adapters were removed
fastqc Zika_1_1.trimmed_PE.fastq.gz --threads 32 --outdir /home/ebarrozo/aagaardhitsclip/recount/results/qc
	## Read length is now: 

## Remove the files that won't be used



############################################################################################################
####################################  Use STAR to align to human reference (refdata-gex-GRCh38-2020-A)  #############################################
############################################################################################################
cd /home/ebarrozo/aagaardhitsclip/recount/data


# Align to custom ZIKV+hg38 genome (see zikv-hg38_eb.sh)
## Align to GCA_003118495.1 using mostly STAR with the standard ENCODE options (see pg 8 in the STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

for f in *_1.trimmed_PE.fastq.gz; do
STAR --runMode alignReads --runThreadN 32 \
--genomeDir /home/ebarrozo/aagaardhitsclip/docs/ref/STAR \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignEndsType Extend5pOfRead1 \
--sjdbGTFfile /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf \
--outFilterType BySJout \
--sjdbOverhang 50 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn <(gunzip -c $f) \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## Align to GCA_003118495.1 using STAR with the standard ENCODE options (see pg 8 in the STAR manual: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
#for f in *_1.trimmed_PE.fastq.gz; do
#STAR --runMode alignReads --runThreadN 40 \
#--genomeDir /home/ebarrozo/aagaardhitsclip/recount/docs \
#--outFilterType BySJout \
#--outFilterMultimapNmax 20 \
#--alignSJoverhangMin 8 \
#--alignSJDBoverhangMin 1 \
#--outFilterMismatchNmax 999 \
#--outFilterMismatchNoverReadLmax 0.04 \
#--alignIntronMin 20 \
#--alignIntronMax 1000000 \
#--alignMatesGapMax 1000000 \
#--sjdbGTFfile /media/jochum00/Aagaard_Raid2/ebarrozo/customrefs/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
#--sjdbOverhang 150 \
#--outSAMtype BAM SortedByCoordinate \
#--readFilesIn <(gunzip -c $f) \
#--outFileNamePrefix ${f%.trimmed.fastq}
#done

# Convert logs to txt files and save to local
for f in *Log.final.out; do
cat $f > ${f%.out}.txt
done


############################################################################################################
####################################  Remove PCR duplicates using picard rmdup  #############################################
############################################################################################################
## Remove PCR duplicates
# conda install picard

for f in *Aligned.sortedByCoord.out.bam; do
picard MarkDuplicates \
INPUT=$f \
OUTPUT=${f%Aligned.sortedByCoord.out.bam}.duprm.bam \
METRICS_FILE=${f%Aligned.sortedByCoord.out.bam}.met \
REMOVE_DUPLICATES=true
done

############################################################################################################
####################################  Use htseq2 to count transcripts  #############################################
############################################################################################################
## Make counts files for biological replicates
cd /home/ebarrozo/aagaardhitsclip/recount/data

# conda install samtools

for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done

## htseq requires gene_id in the each row of the .gtf file
# cd /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf
wc -l /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf
	# 3667358 genes total
# grep "exon" genes.gtf > filtered.genes.gtf
# wc -l filtered.genes.gtf
	## 
# grep "CDS" genes.gtf > filtered.genes.gtf
# wc -l filtered.genes.gtf
	## 
# grep "gene_id" genes.gtf > filtered.genes.gtf
# wc -l filtered.genes.gtf	
	## 2765969 lines with gene_id

cd /home/ebarrozo/aagaardhitsclip/recount/data

## removed -n which sorts by name not position
for f in *.duprm.bam; do
samtools sort -@ 32 $f -o ${f%.bam}.sorted.bam
done

for f in *duprm.sorted.bam; do
samtools index $f -@ 32
done

## added -r pos to specify alignmed sorted by position
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no -r pos --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/aagaardhitsclip/docs/ref/STAR/fixed.genes.gtf > ${f%.bam}.htseq.counts.txt
done


for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/aagaardhitsclip/recount/results
done


