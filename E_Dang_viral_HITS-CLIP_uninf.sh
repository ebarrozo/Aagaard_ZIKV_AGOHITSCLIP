# Dang_viral_HITS-CLIP_uninf.sh
## Aligning HITSCLIP data to viral genome only

# login to aagaardlab3
ssh <>

## For details on custom ZIKV genome generation, see zikv-hg38_eb_v3.sh
## Originally, ZIKV genome and annotation from NCBI NC_012532.1.gff3
cd /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta
cp /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/genes/genes.gtf.gz .
gunzip genes.gtf.gz
#sed -i 's/gene-//g' genes.gtf

cd /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta
## 	Generate a STAR genome build
STAR \
--runThreadN 1 \
--runMode genomeGenerate \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--genomeFastaFiles /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--sjdbOverhang 74 \
--genomeSAindexNbases 5

cd /home/ebarrozo/dang/data/uninf
mkdir viral
cd viral
cp /home/ebarrozo/dang/data/uninf/SRR7062770.1.trimmed.fastq .
cp /home/ebarrozo/dang/data/uninf/SRR7062771.1.trimmed.fastq .
cp /home/ebarrozo/dang/data/uninf/SRR7062772.1.trimmed.fastq .
######### New alignment params for v4
cd /home/ebarrozo/dang/data/uninf/viral
for f in *.trimmed.fastq; do
STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNmin 0 \
--outFilterMatchNminOverLread  0.33 \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 2 \
--alignEndsType Local \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--outSJfilterReads All \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## changing --outFilterMatchNmin 0.1 to 0 and --outFilterMultimapNmax 5 to 2 and --outFilterMatchNminOverLread 0.66 to 0.33 and outFilterScoreMinOverLread 0.1 to 0.33

## minimal
for f in *.trimmed.fastq; do
STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--outFilterScoreMinOverLread 0.1 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## This is the short read filter: --outFilterMatchNminOverLread  0.1 \ 0.66 default
## Try these filters on the infected data and see if alignments are improved compared to uninf data

############## Examine STAR final out logs to determine scaling factor used below
cat SRR7062770.1Log.final.out
cat SRR7062771.1Log.final.out
cat SRR7062772.1Log.final.out

STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--outFilterScoreMinOverLread 0.66 \
--outFilterMatchNminOverLread 0.2 \
--outFilterMultimapNmax 3 \
--outFilterMismatchNmax 2 \
--outFilterType BySJout \
--sjdbOverhang 74 \
--alignEndsType Extend5pOfRead1 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn SRR7062770.1.trimmed.fastq \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix SRR7062770.1.min.
## Trying alignment filters optimized in infected samples. Unique reads in inf samples was 22%
cat SRR7062770.1.min.Log.final.out
--outFilterScoreMinOverLread 0.33 \ ;; 3% unique in inf samples. here---- 3.16%... 


picard MarkDuplicates \
INPUT=SRR7062770.1.min.Aligned.sortedByCoord.out.bam \
OUTPUT=SRR7062770.1.min.duprm.bam \
METRICS_FILE=SRR7062770.1.min..met \
REMOVE_DUPLICATES=true
samtools sort -@ 8 -n SRR7062770.1.min.duprm.bam  -o SRR7062770.1.min.sorted.duprm.bam 
htseq-count -f bam -s no SRR7062770.1.min.sorted.duprm.bam /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf > SRR7062770.1.min.htseq.counts.txt


STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--outFilterType BySJout \
--sjdbOverhang 74 \
--alignEndsType Extend5pOfRead1 \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn SRR7062770.1.trimmed.fastq \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix SRR7062770.1.

cat SRR7062770.1.Log.final.out

picard MarkDuplicates \
INPUT=SRR7062770.1.Aligned.sortedByCoord.out.bam \
OUTPUT=SRR7062770.1.duprm.bam \
METRICS_FILE=SRR7062770.1.met \
REMOVE_DUPLICATES=true
samtools sort -@ 8 -n SRR7062770.1.duprm.bam  -o SRR7062770.1.sorted.duprm.bam 
htseq-count -f bam -s no SRR7062770.1.sorted.duprm.bam /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf > SRR7062770.1.htseq.counts.txt

cat SRR7062770.1.htseq.counts.txt

## Remove PCR duplicates
for f in *Aligned.sortedByCoord.out.bam; 
do
picard MarkDuplicates \
INPUT=$f \
OUTPUT=${f%Aligned.sortedByCoord.out.bam}.duprm.bam \
METRICS_FILE=${f%Aligned.sortedByCoord.out.bam}.met \
REMOVE_DUPLICATES=true
done

## generate counts matrices
for f in *.duprm.bam; do
samtools sort -@ 8 -n $f -o ${f%.bam}.sorted.bam
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no $f /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf > ${f%.bam}.htseq.counts.txt
done
cp *.htseq.counts.txt /home/ebarrozo/dang/results/viral

## Merge BAM files and index 
ls *.duprm.bam > BAM.list
samtools merge -f uninf.merged.bam -b BAM.list
samtools index uninf.merged.bam

## change viral chromosome to chromosome 1
# Make the viral chromosome chr1, so it can be displayed on hg38
cp /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa .
sed -i 's/NC_012532.1/chr1/g' genome.fa
samtools faidx genome.fa
head genome.fa

cd /home/ebarrozo/dang/results/viral
cp /home/ebarrozo/dang/data/uninf/viral/uninf.merged.bed .

# Make the viral chromosome chr1, so it can be displayed on hg38
head uninf.merged.rpm.plus.sorted.bedgraph
sed -i 's/NC_012532.1/chr1/g' uninf.merged.bed

## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 7.9 -bg -strand + -5 -i uninf.merged.bed -g /home/ebarrozo/dang/results/viral/genome.fa.fai > uninf.merged.rpm.plus.bedgraph 
head uninf.merged.rpm.plus.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n uninf.merged.rpm.plus.bedgraph > uninf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig uninf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/dang/results/viral/genome.fa.fai uninf.merged.rpm.plus.sorted.bw

bedtools genomecov -scale 7.9 -bg -strand + -5 -i uninf.merged.bed -g /home/ebarrozo/dang/results/viral/genome.fa.fai > uninf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n uninf.merged.rpm.minus.bedgraph > uninf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig uninf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/dang/results/viral/genome.fa.fai uninf.merged.rpm.minus.sorted.bw


####################### Call peaks using piranha on combined bed
for f in *merged.bed;
do
Piranha $f -v -s -b 50 -p 0.05 -o ${f%.bed}.piranha.bed
done


