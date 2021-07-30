# Dang_viral_HITS-CLIP.sh
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

cd /home/ebarrozo/dang/data/inf
mkdir viral
cd viral
cp /home/ebarrozo/dang/data/inf/SRR7062773.1.trimmed.fastq .
cp /home/ebarrozo/dang/data/inf/SRR7062774.1.trimmed.fastq .
cp /home/ebarrozo/dang/data/inf/SRR7062775.1.trimmed.fastq .
######### New alignment params for v4
cd /home/ebarrozo/dang/data/inf/viral
for f in *.trimmed.fastq; do
STAR --runMode alignReads --runThreadN 48 \
--genomeDir /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta \
--outFilterScoreMinOverLread 0.1 \
--outFilterMatchNmin 0.1 \
--outFilterMatchNminOverLread  0.1 \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterType BySJout \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMultimapNmax 5 \
--alignEndsType Local \
--sjdbGTFfile /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf \
--outSJfilterReads All \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn $f \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix ${f%.trimmed.fastq}
done

## minimal

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
--readFilesIn SRR7062773.1.trimmed.fastq \
--limitBAMsortRAM 1376155254 \
--outFileNamePrefix SRR7062773.1.min.

cat SRR7062773.1.min.Log.final.out
## Minimal looks very much like uninf, throwing out reads bc they're too short
## Adding bySJout filter. 1.38% uniquely mapped reads. 96.45 % thrown out bc they're too short. 
# --sjdbOverhang 74 reduced from 1.38 to 0.34%;; try to figure out too short filter again
--outFilterScoreMinOverLread 0.1 \ result- up to 1.38 % unique
--outFilterMatchNmin 0.1 \ no change
--outFilterMatchNminOverLread 0.1 \ ******** too short went to zero. 
try higher? --outFilterMatchNminOverLread 0.33 \ went to 64% too short. 
--outFilterMatchNminOverLread 0.25 \ --outFilterScoreMinOverLread 0.66 \ ;; went back to 0.34%. 
--outFilterMatchNminOverLread 0.1 \ --outFilterScoreMinOverLread 0.66 \ ;; went back to 0.34%. 
--outFilterMatchNminOverLread 0.2 \ --outFilterScoreMinOverLread 0.1 \ ;; too short went from 98 to 17%
--outFilterMultimapNmax 5 \ from 22 to 17% unique and 17% mapped to too many loci
--outFilterMismatchNoverReadLmax 0.04 \ from 0.16% too many mismatches to 0
--outFilterMismatchNmax 5 default 10 changed to 5; rm 0.04 filter too. results - no change
--outFilterMismatchNmax 2 ;; still 0
--alignEndsType Extend5pOfRead1 \ too many mismatches went to 14.8%. but too short went up to 44% from 23%. but 13% unique. Could this be it? 
--outFilterMatchNminOverLread 0.1 \--outFilterMultimapNmax 3 \ too short went to 3.2% and mismatch went to 16%. unique is 22%. 
## trying this filter in uninf. if there is ~20% in uninf alignment, increase 0.1 filter above back to 0.2
## it was the same in uninf... try increasing score filter. mapped length was 10.24 and should be higher
--outFilterScoreMinOverLread 0.33 \ 3.25% unique. 87% too short. 
--limitBAMsortRAM 1376155254 \	--outFilterScoreMinOverLread 0.66 \

picard MarkDuplicates \
INPUT=SRR7062773.1.min.Aligned.sortedByCoord.out.bam \
OUTPUT=SRR7062773.1.min.duprm.bam \
METRICS_FILE=SRR7062773.1.min..met \
REMOVE_DUPLICATES=true
samtools sort -@ 8 -n SRR7062773.1.min.duprm.bam  -o SRR7062773.1.min.sorted.duprm.bam 
htseq-count -f bam -s no SRR7062773.1.min.sorted.duprm.bam /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genes.gtf > SRR7062773.1.min.htseq.counts.txt


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

## Merge BAM files and index for pureclip peak calling
ls *.duprm.bam > BAM.list
samtools merge -f inf.merged.bam -b BAM.list
samtools index inf.merged.bam

## Make a bigWig file from biological replicates merged
bedtools bamtobed -i inf.merged.bam > inf.merged.bed
bedtools genomecov -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.plus.bedgraph > inf.merged.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.plus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.plus.sorted.bw

bedtools genomecov -bg -strand - -5 -i inf.merged.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.minus.bedgraph > inf.merged.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.minus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.minus.sorted.bw

############## Examine STAR final out logs to determine scaling factor used below
cat SRR7062773.1Log.final.out
cat SRR7062774.1Log.final.out
cat SRR7062775.1Log.final.out
##################### To normalize to reads per million (uniquely mapped reads/1,000,000)
# Determine the number of uniquely mapped reads, which can be extracted from the STAR [15] output log files.
## see Log.final.out Uniquely mapped reads number 
## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.rpm.plus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.plus.bedgraph > inf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.rpm.plus.sorted.bw

bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.minus.bedgraph > inf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.rpm.minus.sorted.bw

####################### Call peaks using piranha on combined bed
for f in *merged.bed;
do
Piranha $f -v -s -b 50 -p 0.05 -o ${f%.bed}.piranha.bed
done

## Make UCSC tracks of Piranha peaks using RPM scaling.
bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.piranha.rpm.plus.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.plus.bedgraph > inf.merged.piranha.sorted.bed
bedGraphToBigWig inf.merged.piranha.sorted.bed /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.piranha.rpm.sorted.bw
cp inf.merged.piranha.rpm.sorted.bw /home/ebarrozo/dang/results/viral

bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.piranha.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.piranha.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.piranha.rpm.minus.bedgraph > inf.merged.piranha.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.piranha.rpm.minus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.piranha.rpm.minus.sorted.bw
cp inf.merged.piranha.rpm.minus.sorted.bw /home/ebarrozo/dang/results/viral


## PureCLIP requires an older version of python. I have it installed on another conda "snakes".
## It also takes a reallly long time, so I run it in the background in a docker container

cd /home/ebarrozo/dang
docker stop dang-pureclip
docker rm dang-pureclip
docker run -it -v $PWD:$PWD --name "dang-pureclip"  zavolab/pureclip:1.0.5 bash 
pureclip -i /home/ebarrozo/dang/data/inf/viral/inf.merged.bam -bai /home/ebarrozo/dang/data/inf/viral/inf.merged.bam.bai -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa -ld -nt 24 -o /home/ebarrozo/dang/data/inf/viral/inf.merged.PureCLIP.crosslink_sites.bed -or /home/ebarrozo/dang/data/inf/viral/inf.merged.PureCLIP.crosslink_regions.bed
# ctrl p+q
## if it is running- check logs as they come out
docker logs -f dang-pureclip
#ctrl + z or c to exit
#or 
docker logs dang-pureclip

#### Remove 7th column of PureCLIP output file
cd /home/ebarrozo/dang/data/inf/viral
cat inf.merged.PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6 > inf.merged.PureCLIP.crosslink_sites_short.bed

## Make UCSC tracks of pureCLIP peaks using RPM scaling. 
bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.PureCLIP.crosslink_sites.rpm.plus.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.PureCLIP.crosslink_sites.rpm.plus.bedgraph > inf.merged.PureCLIP.crosslink_sites.sorted.bed
bedGraphToBigWig inf.merged.PureCLIP.crosslink_sites.sorted.bed /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.PureCLIP.crosslink_sites.rpm.sorted.bw
cp inf.merged.PureCLIP.crosslink_sites.rpm.sorted.bw /home/ebarrozo/dang/results/viral

bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.PureCLIP.crosslink_sites.bed -g /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai > inf.merged.PureCLIP.crosslink_sites.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.PureCLIP.crosslink_sites.rpm.minus.bedgraph > inf.merged.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.PureCLIP.crosslink_sites.rpm.minus.sorted.bedgraph /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa.fai inf.merged.PureCLIP.crosslink_sites.rpm.minus.sorted.bw
cp inf.merged.PureCLIP.crosslink_sites.rpm.minus.sorted.bw /home/ebarrozo/dang/results/viral

## change viral chromosome to chromosome 1
# Make the viral chromosome chr1, so it can be displayed on hg38
cp /home/ebarrozo/dang/docs/ref/NC_012532.1_genome/fasta/genome.fa .
sed -i 's/NC_012532.1/chr1/g' genome.fa
samtools faidx genome.fa
head genome.fa

cd /home/ebarrozo/dang/results/viral
cp /home/ebarrozo/dang/data/inf/viral/inf.merged.bed .

# Make the viral chromosome chr1, so it can be displayed on hg38
head inf.merged.rpm.plus.sorted.bedgraph
sed -i 's/NC_012532.1/chr1/g' inf.merged.bed

## ## 1.25.21_mergedagoclipseq_RPM.xlsx
bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/dang/results/viral/genome.fa.fai > inf.merged.rpm.plus.bedgraph 
head inf.merged.rpm.plus.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.plus.bedgraph > inf.merged.rpm.plus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.plus.sorted.bedgraph /home/ebarrozo/dang/results/viral/genome.fa.fai inf.merged.rpm.plus.sorted.bw

bedtools genomecov -scale 7.9 -bg -strand + -5 -i inf.merged.bed -g /home/ebarrozo/dang/results/viral/genome.fa.fai > inf.merged.rpm.minus.bedgraph 
export LC_COLLATE=C
sort -k1,1 -k2,2n inf.merged.rpm.minus.bedgraph > inf.merged.rpm.minus.sorted.bedgraph
bedGraphToBigWig inf.merged.rpm.minus.sorted.bedgraph /home/ebarrozo/dang/results/viral/genome.fa.fai inf.merged.rpm.minus.sorted.bw

