## dang-counts_v2.sh

cd /home/ebarrozo/dang/data
mkdir agoclip
cd uninf
ls *.duprm.bam
## 70-72 uninf
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/agoclip
done
cd ../inf
ls *.duprm.bam
## 73-75
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/agoclip
done
cd /home/ebarrozo/dang/data/agoclip

for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done

head /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf
tail /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf

for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/dang/results/clipseq
done


cd /home/ebarrozo/dang/data/bulkrnaseq
mkdir all
cd MR766
ls *.duprm.bam
## 66-67
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/bulkrnaseq/all
done
cd ../Paraiba
ls *.duprm.bam
## 68-69
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/bulkrnaseq/all
done
cd ../uninf
ls *.duprm.bam
## 64-65
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/bulkrnaseq/all
done
cd /home/ebarrozo/dang/data/bulkrnaseq/all
for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/dang/results/bulkcounts
done



cd /home/ebarrozo/dang/data/mirnaseq

mkdir all
cd MR766
ls *.duprm.bam
## 76-79
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/mirnaseq/all
done
cd ../Paraiba
ls *.duprm.bam
## 80-81
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/mirnaseq/all
done
cd ../uninf
ls *.duprm.bam
## 76-77
for f in *.duprm.bam; do
cp $f /home/ebarrozo/dang/data/mirnaseq/all
done

cd /home/ebarrozo/dang/data/mirnaseq/all
ls *.duprm.bam
for f in *.duprm.bam; do
samtools sort -@ 32 -n $f -o ${f%.bam}.sorted.bam
done
for f in *duprm.sorted.bam; do
htseq-count -f bam -s no --max-reads-in-buffer=160000000000 --additional-attr=gene_name $f /home/ebarrozo/dang/docs/ref/STAR/mirbase.genes2.filtered.gtf > ${f%.bam}.htseq.counts.txt
done
for f in *.htseq.counts.txt; do
cp $f /home/ebarrozo/dang/results/mircounts
done

