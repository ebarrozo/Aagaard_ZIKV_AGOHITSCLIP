## aagaardhitsclip_de_v2.R
# Enrico Barrozo, BCM, Aagaard Lab, 2021

## A script analyzing unpublished hits-clip data generated in the aagaard lab

#Set the randomnness to a particular seed for reproducibility downstream
set.seed(seed=1)

# Load the required packages
library(edgeR)
library(DESeq2)
library(limma)
library(statmod)

##############################################################################################
##########################		Formatting Input Data		###########################
##############################################################################################

## Reading data
##Set the working directory to the folder containing the HTSeq .txt output files
setwd("/Users/enricobarrozo/Box/AagaardLab/aagaardhitsclip/results/v2")

#Manually change names from SRRXXXXXXXX.duprm.sorted.htseq.counts.txt to Mock_1.txt thru 5 and ZIKV 1 thru 4
## Manually fix gene names and add column names gene and counts; add miRNAs to column b, remove counts summary at the bottom, delete column A, add column names "gene" "counts"

#Calling the names of each files
Filenames <- dir( pattern = ".txt")

Names <- strsplit(x = Filenames, split = ".txt"); Names
Names <- unlist(Names)
Names

#Creating a list with all the data from each file in a vector
myfiles = (lapply(Filenames, function(x) read.delim(x, header = TRUE)))

#Merging all the data in df.merged
temp = myfiles[[1]]
for (i in 2:10) {
        temp = merge(temp, myfiles[[i]], by.x = "gene", by.y = "gene")
}
df.merged <- temp[-c(1:9),]
#Change the column names to the sample IDs
colnames(df.merged)[c(2:10)] <- Names

library(dplyr)
## Consolidate and sum duplicated rownames
df.merged <- df.merged %>% 
     group_by(gene) %>% 
     summarise_all(funs(sum))

row.names(df.merged) <- df.merged$gene
mycount <- data.frame(df.merged[,c(1:10)])
row.names(mycount) <- mycount$gene
mycount$gene = NULL

#Removing low counts
dim(mycount)
## 39222 genes
mycount <- mycount[rowSums(mycount) > 0, ]
dim(mycount)
## 8935 genes
mycount <- mycount[rowSums(mycount) > 10, ]
dim(mycount)
## 407 genes
mycount <- mycount[rowSums(mycount) > 30, ]
dim(mycount)
## 178 genes
#mycount <- mycount[rowSums(mycount) > 30, ]
#dim(mycount)
## 42 genes
#mycount <- mycount[rowSums(mycount) > 100, ]
#dim(mycount)
## 3 genes: 1-Mar, LINC01505, TBCE
#          Mock_1.2 Mock_1 Mock_2.1 Mock_3.1 Mock_3 Zika_1.2 Zika_1 Zika_2.1 Zika_2
#1-Mar            0    248        0        0      0        0      0        0      0
#LINC01505        0      0        0        0    256        0      0        0      0
#TBCE           256    256      512        0    256      512    512      512    768
#### Let's try removing genes with < 247 counts
mycount <- mycount[rowSums(mycount) < 247, ]
dim(mycount)
## 8932 genes excluding the 3 outliers above 

## Create annotation for colData
##	https://informatics.fas.harvard.edu/differential-expression-with-deseq2.html
sampleNames <- Names
condition <- c('Mock', 'Mock','Mock', 'Mock', 'Mock','ZIKV', 'ZIKV', 'ZIKV', 'ZIKV')
colData <- data.frame(row.names=colnames(mycount), condition=factor(condition, levels=c('ZIKV','Mock')))

dataset <- DESeqDataSetFromMatrix(countData = mycount,
                                  colData = colData,
                                  design = ~condition)

dds <- DESeq(dataset)
result <- results(dds, contrast=c('condition','ZIKV','Mock'))
result <- result[complete.cases(result),]  #remove any rows with NA
resultsNames(dds)

resOrdered <- result[order(result$pvalue),]
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

head(resOrdered)
n = 50
#n = 10
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 0, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < 0, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]
write.csv(topResults, file = "ZIKV_Mock_DESeq2_topResults.csv")
write.csv(resOrdered, file = "ZIKV_Mock_DESeq2_allResults.csv")

png(file=paste0("ZIKV_Mock_MAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotMA(resOrdered, main='RNA-seq', ylim=c(-2,2))
dev.off()

rld <- rlogTransformation(dds, blind=TRUE)

png(file=paste0("ZIKV_Mock_PCAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotPCA(rld)
dev.off()

png(file=paste0("ZIKV_Mock_sparsityplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotSparsity(dds)
dev.off()

png(file=paste0("ZIKV_Mock_countsplot.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(resOrdered$padj), intgroup='condition', pch = 19)
dev.off()

#BiocManager::install("RColorBrewer")
library('RColorBrewer')
hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
png(file=paste0("ZIKV_Mock_heatmap1.png"),
                res=300, 
                width=2500, 
                height=1500)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))
dev.off()

png(file=paste0("ZIKV_Mock_heatmap2.png"),
                res=300, 
                width=2500, 
                height=1500)
hmcol <- brewer.pal(11,'RdBu')
nCounts <- counts(dds, normalized=TRUE)
heatmap(as.matrix(nCounts[ row.names(topResults), ]), Rowv = NA, col = hmcol, mar = c(8,2))
dev.off()

png(file=paste0("ZIKV_Mock_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

png(file=paste0("ZIKV_Mock_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

png(file=paste0("ZIKV_Mock_p-hist_smooth.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue[resOrdered$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()

#BiocManager::install("genefilter")
library("genefilter")
library(gplots)
topVarGenes <- head(order(-rowVars(assay(rld))),35)

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$Mock ]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$Mock,"-",rld$condition)

png(file=paste0("ZIKV_Mock_clust-heatmap.png"),
                res=300, 
                width=2500, 
                height=1500)
#heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")
heatmap(mat, col=colors,
          labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

png(file=paste0("ZIKV_Mock_clust-heatmap.2.png"),
                res=300, 
                width=1500, 
                height=1500)
heatmap.2(mat, col=colors,
           labRow=FALSE, mar=c(10,2), scale="row")
dev.off()

res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$pvalue < 0.05, na.rm=TRUE)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
png(file=paste0("ZIKV_Mock_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

#BiocManager::install("pheatmap")
# this gives log2(n + 1):: transformation
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
ntd <- normTransform(dds)
#BiocManager::install("vsn")
library("vsn")
png(file=paste0("ZIKV_Mock_ntd.png"),
                res=300, 
                width=1500, 
                height=1500)
meanSdPlot(assay(ntd))
dev.off()
png(file=paste0("ZIKV_Mock_vsd.png"),
                res=300, 
                width=1500, 
                height=1500)
                meanSdPlot(assay(vsd))
dev.off()
png(file=paste0("ZIKV_Mock_rld.png"),
                res=300, 
                width=1500, 
                height=1500)
meanSdPlot(assay(rld))
dev.off()

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
png(file=paste0("ZIKV_Mock_vsd-heatmap.png"),
                res=300, 
                width=1500, 
                height=1500)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE)
dev.off()


png(file=paste0("ZIKV_Mock_vsd_heatmap1.png"),
                res=300, 
                width=2500, 
                height=1500)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE)
dev.off()
