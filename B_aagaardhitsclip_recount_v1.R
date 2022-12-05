## aagaardhitsclip_recount_v1.R
	## bulkRNA-seq analysis optimized in tropho-PAH_RNA-seq_script_v2.R

# Enrico Barrozo, Ph.D. 
# Postdoctoral Associate | Aagaard Lab
# Baylor College of Medicine | Texas Children’s Hospital
# Department of Obstetrics & Gynecology, Division of Maternal-Fetal Medicine

## Analyze data on AagaardLab2
# 


## aagaardhitsclip recounted htseq-outs are located here /home/ebarrozo/aagaardhitsclip/results
## see recount_job_v1.sh and recount_script_v1.sh

################################################################################
##########################    Data Wrangling in Terminal      ######################################################
################################################################################
## in Terminal make a copy of the data in case there is an error
cd /home/ebarrozo/aagaardhitsclip/recount/results
rm -r DESeq2_analysis_v1
mkdir DESeq2_analysis_v1
cp *.txt /home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1

cd /home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1
ls
## Run Sample_regex.sh in terminal 
  ## removes last 5 lines that are the counts summary
  ## Removes coaagaardhitsclipn A
for n in $(ls *.txt);
do
  echo "running $n";
  cat $n|cut -f2-3>$n.tmp.txt;
  head -n-5 $n.tmp.txt>$n;
  rm $n.tmp.txt;
done
  ## now there are 2 columns: genes and counts

## need to fix columnnames
	## copy A23_huTropho_Mock_1_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt to local
		# open excel, ctrl + o and select A23_huTropho_Mock_1_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt
			# while opening, make sure to select _text for each COLUMN so gene names won't be mixed up
			## create new.names.csv 

	## see new.names.csv for what I want to copy and paste into B[36373:39235]

## make a list of the file names we can call on later
ls *.txt > list.tsv

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 0_meta.R
## Go back to RStudio 
## Set the working directory
setwd("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1")

## load the libraries
library(tidyverse)
# library(mosaic) # install.packages("mosaic")
library(pheatmap)
library(vegan) # install.packages("vegan")
library(ggpubr) # install.packages("ggpubr") 	## on terminal sudo apt install cmake 
library(ggsci) # install.packages("ggsci")

## load the data using the list
df<-as_tibble(read.table("list.tsv",header = F,sep = "\t"))

# View the data
df  ## confirmed all 9 samples were processed
# write.table(x = df,file = "samples.tsv",sep = "\t",row.names = T)

## RNA-seq filenames    SRR10318591__1.trimmed_PE.fastq.duprm.sorted.htseq.counts.txt
## RNA-seq filenames    A23_huTropho_Mock_1-2_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt


# Clip file names ;; this was to use the prefix as the filenames. Does not apply to this project. 
df<-df%>%mutate(name=gsub("_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt","",V1))%>%select(-V1)
df  

## Need to do something about the NAs in the replicate numbers
# tally(~name,df)

#data.frame(df)

# df
## save the metadata
# write.table(df,"meta.tsv",sep = "\t",row.names = F)

################################################################################
##########################    Loading Data into RStudio and making Metadata Table   ######################################################
################################################################################ EB modified from MJ 1_import.R
#import the metadata df that we made in the previous script
getwd()
#meta<-as_tibble(read.table("meta.tsv",sep = "\t",header = T))

## from metadata using sampleID and Treatment. Added an X in front of 6000 samples for R
meta<-as_tibble(read.csv("/home/ebarrozo/aagaardhitsclip/recount/docs/meta.csv",sep = ",",header = T))
meta  ## 42 × 3
# Meta has SampleID, Treatment, rep
    ## Treatments are Uninfected, ZIKV, BaP, BbF, DBA, and Mixture
    ## 7 replicates each

#make a dataframe with a coaagaardhitsclipn that matches the names of all the count files in the results folder
df<-data.frame(file=list.files("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1"))
df
## remove list.tsv
df<-df%>%filter(file!="list.tsv")


# df<-df%>%filter(file!="all-PCA")

 # df<-df%>%separate(col = file,into = c("file"),sep = "_",remove = F)
## View df
df

dim(df) ## 42 x 2

colnames(df)<- "Files"
df
meta
	## Edit meta to make sure Files=Names order or make a file.names csv where names match order of df
	file.names <- read.csv(file="/home/ebarrozo/aagaardhitsclip/recount/docs/filenames.csv", header=F)
	file.names
colnames(file.names)<-"Names"
df
file.names
df <- bind_cols(df,file.names)
  ## Make sure rows match

df <-full_join(df,meta)

## Confirm merge worked ; rows must match dims above
dim(df) #  [1] 8 x  6


#start up a dataframe that can be used as a matrix for joining
# tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[1],"_1.trimmed_PE.fastq.gz.duprm.sorted.htseq.counts.txt"),sep = "\t",header = F))
  tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[1]),sep = "\t",header = F))

colnames(tmp)<-c("gene",df$SampleID[1])
# colnames(tmp)<-c("gene",df$Sample[1])
tmp
df2<-tmp

new.names <- read.csv(file="/home/ebarrozo/aagaardhitsclip/recount/docs/new.names.csv", header=F)
colnames(new.names)<-"gene"

tmp<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[1]),sep = "\t",header = F))
tmp2<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[2]),sep = "\t",header = F))
tmp3<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[3]),sep = "\t",header = F))
tmp4<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[4]),sep = "\t",header = F))
tmp5<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[5]),sep = "\t",header = F))
tmp6<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[6]),sep = "\t",header = F))
tmp7<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[7]),sep = "\t",header = F))
tmp8<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[8]),sep = "\t",header = F))
tmp9<-as_tibble(read.table(file = paste0("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/",df$Files[9]),sep = "\t",header = F))

colnames(tmp)<-c("gene",df$SampleID[1])
colnames(tmp2)<-c("gene",df$SampleID[2])
colnames(tmp3)<-c("gene",df$SampleID[3])
colnames(tmp4)<-c("gene",df$SampleID[4])
colnames(tmp5)<-c("gene",df$SampleID[5])
colnames(tmp6)<-c("gene",df$SampleID[6])
colnames(tmp7)<-c("gene",df$SampleID[7])
colnames(tmp8)<-c("gene",df$SampleID[8])
colnames(tmp9)<-c("gene",df$SampleID[9])

tmp$gene<- new.names$gene
tmp2$gene<- new.names$gene
tmp3$gene<- new.names$gene
tmp4$gene<- new.names$gene
tmp5$gene<- new.names$gene
tmp6$gene<- new.names$gene
tmp7$gene<- new.names$gene
tmp8$gene<- new.names$gene
tmp9$gene<- new.names$gene

df3 <- full_join(tmp,tmp2)
df3

df3 <- full_join(df3,tmp3)
df3 <- full_join(df3,tmp4)
df3 <- full_join(df3,tmp5)
df3 <- full_join(df3,tmp6)
df3 <- full_join(df3,tmp7)
df3 <- full_join(df3,tmp8)
df3 <- full_join(df3,tmp9)

df3 	# 54,691 × 10


#remove the redundant gene coaagaardhitsclipns
df4<-df3%>%select(gene,!contains("gene"))
#convert it to a matrix with the rownames being the gene
df5<-as.matrix(df4%>%select(-gene))
rownames(df5)<-df4$gene

rownames(df5)<-df3$gene

dir.create("all-PCA")
setwd("all-PCA")


#convert it to a matrix with the rownames being the gene
df4.1<-as.matrix(df4%>%select(-gene))
rownames(df4.1)<-df4$gene

colSums(df4.1)

# 10-7-22; looks like too few reads. Checked the STAR alignment outputs and 43% of reads unmapped because they were too short. 
	## addressing in recount_job_v2.sh and will check back later
	
# Mock2 Mock1 Mock3 Mock4 Mock5 ZIKV2 ZIKV1 ZIKV4 ZIKV3 
#  734   572   200   473   369   680   527   597   489 

col.sums <- colSums(df4.1)

#write the count matrix to a table
 write.table(x = df4.1,file = "aagaardhitsclip.count_matrix.tsv",sep = "\t",row.names = T)
  ## YAY we have a counts matrix :)
   write.table(x = col.sums,file = "aagaardhitsclip.sum.counts.tsv",sep = "\t",row.names = T)

################################################################################
#setwd("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v2")
## load the libraries
#library(tidyverse)
# library(mosaic) # install.packages("mosaic")
# library(pheatmap)
# library(vegan) # install.packages("vegan")
# library(ggpubr) # install.packages("ggpubr") 	## on terminal sudo apt install cmake 
# library(ggsci) # install.packages("ggsci")
# df4.1<-as.matrix(read.table(file = "aagaardhitsclip.count_matrix.tsv",sep = "\t",header = F, ,row.names = T))

################################################################################
##########################    DESeq2 Analysis and PCA Plot  ######################################################
################################################################################ EB modified from MJ 2_analysis.R

df5 <- df4.1
#convert the counts to a distance matrix
d<-vegdist(x = t(df5),method = "euclidean")

# Draw a plot for a non-vegan ordination (cmdscale).
dis <- vegdist(t(df5),method = "gower")
mds <- cmdscale(dis, eig = TRUE)

mds<-data.frame(mds$points)
mds
mds<-mds%>%mutate(Sample=rownames(mds))
mds.names <- mds$SampleID
df <- bind_cols(df,mds.names)

mds
df$...7
# df$Sample
df<-full_join(df,mds, by=c("...7"="Sample"))
df<-as_tibble(df)%>%
  mutate(Eig1=X1,Eig2=X2)%>%
  select(-c(X1,X2))%>%
  mutate(across(where(is.character),as.factor))


## mutate metadata table
    ### KA-26748, KA.26751, KA.26759, KA.26761 need to be removed
dd_meta<-meta%>%
  mutate(across(where(is.character),factor))
dd_meta




################################################################################
library(ggpubr)
ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Treatment",
          shape = "Treatment",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Treatment",
          shape = "Treatment",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")
a

## save PCA plot 
png(file=paste0("PCA_plot_default-Treatment.png"),
                res=300, 
                width=2000, 
                height=1500)
a
dev.off()

a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "SampleID",
          shape = "SampleID",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")
a
## save PCA plot 
png(file=paste0("PCA_plot_default-Sample.png"),
                res=300, 
                width=2500, 
                height=2000)
a
dev.off()


a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "TrimesterInfected",
          shape = "TrimesterInfected",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")
a
## save PCA plot 
png(file=paste0("PCA_plot_default-TrimesterInfected.png"),
                res=300, 
                width=2500, 
                height=2000)
a
dev.off()


a<-ggscatterhist(data = df,
          x = "Eig2",
          y = "Eig1",
          color = "Individual",
          shape = "Individual",
          ellipse = F,
          ellipse.alpha = 0,
          ellipse.type = "t",palette = "npg",rug = T,
          title = "PCA plot of sample distances by gene count")
a
## save PCA plot 
png(file=paste0("PCA_plot_default-Individual.png"),
                res=300, 
                width=2500, 
                height=2000)
a
dev.off()



# Add vertical and horizontal line to a ggscatterhist
plots <- ggscatterhist(df, x = "Eig2", y = "Eig1",
                       margin.params = list(theme= theme_pubr()),
                       color="Treatment",
                       palette = "aaas",print = FALSE)
plots
plots$sp<-plots$sp+
  stat_density_2d(aes(colour = Treatment,size=0.05),contour = T,alpha=0.2,size=0.05,
                  contour_var = "count",n = 25)+
  xlim(c(-0.25,0.25))+
  ylim(-0.4,0.4)

plots$xplot<-plots$xplot+theme_pubr()
plots$yplot<-plots$yplot+theme_pubr()
plots

## save custom PCA plot 
png(file=paste0("PCA_custom-plot.png"),
                res=300, 
                width=2000, 
                height=1500)
plots
# ggpar(p = plots,legend.title = "Group") ## this command will change the legend from Treatment to whatever you want
dev.off()


### Make a custom plot 2
# tally(~Treatment,df)
data.frame(df)
# customise plot
customised_plot <- 
  ggplot(df,mapping = aes(x = Eig2,y = Eig1, z=0, color=Treatment)) +
  geom_point()+
  geom_contour() +
scale_colour_brewer(palette = "Set1") +
  theme_pubr(legend = "bottom")+xlim(c(-0.25,0.25))+ylim(-0.3,0.3)+
  coord_fixed(ratio = 0.5, clip = "off")

 customised_plot

png(file=paste0("PCA_plot_simple.png"),
                res=300, 
                width=2000, 
                height=1500)
customised_plot
dev.off()

################################################################################
##########################    DESeq2 Load Object  ######################################################
################################################################################ EB modified from MJ help 7/20/2022
#Set the randomnness to a particular seed for reproducibility downstream
set.seed(seed=1)
# Load the required packages
library(DESeq2) # BiocManager::install("DESeq2")
library(dplyr)

## mutate metadata table
dd_meta<-meta%>%
  mutate(across(where(is.character),factor))
dd_meta

#add up all the duplicate transcripts
    ## remove # KA-26748, KA.26751, KA.26759, KA.26761
df4.1<-df4%>%
  group_by(gene)%>%
  summarise(across(everything(),sum))
df4.2<-df4.1  

df4.2
#convert 
df5.1<-as.matrix(df4.2%>%select(-gene))
df5.1
colSums(df5.1)
  ## 0 samples have 0 reads

## Add gene as a coaagaardhitsclipn name
rownames(df5.1)<-df4.1$gene


## Finally, load data and metadata into DESeq2
deseq_counts <- DESeqDataSetFromMatrix(countData = df5.1, 
#deseq_counts <- DESeqDataSetFromMatrix(countData = df4.1, 
                                       colData =dd_meta,
                                       design = ~Treatment) 

## duplicated rownames

getwd()

## Run a quick analysis w/o collapsing replicates
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Treatment","Uninfected","ZIKV"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
    # 63 up and 190 down;  low counts [2]     : 11273, 45%%

## Make an elbowplot to determine how many PCs in a PCA = > 85% of variance. This will determine cutreen=
## https://support.bioconductor.org/p/83626/
     library(matrixStats)
     
     #How to get PCA plot?
     ##how to obtain d.deseq was described in DESeq2 manual
     
     cds=estimateDispersions(dds)
     vsd=varianceStabilizingTransformation(cds)
     plotPCA(vsd,intgroup=c("Treatment"))
     p1 <-      plotPCA(vsd,intgroup=c("Treatment"))
     p1
     #How to get PCA scree plot?
     
     ## calculate the variance for each gene
     rv <- rowVars(assay(vsd))

     ## select the ntop genes by variance
     select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
     
     ## perform a PCA on the data in assay(x) for the selected genes
     pca <- prcomp(t(assay(vsd)[select,]))
     
     ## the contribution to the total variance for each component
     percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
     
     ##plot the "percentVar"
     scree_plot=data.frame(percentVar)
     scree_plot[,2]<- c(1:9)
     p2 <-      scree_plot[,2]<- c(1:9)
     

     colnames(scree_plot)<-c("variance","component_number")
     ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")
     p3 <-      ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")
     
png(file=paste0("PCAplot-manual.png"),
                res=300, 
                width=2500, 
                height=1500)
p1
dev.off()
png(file=paste0("PCAplot-scree_plot-default.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

png(file=paste0("PCAplot-scree_plot-custom.png"),
                res=300, 
                width=2500, 
                height=1500)
p3
dev.off()

rm(a)
rm(customised_plot)
rm(map)

## elbow plot results -> use 4-5 PCs for >85% of variance
save.image("aagaardhitsclip_data-v1.RData")

setwd("..")
###### DESeq2 will have several iterations 
#########################################################################################################
##########################    Iteration I:    Uninfected-vs-ZIKVTreatment Analysis       ###########################################
#########################################################################################################
########## Determine the transcript counts for all samples to determine if there are any 0's or outliers. 
setwd("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1")
# load("/home/ebarrozo/aagaardhitsclip/recount/results/DESeq2_analysis_v1/all-PCA/aagaardhitsclip_data-v1.RData")

dir.create("Uninfected-vs-ZIKVTreatment")
setwd("Uninfected-vs-ZIKVTreatment")

# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
dds <-DESeq(deseq_counts,parallel = T)
res <- results(dds, contrast=c("Treatment","Uninfected","ZIKV"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res)
  ## Whatever is last in the contrast option is what the FC value is attibuted to. Any DE with this analysis is attributed to the ZIKV 
  ## Wald hypothesis testing results in 63 up and 190 down significant genes. Also consider LTR 
write.csv(res, file = "UninfectedTreatment_ZIKV_DESeq2_allResults-Wald.csv")



# Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution with Likelihood ratio test hypothesis testing, instead of Wald, where we use the estimated standard error of a log2 fold change to test if it is equal to zero
dds.LTR <-DESeq(deseq_counts, test="LRT", reduced=~1, parallel = T)
res.LTR <- results(dds.LTR, contrast=c("Treatment","Uninfected","ZIKV"), alpha=0.05, cooksCutoff=FALSE) ## Filtering adj. p < 0.05
summary(res.LTR)
write.csv(res.LTR, file = "UninfectedTreatment_ZIKV_DESeq2_allResults-LTR.csv")
  ## LTR results in 60 up and 129 down sig genes  
  ## What is LTR? see https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html 
    #The LRT is comparing the full model to the reduced model to identify significant genes. The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula (not log2 fold changes). Essentially the LRT test is testing whether the term(s) removed in the ‘reduced’ model explains a significant amount of variation in the data

########### Visualize DE results
res <- res[complete.cases(res),]  #remove any rows with NA
resOrdered <- res[order(res$pvalue),]
resOrdered

res05<- sum(resOrdered$padj < 0.05, na.rm=TRUE)
res05 # 205 for padj, plenty enough to plot. If not, consider filtering by p-value as done below
# res05<- sum(resOrdered$pvalue < 0.05, na.rm=TRUE)
# res05 ## 2631 sig genes
head(resOrdered)
n = 50
topResults <- rbind( resOrdered[ resOrdered[,'log2FoldChange'] > 2, ][1:n,], 
                    resOrdered[ resOrdered[,'log2FoldChange'] < -2, ][n:1,] )
topResults[c(1:5,(2*n-4):(2*n)), c('baseMean','log2FoldChange','pvalue')]
topResults
write.csv(topResults, file = "UninfectedTreatment_ZIKV_DESeq2_topResults.csv")
write.csv(resOrdered, file = "UninfectedTreatment_ZIKV_DESeq2_allResults-ordered.csv")

########## Generate a heatmap of the counts with only the Uninfected and ZIKV Treatment samples + the significant genes
res2<-as_tibble(res,rownames = "gene")%>%filter(padj<0.05)
res2 ## 205 sig genes; but we want to plot something, so let's filter by pvalue
# res2<-as_tibble(res,rownames = "gene")%>%filter(pvalue<0.05)
# res2  ## 2,631 genes
## let's add fold-change filters
res3<-res2%>%filter(log2FoldChange>2)
res3  ## 19 genes
res3.1 <- head(res3, n=15)
res4<-res2%>%filter(log2FoldChange<(-2))
res4 ## 95 genes
res4.1 <- head(res4, n=15)
res5 <- add_row(res3.1, res4.1)
res5
## top 30 genes filtered

## replace the whole list of 104 sig genes by this list of the top 20
	res2 <- res5

library(corrplot)

res3<-df4%>%
  filter(gene%in%res2$gene)%>%
  pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
res3


## res4<-full_join(res3,df%>%select(Sample,Treatment,rep)%>%distinct_all())
res4<-right_join(res3,df, by=c("SampleID"="...7"))%>%select(SampleID,Treatment)%>%distinct_all()
res4
res5<-res4%>%filter(Treatment%in%c("Uninfected","ZIKV"))%>%distinct_all()
res5

df_log<-df4.2%>%
  mutate(across(.cols = -gene,log1p))
gene_counts<-data.frame(colSums(df4.2%>%select(-gene)))
gene_counts
library(mosaic)
favstats(gene_counts$colSums.df4.2.....select..gene..)
sig<-df_log%>%filter(gene%in%res2$gene)%>%pivot_longer(cols = -gene,names_to = "SampleID",values_to = "count")
# sig2<-full_join(sig,df%>%select(Sample,Treatment)%>%distinct_all())
sig2<-right_join(sig,df, by=c("SampleID"="...7"))%>%select(SampleID,Treatment,gene,count)%>%distinct_all()
sig2
sig3<-sig2%>%filter(Treatment%in%c("Uninfected","ZIKV"))%>%distinct_all()%>%
  group_by(SampleID,gene)%>%summarise(count=mean(count))%>%ungroup()
sig2  
sig3  ## sig3 has some NAs 
tally(~SampleID,sig3)
## we have one of the removed samples still
#  
sig3%>%filter(!is.na(gene))


top.genes <- sig3$gene
top.genes <- unique(top.genes)
top.genes.Uninfected.ZIKV <- top.genes

sig4<-sig3%>%pivot_wider(names_from = SampleID,values_from = count)%>%filter(!is.na(gene))
  # %>%filter(!is.na(gene))%>%filter(!is.na(gene)) %>%filter(!is.na(gene))
##### Make sure this table doesn't have the empty samples

h<-as.matrix(sig4%>%select(-gene))
h
rownames(h)<-sig4$gene
dim(h)  # 20 12
df_slim<-df%>%
  select("...7",Treatment)%>%
  distinct_all()%>%
  filter(Treatment%in%c("Uninfected","ZIKV"))%>%
  distinct_all()
df_slim
df
anot_col<-data.frame(df_slim)%>%select(Treatment)
rownames(anot_col)<-df_slim$"...7"

anot_col
library(RColorBrewer)

####### Default heatmap
pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(13),
         annotation_col = anot_col,
         scale = "none")

p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         annotation_col = anot_col,
         scale = "none")
png(file=paste0("Uninfected-VS-ZIKV_counts-heatmap_row-scaled_notcut.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()


p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 3,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         scale = "row")
p1

p2 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         scale = "row")
p2

## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p2
dev.off()
## Save with the code
png(file=paste0("counts-heatmap_row-scaled_cut3.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()

####### Heatmap cutting clusters
p1 <- pheatmap(mat = h,
         color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdYlBu")))(10),
         # annotation_col = anot_col,cutree_rows = 2,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         annotation_col = anot_col,cutree_rows = 1,cutree_cols = 5,border_color = "black", ## Can change cutree_rows and coaagaardhitsclips
         scale = "row")
p1

## Save with the code
png(file=paste0("counts-heatmap_row-scaled-cut5.png"),
                res=300, 
                width=3000, 
                height=2000)
p1
dev.off()
############# End Analysis with Michael 7/20/22

## Perform VST for visualization
vst <- vst(dds, blind=FALSE)

## PCA plot
pcaData <- plotPCA(vst, intgroup=c("Treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
p2
png(file=paste0("PCAplot.png"),
                res=300, 
                width=2500, 
                height=1500)
p2
dev.off()

############ Plot histogram with 1 gene. Can be modified to compare any gene of interest.
png(file=paste0("UninfectedTreatment_ZIKVTreatment_countsplot-padj.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$padj), intgroup='Treatment')
dev.off()

############ Plot histogram with 1 gene. Can be modified to compare any gene of interest.
png(file=paste0("UninfectedTreatment_ZIKVTreatment_countsplot-pvalue.png"),
                res=300, 
                width=2500, 
                height=1500)
plotCounts(dds, gene=which.min(res$pvalue), intgroup='Treatment')
dev.off()

### Shows dispersion of all data
png(file=paste0("UninfectedTreatment_ZIKVTreatment_DispEsts.png"),
                res=300, 
                width=2500, 
                height=1500)
plotDispEsts(dds)
dev.off()

### Shows how many significant genes there are based on p value for all data
png(file=paste0("UninfectedTreatment_ZIKVTreatment_p-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on padj value for all data
png(file=paste0("UninfectedTreatment_ZIKVTreatment_padj-hist.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$padj, breaks=20, col="grey50", border="white")
dev.off()

### Shows how many significant genes there are based on p value, basemean>1
png(file=paste0("UninfectedTreatment_ZIKVTreatment_p-hist_smooth.png"),
                res=300, 
                width=2500, 
                height=1500)
hist(resOrdered$pvalue[resOrdered$baseMean > 1], breaks=20, col="grey50", border="white")
dev.off()


# Make a basic volcano plot
#BiocManager::install("genefilter")
library("genefilter")
library(gplots) # BiocManager::install("gplots")
topVarGenes <- head(order(-rowVars(assay(vst))),35)
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ vst$ZIKV ]
mat <- assay(vst)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(vst$ZIKV,"-",vst$Treatment)

# res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
#reset par
par(mfrow=c(1,1))

png(file=paste0("UninfectedTreatment_ZIKVTreatment_volcanoplot.png"),
                res=300, 
                width=1500, 
                height=1500)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="", xlim=c(-5,30)))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
p3<- with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

p2 <- p3 + scale_x_break(c(5, 22))

dev.off()


# Enhanced Volcano plot
library(EnhancedVolcano) #   devtools::install_github('kevinblighe/EnhancedVolcano')
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

png(file=paste0("UninfectedTreatment_ZIKVTreatment_Enhancedvolcanoplot.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

## more options : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#only-label-key-variables
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    labSize=4.0,
    x = 'log2FoldChange',
     xlim = c(-7.5, 30),
    drawConnectors = TRUE,
    max.overlaps= 30,
    y = 'pvalue')
p3
png(file=paste0("UninfectedTreatment_ZIKVTreatment_Enhancedvolcanoplot_top9.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

## more options : https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#only-label-key-variables
p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    labSize=4.0,
    x = 'log2FoldChange',
     xlim = c(-5, 5),
    drawConnectors = TRUE,
    max.overlaps= 30,
    y = 'pvalue')
p3
png(file=paste0("UninfectedTreatment_ZIKVTreatment_Enhancedvolcanoplot_top.wo.tmsb15b.png"),
                res=300, 
                width=2500, 
                height=2500)
p3
dev.off()

install.packages("ggbreak")
library(ggbreak)

p3 <-  EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
p3

p2 <- p3 + scale_x_break(c(7, 22))+
    theme(legend.position="none") 

png(file=paste0("UninfectedTreatment_ZIKVTreatment_Enhancedvolcanoplot_clipx.png"),
                res=300, 
                width=2500, 
                height=2500)
p2
dev.off()

setwd("..")
