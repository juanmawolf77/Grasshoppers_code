

rm(list = ls())
library("DESeq2")
library("tximport")  


### Import kallisto pseudo-counts into memory 
# Creating the vector for the files with a list
datakallisto<-list.files("quants_kallisto", pattern = "abundance.tsv", full.names = T, recursive = T) 


# Read in sample meta-data design
listgrasshopperssamples<-read.table("RNAseq_DE_mapping_samples_list.txt", header = TRUE)

# reading the count data with tximport
countskallisto<-tximport(files=datakallisto , type= "kallisto", txOut =TRUE, countsFromAbundance= "lengthScaledTPM")

### Check files which have been read ###
head(datakallisto)
head(countskallisto)
str(listgrasshopperssamples)

### creating a object to give to r the information about the design and the information of each sample
ddsdata<-DESeqDataSetFromTximport(countskallisto, colData =listgrasshopperssamples, design= ~ morph)
ddsdata
summary(ddsdata)
str(ddsdata)
nrow(ddsdata)

###  removing low counts 
keep <- rowSums(counts(ddsdata)) >= 120
ddsdatafil<- ddsdata[keep,]
ddsdatafil
keep
summary(keep)

### using DESEQ2 to analyse the data creating the object 
ddsdataresults<-DESeq(ddsdatafil, parallel = TRUE)


### using the result function to obtain the statistical values  
ddsre<-results(ddsdataresults, parallel = TRUE)
### the results show the base mean where is possible to see average of the normalized count values, ### 
### log2fold is possible to see how much gene expression have change and the adjusted p value that  ###
### shows false discovery rate and using contrast to estimated the comparisons                      ###
ddsre
nrow(ddsre)
str(ddsre)
summary(ddsre)


### Using contrast argument
ddsre<- results(ddsdataresults, contrast=c("morph","G","B"))

summary(ddsre)

### visualization of all the genes
nrow(as.data.frame(ddsre))
mcols(ddsre, use.names=TRUE)




#### Diagnostic plots dispersion nd volcano plot
#Volcanoplot
library("ggrepel")
library("ggplot2") 
library("EnhancedVolcano")
png("Enhanced_VolcanoPlot_1.png", height = 100, width = 100)
EnhancedVolcano(ddsre, 
                lab = rownames(ddsre), 
                x="log2FoldChange", 
                y="padj", 
                FCcutoff = 1.0, 
                pCutoff =0.1,
                pointSize = 1.5, 
                labSize = 3.0, 
                ylim = c(0,10), 
                xlim = c(-7.5,7.5),
                colAlpha=1)
dev.off()


library("ggpubr")
png("Kallisto_maplot_Log2_MeanExp_vs_Log2FC.png", width = 1600, height = 1200)
maplot = ggmaplot(ddsre, fdr = 0.01, fc = 1, size = 1,
                  palette = c("#e55c30", "#84206b", "#f6d746"),
                  genenames = as.vector(row.names(ddsre)),
                  legend="top", top = 15,font.label = c("bold", 11),
                  label.rectangle = TRUE, font.legend = c("bold",12), font.main = "bold",
                  xlab = "Log2 Mean Expression",  ylab="Log2 FC")
show(maplot)
dev.off()

#Dispersionplot
### Plot this DispEst plot into a png or pdf for presentation/publication
plotDispEsts(ddsdataresults)




###using LFC to visualize and ranking the genes using the shrinkage  effect size using apeglm which improves the estimator
library("apeglm")
resLFC <- lfcShrink(ddsdataresults, coef="morph_G_vs_B", type="apeglm", parallel = TRUE)
resLFC
summary(resLFC)
### creating a MA-plot to see the log2fold changes 
plotMA(resLFC, ylim = c(-15, 15))
plotMA(ddsre, ylim = c(-15, 15))
png("MAplotLFC_1.png", height = 10, width = 10)
maplot = ggmaplot(resLFC, fdr = 0.01, fc = 1, size = 1,
                  palette = c("#e55c30", "#84206b", "#f6d746"),
                  genenames = as.vector(row.names(resLFC)),
                  legend="top", top = 15,font.label = c("bold", 11),
                  label.rectangle = TRUE, font.legend = c("bold",12), font.main = "bold",
                  xlab = "Log2 Mean Expression",  ylab="Log2 FC")

show(maplot)
dev.off()
###  
png("Kallisto_maplot_Log2_MeanExp_vs_Log2FC.png", width = 1600, height = 1200)
EnhancedVolcano(resLFC, 
                lab = rownames(resLFC), 
                x="log2FoldChange", 
                y="padj", 
                FCcutoff = 1.0, 
                pCutoff =0.8,
                pointSize = 1.5, 
                labSize = 3.0, 
                ylim = c(0,10), 
                xlim = c(-8,8),
                colAlpha=1)
dev.off()




############### normalizing the data for the pheatmap with the vst method
vsd <- vst(ddsdataresults)



############### crateing a pheatmap 
library(pheatmap)
genes <- order(resLFC$log2FoldChange, decreasing=TRUE)[1:20] ############### selecting 20 genes with the largest postiive logfold change
############### plotting the pheatmap
pheatmap(assay(vsd)[genes, ], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)
############### plotting the pheatmap
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)[genes, ]), decreasing = TRUE), 20)
mat  <- assay(vsd)[ genes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("morph","sex")])
pheatmap(mat, annotation_col = anno)

############### the most significant gene is the transcript_84434
topGene <- rownames(ddsre)[which.min(ddsre$padj)]
plotCounts(ddsdataresults, gene = topGene, intgroup=c("morph"))

##### Loop this to plot all 20 top genes, and make a panel for the presentation ####
top20Genes<-rownames(ddsre)[order(ddsre$padj)[1:20]]
plotCounts(ddsdataresults, gene = top20Genes[2], intgroup = c("morph"))
par(mfrow=c(3, 3))
for(i in 1:10)
  plotCounts(ddsdataresults, gene = top20Genes[i], intgroup = c("morph"))
####### Annotation result making the list of the most expressed genes  

library("AnnotationDbi")

####### 100 list of most expressed genes with reLFC

resOrdered <- resLFC[order(resLFC$padj),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
resOrderedDF
write.csv(resOrderedDF, file = "resultsLFC.csv")

####### 100 list  of most expressed genes with ddsre
resOrdered2 <- ddsre[order(ddsre$padj),]
head(resOrdered2)

resOrderedDF2 <- as.data.frame(resOrdered2)[1:100, ]
resOrderedDF2
write.csv(resOrderedDF2, file = "resultsddsre.csv")


plotMA(ddsre)


#### Next steps!!!!! #####
#### Download the sequence and annotation files from dropbox, run THIS ONLY ONCE! ####


download.file("https://www.dropbox.com/s/kb5y5cyn7slazzf/good.Renamed_Gsib_trans_draft3_c95.fasta?dl=1",destfile = "Gsib_transcriptome_draft3_seqs.fasta")
download.file("https://www.dropbox.com/s/je2qi7iasg4vjzh/Gsib_trans_draft3_annotation_formatted.gff3?dl=1",destfile = "Gsib_transcriptome_draft3.gff3")
download.file("https://www.dropbox.com/s/q5apm1nupwzrps2/Gsib_transcriptome_draft3.gff3?dl=0",destfile = "Gsib_trans_draft3_annotation_formatted.gff3")
#### Read in the annotation file #####

##### 

library(Biostrings)
library(rtracklayer)
library(GenomicRanges)



### Import seqs and annotations
### Create a GRanges object from gff3 file import (of the transcriptome) #####
Gsib_annotation<-import.gff3(con = "Gsib_transcriptome_draft3.gff3")
Gsib_annotation
###

### Read transcriptome sequences into R ###
Gsib_sequences<-readDNAStringSet("Gsib_transcriptome_draft3_seqs.fasta")
Gsib_sequences
### Import annotations as dataframe
### Using brute force to import a GFF3 file into a dataframe, risky, but this is quick and dirty ###
Gsib_annotation_table<-readGFF("Gsib_transcriptome_draft3.gff3")

### Get names of all transcripts ###
names(Gsib_sequences)
seqnames(Gsib_sequences)
###### Write FASTA file ####
writeXStringSet(Gsib_sequences,"your_filename",format="fasta")
Gsib_sequences




###### getting info from annotattion using Genomic ranges

###### basic granges list commands

seqnames(Gsib_annotation)
ranges(Gsib_annotation)
strand(Gsib_annotation)
length(Gsib_annotation)
granges(Gsib_annotation)

names(Gsib_sequences)

elementNROWS(Gsib_sequences)
mcols(Gsib_annotation) ###### Annotation cordinates can be extracted with mcols

coverage(Gsib_annotation)

seqlengths(Gsib_annotation)





      