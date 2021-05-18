
rm(list = ls())
library("DESeq2")
library("tximport")


### Import kallisto pseudo-counts into memory 
# Creating the vector for the files with a list
datakallisto<-list.files("quants_kallisto", pattern = "abundance.tsv", full.names = T, recursive = T) 

# readig the whole files using a anonymus fuction to then use tximport
transciptdatakallisto<-lapply(X=datakallisto, FUN = function(x)read.table(x,header = T))  

# Read in sample meta-data design
listgrasshopperssamples<-read.table("RNAseq_DE_mapping_samples_list.txt", header = TRUE)

# reading the count data with tximport
countskallisto<-tximport(files=datakallisto , type= "kallisto", txOut =TRUE, countsFromAbundance= "lengthScaledTPM")

# sanity check to see if my column names and row names are OK, also if our data is read into R fine. 
all(colnames(countskallisto$counts)==row.names(datakallisto))
View(datakallisto)
head(datakallisto)
head(countskallisto)
str(listgrasshopperssamples)

### creating a object to give to r the information about the design and the information of each sample
ddsdata<-DESeqDataSetFromTximport(countskallisto, colData =listgrasshopperssamples, design= ~ morph)

summary(ddsdata)
str(ddsdata)
### using DESEQ2 to analyse the data creating the object 
ddsdatare<-DESeq(ddsdata)
### using the result function to obtain the statistical values  
ddsresults<-results(ddsdatare)
### the results show the base mean where is possible to see average of the normalized count values, log2fold is possible to see how much gene expression have change and the adjusted p value that shows false discovery rate and using contrast to estimated the comparisons  
ddsresults
ddsresults <- results(ddsdatare, name="morph_G_vs_B")
ddsresults <- results(ddsdatare, contrast=c("morph","G","B"))
summary(ddsresults)
resultsNames(ddsdatare)
### using LFC to visualize and ranking the genes using the shrinkage  effect size using apeglm which improves the estimator
library(apeglm)
resLFC <- lfcShrink(ddsdatare, coef="morph_G_vs_B", type="apeglm")
resLFC
### creating a MA-plot to see the log2fod changes 
plotMA(ddsresults, ylim=c(-5,5))
plotMA(resLFC, ylim=c(-2,2))





           
                                

