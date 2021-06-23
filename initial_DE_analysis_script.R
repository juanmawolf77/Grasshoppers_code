
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

### The summary shows lots of transcripts with low counts. Can we remove them ?
### I used this function to reduce the low counts to 7% 
keep <- rowSums(counts(ddsdatare)) >= 10
ddsdatare <- ddsdatare[keep,]

### using DESEQ2 to analyse the data creating the object 
ddsdatare<-DESeq(ddsdata)
### using the result function to obtain the statistical values  
ddsresults<-results(ddsdatare)
### the results show the base mean where is possible to see average of the normalized count values, ### 
### log2fold is possible to see how much gene expression have change and the adjusted p value that  ###
### shows false discovery rate and using contrast to estimated the comparisons                      ###
ddsresults

### In the case of continuous variables, use the name argument
#ddsresults <- results(ddsdatare, name="morph_G_vs_B")
### Using contrast argument
ddsresults <- results(ddsdatare, contrast=c("morph","G","B"))
summary(ddsresults)


### using LFC to visualize and ranking the genes using the shrinkage  effect size using apeglm which improves the estimator
library(apeglm)
resLFC <- lfcShrink(ddsdatare, coef="morph_G_vs_B", type="apeglm")
resLFC
summary(resLFC)
### creating a MA-plot to see the log2fold changes 
### summary(resLFC) and summary(ddsresults) are conflicting! we need figure out where we are going wrong ###
par(mfrow=c(2,1))
plotMA(ddsresults, ylim=c(-10,10))
plotMA(resLFC, ylim=c(-10,10))
DESeq2::plotMA(object = ddsresults, ylim = c(-5, 5))
plotMA(ddsdatare)

sum(ddsresults$padj < 0.1, na.rm=TRUE)
print(ddsdata)


#### Next steps!!!!! #####
#### Download the sequence and annotation files from dropbox, run THIS ONLY ONCE! ####


download.file("https://www.dropbox.com/s/kb5y5cyn7slazzf/good.Renamed_Gsib_trans_draft3_c95.fasta?dl=1",destfile = "Gsib_transcriptome_draft3_seqs.fasta")
download.file("https://www.dropbox.com/s/je2qi7iasg4vjzh/Gsib_trans_draft3_annotation_formatted.gff3?dl=1",destfile = "Gsib_transcriptome_draft3.gff3")

#### Read in the annotation file #####

##### i am having issues in this part i can not install the packages
BiocManager::install("Biostrings")
library(Biostrings)
BiocManager::install("rtracklayer")
library(rtracklayer)
BiocManager::install("GenomicRanges")
library(GenomicRanges)

### Import seqs and annotations
Gsib_annotation<-import.gff3(con = "Gsib_transcriptome_draft3.gff3")

Gsib_sequences<-readDNAStringSet("Gsib_transcriptome_draft3_seqs.fasta")

### Import annotations as dataframe

### it does not let me to import the  annotations
Gsib_annotation_table<-readGFF("Gsib_trans_draft3_annotation_formatted.gff3")



           
                                

