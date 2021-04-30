
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

ddsdatare<-DESeq(ddsdata)
ddsresults<-results(ddsdatare)
ddsresults
ddsresults <- results(ddsdatare, name="morph_G_vs_B")
resultsNames(ddsdatare)
