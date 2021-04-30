library("DESeq2")
library("tximport")
#IMPORTANTING THE TRANSCRIPT ABUNDANCE DATA
datakallisto<-list.files("quants_kallisto", pattern = "abundance.tsv", full.names = T, recursive = T) # Creating the vector for the files with a list
transciptdatakallisto<-lapply(X=datakallisto, FUN = function(x)read.table(x,header = T))  # readig the whole files using a anonymus fuction to then use tximport
listgrasshopperssamples<-read.table("RNAseq_DE_mapping_samples_list.txt", header = TRUE)
#reading the data and then generating counts with tximport
countskallisto<-tximport(files=datakallisto , type= "kallisto", txOut =TRUE, countsFromAbundance= "lengthScaledTPM")
all(colnames(countskallisto$counts)==row.names(datakallisto))
View(datakallisto)
head(datakallisto)
head(countskallisto)
str(listgrasshopperssamples)
ddsdata<-DESeqDataSetFromTximport(countskallisto, colData =listgrasshopperssamples, design= ~ morph)

summary(ddsdata)
str(ddsdata)

ddsdatare<-DESeq(ddsdata)
ddsresults<-results(ddsdatare)
ddsresults
ddsresults <- results(ddsdatare, name="morph_G_vs_B")
resultsNames(ddsdatare)
