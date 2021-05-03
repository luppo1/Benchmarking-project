# Load all required packages
library("tximport")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("rhdf5")
library("gplots")
library("ggplot2")
library("GenomicFeatures")
library("affy")
library("REBayes")
library("reliaR")
library("Rmosek")
library("ggfortify")
library("sva")
library("apeglm")

#Set working directory
setwd("C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Hisat")

# Import counts data from kallisto into R.# path to the dataset.
dir = "C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Hisat"

#Import metadata
metadatah = read.csv(file= "metadata.csv", header=T, sep = ",")
head(metadatah)

# Reading counts data from featureCounts
counts <- read.csv("hisat2.csv", header=TRUE, sep = ",")
head(counts)

# Selecting columns of interest and removing irrelevant ones
countdata = counts[, c(7:32)]
countdata = countdata[, -c(2)]
head(countdata)

# Remove the Gene ID column
countdata <- countdata[, -c(1)]
head(countdata)

# Making "Geneid" column the rownames
rownames(countdata) <- counts[,1]
head(countdata)
colnames(countdata)

# Check if the metadata and samples have the same names
table(colnames(countdata)==metadatah$SampleID)

# Create the DESeqDataSet object from Matrix of counts and metadata
dds <- DESeqDataSetFromMatrix(countData = round(countdata), 
                              colData = metadatah, 
                              design = ~Condition)
nrow(dds) 

#Remove Genes with low counts
dds1 <- dds[ rowSums(counts(dds)) >= 100, ]
nrow(dds1)


#Run DESeq function on the data to perform differential gene expression analysis
dds1 <- DESeq(dds1)
head(assay(dds1))

#Building out results table
res_table <- results(dds1)
dim(res_table)
summary(res_table)

#MA plots
plotMA(res_table)

#How many adjusted p-values were less than 0.05?
sum(res_table$padj < 0.05, na.rm=TRUE)

#Working with alpha 0.05
res2 <- results(dds1, alpha=0.05)
dim(res2)
summary(res2)

#How many adjusted p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE)

# Select genes with p value less than 0.05
res_sig <- subset(res2, padj < 0.05)
dim(res_sig)

#We order our results table by the smallest p value:
res_small_p <- res_sig[order(res_sig$pvalue),]
dim(res_small_p)

#Write final list to file
write.csv(as.data.frame(res_small_p), "DESeq2_hisat.csv")
