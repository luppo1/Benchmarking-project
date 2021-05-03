# Load all required packages
library("tximport")
library("pheatmap")
library("RColorBrewer")
library("limma")
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
setwd("C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Star")

# Import counts data from kallisto into R.# path to the dataset.
dir = "C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Star"

#Import metadata
metadatah = read.csv(file= "metadata.csv", header=T, sep = ",")
head(metadatah)

# Reading counts data from featureCounts
counts <- read.delim("star.txt", header=TRUE, comment.char="#")
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
dds1 <- dds[ rowSums(counts(dds)) > 5, ]
nrow(dds1)

#Run DESeq function on the data to perform differential gene expression analysis
dds1 <- DESeq(dds1)
head(assay(dds1))

#Building out results table
res_table <- results(dds1)
summary(res_table)

#MA plots
plotMA(res_table)

#p-values and adjusted p-values
#We can order our results table by the smallest p value:
order_results <- res_table[order(res_table$pvalue),]
head(order_results)

#How many adjusted p-values were less than 0.1?
sum(order_results$padj < 0.1, na.rm=TRUE)

#Working with alpha 0.05
res2 <- results(dds1, alpha=0.05)
summary(res2)

#How many adjusted p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE)

#We order our results table by the smallest p value:
res_small_p <- res2[order(res2$pvalue),]

# Select genes with p less than 0.05
res_sig <- subset(res_small_p, padj < 0.05)
dim(res_sig)

#Write final list to file
write.table(as.data.frame(res_sig),"star_DESeq2.csv", row.names = T,sep=",",col.names=T)
