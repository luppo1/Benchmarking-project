# Load all required packages
library("tximport")
library("pheatmap")
library("RColorBrewer")
library("limma")
library("edgeR")
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

#Creating the edgeR gene list
dge <- DGEList(countdata)
dim(dge)

#Create a design matrix
design <- cbind("1"=1,"1vs2"=rep(c(1,2), each = nrow(metadatah)/2))

# Removing genes that are lowly expressed
# Number of genes with 0 count in all samples 
table(rowSums(dge$counts==0)==24)

# The filterByExpr function in the edgeR package provides an automatic way to 
# filter genes, while keeping as many genes as possible with worthwhile counts.
keep.exprs <- filterByExpr(dge, group=metadatah$Condition)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)

# Apply scale normalization to counts, TMM normalization method performs 
# well for comparative studies.
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors

# Data transformation with CPM
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
#logCPM <- cpm(dge, log=FALSE, prior.count=2)

# Running the limma voom function
v <- voom(dge, design, plot=TRUE, normalize="quantile")

# After this, the usual limma pipelines for differential expression is be applied.
fit <- lmFit(v, design)
fit <- eBayes(fit)
res <- topTable(fit, coef=ncol(design),number=Inf)
dim(res)
res_pvalue <- as.data.frame(subset(res, adj.P.Val < 0.05))
dim(res_pvalue)

#We order our results table by the smallest p value:
order_res <- res_pvalue[order(res_pvalue$adj.P.Val),]
dim(order_res)

#Write final list to file
write.csv(order_res,"limma_star.csv")
