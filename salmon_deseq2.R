# First load all the libraries that will be needed for the analysis.
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

# Set working directory
setwd("C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Salmon")

# Import counts data from kallisto into R.# path to the dataset.
dir = "C:/Users/Albert Doughan/Desktop/Bioinformatics review/R Data/Salmon"

#Import metadata
metadatah <- read.csv(file= "metadata.csv", header=T, sep = ",")
head(metadatah)

# Import the counts data
files = file.path(dir, metadatah$SampleID, "quant.sf")

# check if the path to the files is valid and the files exist
all(file.exists(files))

# Generate the column names
b <- vector(mode = "list", length = nrow(metadatah))
x=1
for (i in metadatah$SampleID){
  b[[x]]<-i
  x=x+1
}
names(files) <- paste0("", b)

# Import file with annotations available as a GTF file
gtfFile <- file.path(dir, "Homo_sapiens.GRCh38.95.gtf.gz")

# Check if the path to the gffFile is valid
file.exists(gtfFile)

# Make a TxDb object from annotations available as a GTF file 
txdb <- makeTxDbFromGFF(file=gtfFile, 
                        dataSource="Homo sapiens gtf file obtained from ensembl",
                        organism="Homo sapiens")

# Using tx2gene link transcript level information to gene level counts
k <- keys(txdb, keytype = "TXNAME","GENEID")
tx2gene <- select(txdb, keys = k,  columns = "GENEID", keytype = "TXNAME") 
head(tx2gene, 10)

# Import transcript-level abundances and counts for transcript- and gene-level analysis 
txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance = "no",
                       ignoreAfterBar = TRUE, ignoreTxVersion= TRUE)

countdata <- txi.salmon$counts
head(countdata)

## Create DESeq object
dds <- DESeqDataSetFromTximport(txi = txi.salmon, 
                                colData = metadatah, 
                                design = ~Condition)
nrow(dds) 


#Remove Genes with low counts
dds1 <- dds[rowSums(counts(dds)) > 10, ]
nrow(dds1)

## Run DE analysis
dds2 <- DESeq(dds1)
colData(dds2)

#Building out results table
res_table <- results(dds2)
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
res2 <- results(dds2, alpha=0.05)
summary(res2)

#How many adjusted p-values were less than 0.05?
sum(res2$padj < 0.05, na.rm=TRUE)


#We order our results table by the smallest p value:
res_small_p <- res2[order(res2$pvalue),]
sum(res_small_p$padj < 0.05, na.rm=TRUE)

# Select genes with p less than 0.05
res_sig <- subset(res_small_p, padj < 0.05)
dim(res_sig)

#Write final list to file
write.table(as.data.frame(res_sig),"salmon_DESeq2.csv", row.names = T,sep=",",col.names=T)





