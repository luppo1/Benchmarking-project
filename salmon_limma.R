# First load all the libraries that will be needed for the analysis.
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
write.csv(order_res,"limma_salmon.csv")
