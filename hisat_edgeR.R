# Load all required packages
library("tximport")
library("pheatmap")
library("RColorBrewer")
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

#Creating the edgeR gene list
d <- DGEList(counts=countdata,group=factor(metadatah$Condition))
dim(d)

#Filter out low count genes
head(d$counts)
head(cpm(d))

#Total gene counts per sample
apply(d$counts, 2, sum) 

#Remove low CPM genes
#The filterByExpr function keeps rows that have worthwhile counts in a minimum number
#of samples. The function #accesses the group factor contained in d in order to compute the minimum group size, but
#the filtering is performed independently of which sample belongs to which group so that no
#bias is introduced.
keep <- rowSums(cpm(d)>10) >= 2
d <- d[keep,]
dim(d)

#After filtering, it is a good idea to reset the library sizes:
d$samples$lib.size <- colSums(d$counts)
d$samples

#Normalizing the data
# The calcNormFactors() function normalizes for RNA composition by finding a set 
#of scaling factors for the library sizes that minimize the log-fold changes between 
#the samples for most genes. The default method for computing these scale factors 
#uses a trimmed mean of M-values (TMM) between each pair of samples.
d <- calcNormFactors(d)
d

# Exploratory analysis
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

#Estimating the Dispersion
d1 <- estimateCommonDisp(d, verbose=T) # estimating common dispersions
names(d1)

d1 <- estimateTagwiseDisp(d1) # estimating tagwise dispersions
names(d1)

#plotBCV() plots the tagwise biological coefficient of variation (square root of 
#dispersions) against log2-CPM.
plotBCV(d1)

##Testing for DE genes
# Classical Exact test approach
et <- exactTest(d1)
topTags(et, n=10)

#The total number of differentially expressed genes at FDR< 0:05 is:
de1 <- decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(de1)

#The function plotSmear generates a plot of the tagwise log-fold-changes against log-cpm 
#(analogous to an MA-plot for microarray data). DE tags are highlighted on the plot:
# differentially expressed tags from the naive method in d1
de1tags <- rownames(d1)[as.logical(de1)] 
plotSmear(et, de.tags=de1tags)
abline(h = c(-1, 1), col = "blue")  

#Order by P-value  
res = et$table
order_res <- res[order(res$PValue),]
dim(order_res)

# Select only genes with P-value less that 0.05
sig_pvalue <- subset(order_res, PValue < 0.05)
dim(sig_pvalue)

#Write final list to file
write.csv(as.data.frame(sig_pvalue), "edgeR_hisat2.csv")
