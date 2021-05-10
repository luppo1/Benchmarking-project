# Benchmarking-project
Here, you can find all the bash scripts and R codes used in our paper on hematological malignancies

---
title: "RNA-seq analyses pipeline for cancer research"
author: "**Albert Doughan**"
#date: "May 05, 2021"
output: github_document
---


# **Introduction**

This tutorial describes in details how to perform differential gene expression analysis in
RNA-seq studies. Here we only focused on **DESeq2** using count data from **HISAT2** aligner. 

#### Data description
Our dataset was generated from patients with Acute lymphoid leukemia, Chronic lymphocytic leukemia, Acute Myeloid leukemia and Burkitt lymphoma. Our control groups were obtained from healthy non-cancer participants of the 1000 genomes project. in all, we had 12 cases and 12 controls.

#### Quality Control, Trimming, Alignment and quantification
As with all NGS data analyses, we checked the quality of the rwa RNA-seq files with FastQC/MultiQC. Low quality bases, adapter sequences and short reads were trimmed with Trimmomatic. The reads were then aligned to the human reference genome with HISAT2. Finally, the the rate of alignment was quantified with featureCounts.

The tools above can be installed by following their respective documentations and manuals below 

1. FastQC: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> 
2. MultiQC: <https://multiqc.info/> 
3. Trimmomatic: <http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf> 
4. HISAT2: <http://www.ccb.jhu.edu/software/hisat/manual.shtml> 
5. featureCounts: <https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts> 
6. Samtools: <http://www.htslib.org/doc/samtools.html> 

**NOTE**: In all the scripts below, replace the `/path/to/` with the actual path to your respective files.
Below are the bash script used to perform the above steps:

#### Quality control
```{ }
fastqc -t 2 -o /path/to/output/folder *.fastq.gz
```
The `-t 2` specifies the number processors to be used. If unsure, remove it.

The output from FastQC can be summarized with MultiQC:
```{}
multiqc .
```

#### Trimming
```{}
#!/bin/bash
output=/path/to/output/folder
input=/path/to/input/folder
for i in $input/*_1.fastq.gz;
do
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*.fastq.gz}"
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'` 
trimmomatic PE -threads 2 -trimlog $output/"${base}".log $input/"${base}"_1.fastq.gz $input/"${base}"_2.fastq.gz $output/"${base}"_1.trimmed_PE.fastq.gz $output/"${base}"_1.trimmed_SE.fastq.gz $output/"${base}"_2.trimmed_PE.fastq.gz $output/"${base}"_2.trimmed_SE.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
done
```

The script above is for pair-end RNA-seq reads. For more information on the meaning of parameters and script for single-end reads, visit the manual/documentation.

#### Recheck the quality of trimmed reads to ensure trimming was effective
```{}
fastqc -t 2 -o /path/to/output/folder *PE.fastq.gz
```

#### Summarize output with MultiQC
```{}
multiqc .
```

#### Alignment of reads to the human reference genome
Prior to the actual alignment stage, we need to create genome index files from the reference genome. This can be done with the code below:
```{}
hisat2-build -p 2 /path/to/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa hisatindex
```
The reference genome file can be downloaded from <>

**Note**: Depending on the memory capacity of your PC or server, this step can take a while to complete.

#### Actual alignment stage
Here, we've combine two analysis steps. The first part aligns our trimmed data to the human reference genome using HISAT2. The second portion converts the output SAM format to BAM, and finally, the BAM files are sorted by coordinates. 
```{}
#!/bin/bash
input=/path/to/input/folder
output=/path/to/output/folder
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}" filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
hisat2 -p 2 -x /path/to/hisatindex -1 $input/"${base}"*_1.trimmed_PE.fastq.gz -2 $input/"${base}"*_2.trimmed_PE.fastq.gz -S $output/"${base}".hisat.sam --summary-file $output/"${base}".txt 
echo "$sample_name done!"
samtools view -@ 2 -m 2G -ub $output/"${base}".hisat.sam -o $output/"${base}".hisat.bam
echo "${sample_name} hisat.sam change to bam done!"
samtools sort -n -@ 2 -m 2G -T /tmp/ $output/"${base}".hisat.bam -o $output/"${base}".hisat.sorted.bam
echo "${sample_name} hisat.sorted.bam sort done!"
echo "$base done!"
done
```
Explanation of the parameters can be found in the tool's documentation

#### Quantification step
This step counts the number of times a read maps to a genetic feature in the human reference genome. The output is a countdata file which will be used for differential expression analysis in R in the next stage.
```{}
featureCounts -p -T 2 -t exon -g gene_id --extraAttributes gene_id,gene_biotype -a /path/to/Homo_sapiens.GRCh38.92.gtf -o hisat2.txt *.sorted.bam
```
<br />
## Differential expression analysis in R
Aside the installation of all required packages in R, two files will be needed:

1. Coundata: this is obtained from the alignment stage
2. Meta data file: This is a user-defined CSV files with two columns (SampleID and Condition).   

The countdata and meta data files used in this tutorial can be downloaded [here](https://drive.google.com/drive/folders/183mQmsTIQ0xyClEFVSunscr3MZmsJU-z?usp=sharing)  

### Setting up R and RStudio environment
If you do not already have [R statistical tool](https://cloud.r-project.org/bin/windows/base/R-4.0.5-win.exe) and [RStudio](https://download1.rstudio.org/desktop/windows/RStudio-1.4.1106.exe) installed, these can be downloaded from their respective websites. 

### Installing all required packages
Before testing for differential expression, we need to install some R packages. The packages are maintained by [Bioconductor](https://bioconductor.org/). To install a package, copy and paste the name of the package in the search box on the top right corner of the Bioconductor page, select the first result and you should be presented with a page with a code like the one below:


```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

Repeat the step above by replacing the name of the package *"tximport"* with the next package name. Do this for all packages.

### Setting up working directory
Next, we set our working directory, which will enable us to circumvent specifying long paths to our input files. We do this by following the steps:

1. Create a folder on your desktop and name it _RNA-seq_.
2. Move the countdata and metadata files you downloaded into the _RNA-seq_ folder.
3. The path to the _RNA-seq_ folder will be something like *"C:/Users/Albert Doughan/Desktop/RNA-seq"*

Remember to replace "Albert Doughan" with your username.
Run the following lines by pressing **Ctrl + Enter**

```{r, eval=TRUE}
setwd("C:/Users/Albert Doughan/Desktop/RNA-seq")
dir = "C:/Users/Albert Doughan/Desktop/RNA-seq"
```

### Loading all required packages into R
We will load all the packages we installed above into R with the following:

```{r message=FALSE, warning = FALSE}
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library("ggplot2")
library("affy")
library("ggfortify")
```


The codes below will import the metadata and countdata into R

```{r}
metadatah = read.csv(file= "metadata.csv", header=TRUE, sep = ",")
head(metadatah)
counts <- read.csv("hisat2.csv", header=TRUE, sep = ",")
```

### Data preprocessing
The countdata file contains some columns which are irrelevant to our analysis and will be removed by the following code:
```{r}
countdata = counts[, c(7:32)]
countdata = countdata[, -c(2)]
countdata <- countdata[, -c(1)]
rownames(countdata) <- counts[,1] # Make GeneID the row names
```

Display all the sample names
```{r}
colnames(countdata) 
```

Inspect the first 6 rows of the dataset
```{r}
head(countdata)
```


We can check if the metadata and samples have the same names. This will return **TRUE** if the names in both files are the same. In all there are **24** samples in countdata.
```{r}
table(colnames(countdata)==metadatah$SampleID)
```

### Running differential expression analysis with DESeq2
#### Create the DESeqDataSet object from the countdata and metadata

```{r message=FALSE, warning = FALSE}
dds <- DESeqDataSetFromMatrix(countData = round(countdata),
                              colData = metadatah,
                              design = ~Condition)
```

Check the number of genes in the dataset
```{r}
nrow(dds)
```


There are 58395 genes in our samples

#### Remove Genes with low counts
Removing rows with low count will reduce the memory size of the `dds` object and increase the speed of subsequent steps. Moreover, genes with very low counts will not be interesting biologically. The code below will keeps genes with at least 10 counts.

```{r}
dds1 <- dds[ rowSums(counts(dds)) >= 10, ]
nrow(dds1)
```
Removing low count genes reduces the number from 58395 to 34456. 

#### Run the DESeq function to perform differential expressiona analysis
```{r message=FALSE, warning = FALSE}
dds1 <- DESeq(dds1)
```

Inspect the first 6 rows
```{r}
head(assay(dds1)) 
```

#### Building the results table
This is a dataframe that contains our differentially expressed genes, their p-values, etc. 
`summary` gives us the genes that are up and down regulated in our condition under study, as well as low count genes and outliers. 
```{r}
res_table <- results(dds1)
summary(res_table)
```
Here, we found no outliers, and had 9.8% low counts genes. 

#### Working with alpha 0.05
We only consider genes with p values less than 0.05

```{r}
res2 <- results(dds1, alpha=0.05)
dim(res2)
```

How many adjusted p-values were less than 0.05?
```{r}
sum(res2$padj < 0.05, na.rm=TRUE)
```

Select genes with p value less than 0.05
```{r}
res_sig <- subset(res2, padj < 0.05)
dim(res_sig)
```
Out the 34456 total number of genes, 16141 had p-values less than 0.05

We then order our results table by the smallest p value:
```{r}
res_small_p <- res_sig[order(res_sig$pvalue),]
dim(res_small_p)
```

Finally, we Write the differentially expressed gene list to file
```{r}
write.csv(as.data.frame(res_small_p), "DESeq2_hisat.csv")
```

<br />
<br />

### Data visualization and transformation
So far, we've not been able to visually appreciate our countdata. In the ensuing steps, we will explore our data by plotting PCA, density plots and heatmaps.

#### Principal component analysis

```{r}
PCAdata <- prcomp(t(assay(dds1)))
autoplot(PCAdata, data = metadatah,colour = "Condition", label = FALSE)
```

Here, we observe that there is no clear clusters formed as expected for our cases and control groups. Let's try a hierarchical clustering to see if the problem persists.

#### Hierarchical clustering
```{r}
clusters2 <- hclust(dist( t( assay(dds1) ) ), method ="ward.D")
plot(clusters2, labels = FALSE)
```

The tree above also has two main branches, which a better, however, the sub-branches do not match our number of cases and controls. Let's try one last exploratory plot.

#### Density plot
```{r}
plotDensity(assay(dds1), col=1:24,lwd=2,lty=1,xlab("Density"),ylab("Counts"))
```

Yeah. It seems there is something wrong with our data. The good news is we can fix it through **normalization**!!.

#### Normalization
There are several normalization methods that can be applied to our data. However, due the large number of samples, we recommend Variance Stabilizing transformation (VST) method as this performs quite well on large samples. 

```{r}
vst = vst(dds1, blind=FALSE)
```
`blind=FALSE` greatly reduces the run time.

Now let's generate the plots again and assess if normalization made any difference.

#### Principal component analysis after normalization

```{r}
PCAdata <- prcomp(t(assay(vst)))
autoplot(PCAdata, data = metadatah,colour = "Condition",label = FALSE, main="PCA")
```

It seems VST normalization has greatly improved our PCA plot. The RNA-seq data used in this project were obtained from 3 cancer types, which have been clustered as such (in red). The normal control samples are also clustered nicely.

#### Hierarchical clustering after normalization
```{r}
clusters2 <- hclust(dist( t( assay(vst) ) ),method ="ward.D")
plot(clusters2, main = "Dendrogram", label = metadatah$Condition)
```

Wow!! This also shows an excellent cluster for both cases and control groups.

#### Density plot after normalization
```{r}
plotDensity(assay(vst), lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot")
```

This is better than the previous and is typical of RNA-seq data. However, the unevennes of the plot may suggets the presence of batch effects. More details on this is available in [DESeq2's manual](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results)

<br />

#### Session information
```{r}
sessionInfo()
```

For detailed explanation on all the steps, kindly read the [DESeq2 paper](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results)

<br />


**CONGRATULATIONS!!!**. You are now a master of RNA-seq analysis with HISAT2 and DESeq2.

<br />

##### In case you encounter any challenges, feel free to contact me via
[<img align="left" alt="adoughan | Twitter" width="50px" src="C:/Users/Albert Doughan/Downloads/twitter.svg" />][twitter]  
<br />
[<img align="left" alt="adoughan | LinkedIn" width="50px" src="C:/Users/Albert Doughan/Downloads/linkedIn.svg" />][linkedin]  
<br />
[<img align="left" alt="adoughan | Instagram" width="50px" src="C:/Users/Albert Doughan/Downloads/facebook.svg" />][facebook]



[twitter]: https://twitter.com/adoughan1
[linkedin]: https://www.linkedin.com/in/albert-doughan-4a7564aa/
[facebook]: https://web.facebook.com/alberto.doughani



