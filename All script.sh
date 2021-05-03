## Inital Quality checks
fastqc -t 27 -o ~/albert/benchmark/QC/before *.fastq.gz

## Trimming with Trimmomatic
#!/bin/bash
output=~/albert/benchmark/trimmed
input=/home/knuthstore/albert/AML
for i in $input/*_1.fastq.gz;
do
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*.fastq.gz}"
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'` 
trimmomatic PE -threads 30 -trimlog $output/"${base}".log.gz $input/"${base}"_1.fastq.gz $input/"${base}"_2.fastq.gz $output/"${base}"_1.trimmed_PE.fastq.gz $output/"${base}"_1.trimmed_SE.fastq.gz $output/"${base}"_2.trimmed_PE.fastq.gz $output/"${base}"_2.trimmed_SE.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
done


## Quality checks after trimming
fastqc -t 28 -o /home/abdul/albert/benchmark/QC/after *PE.fastq.gz


## Alignment with Kallisto
## Indexing
kallisto index -i /home/abdul/albert/Reference/kallistoindex.idx -k 31 Homo_sapiens.GRCh38.cdna.all.fa

#Mapping
#!/bin/bash 
input=~/albert/benchmark/trimmed
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}"
filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}" 
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
kallisto quant -i /home/abdul/albert/Reference/kallistoindex.idx -o ~/albert/benchmark/mapping/kallisto/"${base}" --genomebam --bias --gtf /home/abdul/albert/Reference/Homo_sapiens.GRCh38.99.gtf.gz --chromosomes /home/abdul/albert/Reference/chrom.txt -b 50 -t 27 $input/"${base}"_1.trimmed_PE.fastq.gz $input/"${base}"_2.trimmed_PE.fastq.gz &> ~/albert/benchmark/mapping/kallisto/"${base}".log
done



## Alignment with Salmon
## Builiding Salmon Index
grep "^>" <(gunzip -c /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat /home/abdul/albert/Reference/Homo_sapiens.GRCh38.cdna.all.fa.gz /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > gentrome.fa.gz
~/albert/benchmark/mapping/salmon/salmon-latest_linux_x86_64/bin/./salmon index -t /home/abdul/albert/Reference/gentrome.fa.gz -i salmon_index --decoys /home/abdul/albert/Reference/decoys.txt --keepDuplicates -p 1 -k 29
salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31


## Mapping
#!/bin/bash
input=~/albert/benchmark/trimmed
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}" filename=${withpath##*/}
base="${filename%*_*1.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
~/albert/benchmark/mapping/salmon/salmon-latest_linux_x86_64/bin/./salmon quant -i ~/albert/benchmark/mapping/salmon/salmon_index -l A -1 $input/"${base}"*_1.trimmed_PE.fastq.gz -2 $input/"${base}"*_2.trimmed_PE.fastq.gz --validateMappings -p 27 -o ~/albert/benchmark/mapping/salmon/"${base}" –-rangeFactorizationBins 4 --seqBias --numBiasSamples 2000000 --gcBias --numBootstraps 50
done



## Alignment with HISAT2
## Building Index
hisat2-build -p 27 /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa hisatindex

## Mapping
#!/bin/bash
input=~/albert/benchmark/trimmed
output=~/albert/benchmark/mapping/hisat2
for i in $input/*_1.trimmed_PE.fastq.gz;
do 
withpath="${i}" filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}"
sample_name=`echo "${base}" | awk -F "1.trimmed_PE.fastq.gz" '{print $1}'`
hisat2 -p 27 -x ~/albert/benchmark/mapping/hisat2/hisatindex -1 $input/"${base}"*_1.trimmed_PE.fastq.gz -2 $input/"${base}"*_2.trimmed_PE.fastq.gz -S $output/"${base}".hisat.sam --summary-file $output/"${base}".txt 
echo "$sample_name done!"
samtools view -@ 26 -m 26G -ub $output/"${base}".hisat.sam -o $output/"${base}".hisat.bam
echo "${sample_name} hisat.sam change to bam done!"
samtools sort -n -@ 16 -m 2G -T /tmp/ $output/"${base}".hisat.bam -o $output/"${base}".hisat.sorted.bam
rm $output/"${base}".hisat.sam $output/"${base}".hisat.bam
echo "${sample_name} hisat.sorted.bam sort done!"
echo "$base done!"
done

## Quantification step
featureCounts -p -T 27 -t exon -g gene_id --extraAttributes gene_id,gene_biotype -a ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf -o hisat2.txt *.sorted.bam


## Alignment with STAR
## Indexing
#!/bin/bash 
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir ~/albert/benchmark/mapping/star --genomeFastaFiles /home/abdul/reference_genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf --genomeSAindexNbases 12 --genomeSAsparseD 2 --sjdbOverhang 100


## Mapping
#!/bin/bash 
input=~/albert/benchmark/trimmed
for i in $input/*_1.trimmed_PE.fastq.gz; 
do 
withpath="${i}"
filename=${withpath##*/} 
base="${filename%*_*1.trimmed_PE.fastq.gz}" 
sample_name=`echo "${base}" | awk -F ".fastq.gz" '{print $1}'`
mkdir ~/albert/benchmark/mapping/star/"${base}"
STAR --runThreadN 12 --genomeDir ~/albert/benchmark/mapping/star/index --readFilesIn $input/"${base}"*1.trimmed_PE.fastq.gz $input/"${base}"*2.trimmed_PE.fastq.gz --outFileNamePrefix ~/albert/benchmark/mapping/star/"${base}"/"${base}"_ --quantMode GeneCounts --sjdbOverhang 100 --genomeSAsparseD 2 --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted SortedByCoordinate --outSAMmapqUnique 60
done


## Quantification
featureCounts -p -T 27 -t exon -g gene_id --extraAttributes gene_id,gene_biotype -a ~/reference_genomes/ensembl/Homo_sapiens.GRCh38.92.gtf -o star.txt *sortedByCoord.out.bam
