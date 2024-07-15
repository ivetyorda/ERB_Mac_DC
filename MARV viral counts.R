#Explanation of DESeq2
DESeq2 performs the following steps for analysis:
Import and organize count data: This is where the counts for each gene (or other feature) in each sample are loaded into a DESeq2 object.
Normalization: DESeq2 performs a normalization step to correct for differences in sequencing depth between different samples.
Dispersion estimation: DESeq2 estimates the dispersion of counts for each gene, accounting for both biological and technical variability.
Statistical testing: DESeq2 uses the negative binomial distribution to model count data and performs a hypothesis test to identify differentially expressed genes.
The Negative Binomial distribution is a discrete probability distribution that models the number of successes in a sequence of independent and identically distributed Bernoulli trials before a specified (non-random) number of failures occurs.
Multiple testing correction: Finally, DESeq2 applies a correction to control the False Discovery Rate (FDR), which is a statistical method to correct for multiple hypothesis testing.
The result of DESeq2 analysis is a list of genes (or other features), along with log2 fold changes and adjusted p-values indicating the significance of differential expression between conditions.
Note: DESeq2 should only be used with raw counts data, as it uses the raw counts to estimate dispersions and perform statistical testing. It should not be used with data that has already been normalized or transformed, such as by converting counts to Transcripts Per Million (TPM) or Reads Per Kilobase Million (RPKM/FPKM).

## Open packages
library(tximport)
library(readr)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
BiocManager::install("apeglm")
library(apeglm)
library("BiocParallel")
register(SnowParam(4))
library("pheatmap")
library("tidyverse")
library(edgeR)
library(RColorBrewer)

## Set directory to location of TSV values----
setwd("C:/Users/user/folder")
files = list.files()

## Reorder files, eg. when 4 files total
files = files[c(1,2,3,4)]
list.files()
all(file.exists(files))

## Use these to set sample names in order 
samples=c("MARV1_1", "MARV1_2", "MARV1_3", "MARV1_4")

## In advance, download the GTF genome file into your working folder from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/857/325/GCF_000857325.2_ViralProj15199/
## Load the GTF/GFF file gtf_file
library(rtracklayer)
gtf_file <- "C:/Users/user/folder/GCF_000857325.2_ViralProj15199_genomic.gtf"
gtf <- import(gtf_file)

## Extract transcript-to-gene mappings

#tx2gene <- gtf[gtf$type == "transcript"]: This line of code creates a new data frame, tx2gene, which is a subset of the gtf data frame, containing only the rows where the type column is equal to "transcript". In the context of a GTF file, a "transcript" generally refers to the mRNA sequence that is transcribed from a gene.

#tx2gene <- data.frame(transcript_id = tx2gene$transcript_id, gene_id = tx2gene$gene_id): Here, tx2gene is being overwritten as a new data frame with only two columns - transcript_id and gene_id from the original tx2gene data frame. transcript_id corresponds to the unique identifiers of the transcripts, and gene_id corresponds to the unique identifiers of the genes.

tx2gene <- gtf[gtf$type == "transcript"]
tx2gene <- data.frame(transcript_id = tx2gene$transcript_id, gene_id = tx2gene$gene_id)

## Import transcript abundance data # removed "ignoreTxVersion = TRUE" which led to an error in annotation, in this case, there are multiple versions for each transcript.
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

## Set up metadata. Example namings provided for 3 replicates per bat, for 2 bats.
cell_type = as.factor(c(rep("bat1",3), rep("bat2",3)))

#A data frame is created, meta_data, with two columns: sample and cell_type. It is assumed that samples and cell_type are already defined elsewhere in your code.
meta_data = data.frame(sample = samples, cell_type = cell_type)

#general dds----
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta_data,
                                design = ~ cell_type)
as.data.frame(colData(dds))

dds<- DESeq(dds)
class(dds)

rld<-rlog(dds)
write.csv(as.data.frame(assay(rld)), file="filename.csv" )

