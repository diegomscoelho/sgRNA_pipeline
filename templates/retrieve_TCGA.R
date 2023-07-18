#!/usr/bin/env Rscript

library(TCGAbiolinks)
library(SummarizedExperiment)

gene_name_match <- readRDS("${gene_names}")

# You can define a list of samples to query and download providing relative TCGA barcodes.
listSamples <- strsplit("${params.tcga_samples}", ",")[[1]]

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  barcode = listSamples)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

data <- SummarizedExperiment::assay(BRCARnaseqSE, "${params.assay}")
rownames(data) <- BRCARnaseqSE@rowRanges\$gene_name

sub_BRCARNAseq <- data[rownames(data) %in% gene_name_match,]

write.csv(sub_BRCARNAseq, file = "TCGA_retrieved.csv", quote = F)
