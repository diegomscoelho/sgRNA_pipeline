# Ryvu_challenge

This repository is to answer a challenge from Ryvu Therapeutics.

## sgRNA library characterization

* Perform mapping of sgRNA sequences (provided multifasta file: library.fa) to GRCh38.
* For each of the mapped sgRNA sequence provide the following information: chromosome name, start and end positions, strand.
* Assign gene names based on mapping results and compare obtained annotation with gene names provided within description lines of fasta file.
* Retrieve matrix of gene expression for the gene IDs from previous step for sample "TCGA-A7-A13D-01A-13R-A12P-07" and "TCGA-E9-A1RH-11A-34R-A169-07" from TCGA-BRCA dataset. 

Please wrap your solution in a simple nextflow pipeline; use docker for virtualization of used software and provide short documentation. 
As an outcome, please provide the link to the git with code, documentation and results.

## How-to

1) Align reads to a short-reads aligner as BWA / BOWTIE; (BWA_ALIGN; BOWTIE_ALIGN)
2) Create a csv with read_name, chr, start, end, strand, gene_name;
3) Retrieve info from TCGA-BRCA based on sample-id and intersect with gene_name from task 2;

## Usage

```sh

nextflow run diegomscoelho/ryvu_challenge 

```