# sgRNA pipeline

This repository is to answer a challenge from Ryvu Therapeutics.

## sgRNA library characterization

* Perform mapping of sgRNA sequences (provided multifasta file: library.fa) to GRCh38.
* For each of the mapped sgRNA sequence provide the following information: chromosome name, start and end positions, strand.
* Assign gene names based on mapping results and compare obtained annotation with gene names provided within description lines of fasta file.
* Retrieve matrix of gene expression for the gene IDs from previous step for sample "TCGA-A7-A13D-01A-13R-A12P-07" and "TCGA-E9-A1RH-11A-34R-A169-07" from TCGA-BRCA dataset. 

Please wrap your solution in a simple nextflow pipeline; use docker for virtualization of used software and provide short documentation. 
As an outcome, please provide the link to the git with code, documentation and results.

## Workflow

1) Align reads to a short-reads aligner as BWA and generate a BED file; (BWA_ALIGN; BWA_INDEX; BEDTOOLS);
2) Parse GTF information to get chr, start, end, strand, gene_name (PARSE_GTF);
3) Intersect reference GTF info witwh alignment BED + Create a CSV with annotated info (ASSIGN_GENES_MATCH);
4) Retrieve info from TCGA-BRCA based on sample-id and intersect with gene_name from task 2 (RETRIEVE_TCGA);

## Usage

```sh
## To get info just run:
nextflow run diegomscoelho/sgRNA_pipeline -r main --help

## Basic usage:
nextflow run diegomscoelho/sgRNA_pipeline -r main \
    --tcga_samples "TCGA-A7-A13D-01A-13R-A12P-07,TCGA-E9-A1RH-11A-34R-A169-07"

## This will automatically download necessary references + align/annotate reads in
## `inputs` folder and get TCGA info; To avoid TCGA step just not use: `--tcga_samples` flag

```