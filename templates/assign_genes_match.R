#!/usr/bin/env Rscript

library(dplyr)

# Import GTF using rtracklayer
full_gtf <- readRDS("${gtf_rds}")
gtf <- full_gtf %>% filter(type == "gene") %>% select(seqnames, start, end, strand, gene_name)
colnames(gtf) <- c("chromosome", "start", "end", "strand", "gene_name")

bed <- read.table("${bed}", header=F, col.names=c("chromosome", "start", "end", "sample_info", "MAPQ", "strand"))

bed <- bed %>% tidyr::separate(sample_info, c("num", "sgName", "gene_name"), "\\\\|")

by <- join_by(chromosome, strand, within(x\$start, x\$end, y\$start, y\$end))
results <- right_join(bed, gtf, by)

gene_name_match = results %>% filter(gene_name.x == gene_name.y) %>% pull(gene_name.x) %>% unique()

saveRDS(gene_name_match, file = "gene_name_match.rds")
