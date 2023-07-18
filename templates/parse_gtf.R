#!/usr/bin/env Rscript

# Import GTF using rtracklayer
full_gtf <- rtracklayer::import("${gtf}")
# Create a dataframe with full GTF info
gtf = as.data.frame(full_gtf)

# Save object
saveRDS(gtf, file = "${gtf.simpleName}.rds")