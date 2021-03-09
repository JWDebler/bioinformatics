#!/usr/bin/env Rscript

# Author: Zaka Yuen, JCSMR, ANU
# Created on July 2020
# Adapted from https://github.com/comprna/METEORE/tree/master/script
# Edited: Johannes Debler

# This script is to:
# - take the default aggregated modified base output (bedMethyl format) after running Megalodon for modified base detection
# - combine the CpG methylation calls from both strands into a single strand for a CpG site
# - Format it in a way that GenomeBrowse (and in the future methplotlib) can import

library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  end("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "megalodon_freq-perCG.tsv"
}

df <- read.table(args[1], header = FALSE, sep = "\t")
# Remove unnecessary columns
df <- df  %>% select(1, 2, 3, 6, 11, 10)
colnames(df) <- c("chromosome", "start", "end", "strand", "methylated_frequency", "Cov")
# Convert %methylation into frequency 
df[,"methylated_frequency"] <- df[,"methylated_frequency"]/100


# Accumulate CpG sites into the +'ve strand
df[df$strand == "-","start"] <- df[df$strand == "-","start"]-1

# Use data.table to compute the mean of duplicated position while keeping non-duplicated sites
df <- data.table(df)
#str(df)
df <- df[,list(methylated_frequency = mean(methylated_frequency), 
                Cov = sum(Cov)),
                list(chromosome,start,end)]
df <- data.frame(df)
df$start <- format(df$start,scientific=FALSE)
write.table(df, file=args[2],  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
