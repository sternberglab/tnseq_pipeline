#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output <- args[1]
out_dir <- args[2]

library(ggplot2)
library(tidyr)
library(tidyverse)

output_csv <- read.csv(output, header=TRUE)
Other <- output_csv$Filtered.Reads - (output_csv$Unique.Genome.Mapping.Reads + output_csv$Undigested.Donor.Reads)
Sample<-output_csv$Sample
Unique_genome_insertion_reads<-output_csv$Unique.Genome.Mapping.Reads
Contaminating_donor_reads<-output_csv$Undigested.Donor.Reads
df_new<-data.frame(Sample, Other, Unique_genome_insertion_reads, Contaminating_donor_reads)
longdf <- gather(df_new, var, value, -Sample)

barplot_mapping<-ggplot(longdf, aes(x=Sample, y=value, fill = var)) +
    geom_bar(stat = "identity", position="fill") + xlab("Sample") + ylab("Fraction of reads") +scale_fill_manual(values = c("Other" = "Dodgerblue", "Unique_genome_insertion_reads" = "Green", "Contaminating_donor_reads" = "Firebrick1")) + theme(plot.margin=unit(c(1,1,1,1),"cm")) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(args[2],"bar_mapping.pdf"))

barplot_mapping

dev.off()