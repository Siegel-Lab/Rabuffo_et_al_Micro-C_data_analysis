#!/usr/bin/env Rscript

# script to bin mappability bed file for mHiC ICEing
# Author: Benedikt G Brink, LMU Munich, 2019
## Changed by Claudia Rabuffo, LMU Munich for use without Snakemake


library(GenomicRanges)
args <- commandArgs(trailingOnly=TRUE)

args_fai_input = args[1]
args_bed_input = args[2]
args_binsize = args[3]
args_output = args[4]

print(args_fai_input)
print(args_bed_input)
print(args_binsize)
print(args_output)


fai <- read.table(args_fai_input, stringsAsFactors = FALSE)
bedfile <- read.table(args_bed_input, stringsAsFactors = TRUE)
binsize <- as.numeric(args_binsize)

bed_ranges <- GRanges()
ice_frame <- data.frame(chr = vector(), midPoint = vector(), anything = vector(), mappabilityValue = vector())

for (chr in levels(bedfile$V1)) {
  qwe <- bedfile[which(bedfile$V1 == chr),]
  asd <- GRanges(seqnames = qwe$V1, IRanges(start = qwe$V10, end = qwe$V11), mappability = qwe$V12)
  chrlength <- fai[which(fai$V1 == chr),2]
  end(asd[length(asd)]) <- chrlength
  bed_ranges <- c(bed_ranges, asd)
}

bins <- slidingWindows(range(bed_ranges), binsize, binsize)
bins <- unlist(bins)
hits <- findOverlaps(bins, bed_ranges)
agg <- aggregate(bed_ranges, hits, score=mean(mappability))
bins$mappability <- agg$score

ice_frame <- data.frame(chr=seqnames(bins), 
                        midPoint=start(bins)+binsize/2, 
                        anything=rep(NA, length(bins)),
                        mappability=bins$mappability)

write.table(ice_frame, file = args_output, sep = "\t", quote = FALSE, row.names = FALSE)

