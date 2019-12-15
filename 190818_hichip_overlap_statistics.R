# Last Updated: 
# Author: Wesley Cai
# Purpose: Can do this in multiple ways: 1. Jaccard statistic, 2. regioneR

library(data.table)
library(ggplot2)
library(regioneR)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- ""
inputfolder <- "input/"
outputfolder <- "output/"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

mda.atac.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_atac.splitbysummit.bed"
mda.h3k27ac.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K27ac_merge_peaks.broadPeak"
mda.h3k4me1.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me1_merge_peaks.broadPeak"
mda.h3k4me3.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me3_merge_peaks.broadPeak"

h.h3k27ac.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/H2030_Par_BM3_H3K27ac.broad_peaks.bed"
hichip.lin.rep1.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/peak/rep1.cut.nochr.filterchrchr.macs_peaks.broadPeak"
hichip.lin.rep2.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/peak/rep2.cut.nochr.filterchrchr.macs_peaks.broadPeak"

hichip.lin.rep1 <- fread(hichip.lin.rep1.file)
plot(density(hichip.lin.rep1$V9), xlim = c(0,10))
quantile(hichip.lin.rep1$V9, prob = c(0,0.75,1))

hichip.lin.rep2 <- fread(hichip.lin.rep2.file)
plot(density(hichip.lin.rep2$V9), xlim = c(0,10))
quantile(hichip.lin.rep2$V9, prob = c(0,0.75,1))

# Can visually compare chr1:63,941,061-64,838,890 between MDAs and H2030s and HiChIP
# Better yet, can use SPARC chr5:150,946,545-151,272,383
# Also can use GATA3 chr10:8,063,865-8,185,917

mda.h3k27ac <- fread(mda.h3k27ac.file)
mda.h3k27ac[, V1 := paste0("chr", V1)]
unique(mda.h3k27ac$V1)
h.h3k27ac <- fread(h.h3k27ac.file)
h.h3k27ac[, V1 := paste0("chr", V1)]
unique(h.h3k27ac$V1)


hichip.lin <- hichip.lin[V9 > 10,.(V1, V2, V3)]
hichip.lin[, V1 := paste0("chr", V1)]
unique(hichip.lin$V1)





hichip.lin.gr <- toGRanges(hichip.lin)
mda.h3k27ac.gr <- toGRanges(mda.h3k27ac)
h.h3k27ac.gr <- toGRanges(h.h3k27ac)
pt <- overlapPermTest(hichip.lin.gr, mda.h3k27ac.gr, alternative = "auto", ntimes = 100)

pt.h <- overlapPermTest(hichip.lin.gr, h.h3k27ac.gr, alternative = "auto", ntimes = 100)

pt.m.h <- overlapPermTest(mda.h3k27ac.gr, h.h3k27ac.gr, alternative = "auto", ntimes = 100)

# genome <- filterChromosomes(getGenome("hg19"), keep.chr="chr1")
# A <- createRandomRegions(nregions=20, length.mean=10000000, length.sd=20000, genome=genome, non.overlapping=FALSE) 
# B <- c(A, createRandomRegions(nregions=10, length.mean=10000, length.sd=20000, genome=genome, non.overlapping=FALSE))
# 
# pt <- overlapPermTest(A=A, B=B, ntimes=10, genome=genome, non.overlapping=FALSE, verbose=TRUE)
# summary(pt)
# plot(pt)
# plot(pt, plotType="Tailed")  
