# Last Updated: 
# Author: Wesley Cai
# Purpose: Get expressed genes and also TSS ALL TSS's in biomart... multiple per gene

library(data.table)
library(ggplot2)
library(HTSFilter)
library(DESeq2)
library(biomaRt)

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/expressed_genes"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#### Get filtered genes and STAR DEG ####
allgenes <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/genomes/genes/output/ensembl_hgnc_map.txt")

dir <- "/Users/wcai/Google_Drive/_Lab/Data/RNA-seq/E190311_MDA-MB-231_RNA-seq/analyses/STAR/counts/"
names <- dir(dir)
files <- file.path(dir, names)
names(files) <- tolower(gsub("ReadsPerGene.out.tab", "", names))
# Remove BoM and BrM3
files <- files[grep("brm2|lm2|par", names(files))]

# Read in counts and generate dds object
tmp <- fread(files[1], skip = 4)
counts <- data.frame(row.names = tmp$V1)
for(f in 1:length(files)){
  tmp <- fread(files[f], skip = 4)
  counts[,names(files)[f]] <- tmp$V2
}
colData <- data.frame(row.names = colnames(counts), line = gsub("-.*", "", colnames(counts)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ line)
dds <- DESeq(dds)

# Generate filter for low genes https://academic.oup.com/bioinformatics/article/29/17/2146/240530
filter <- HTSFilter(dds, s.len=25, plot=FALSE)$filteredData

lines <- list()
for(line in c("brm2", "lm2")){
  tmp <- results(filter, independentFiltering=FALSE, contrast = c("line", line, "par"), tidy = TRUE)
  colnames(tmp)[3] <- "l2fc"
  lines[[line]] <- tmp[,c(3,7)]
}
final <- data.table(ensembl = tmp$row)
final2 <- as.data.table(do.call(cbind, lines))
final <- cbind(final, final2)
final <- merge(final, allgenes, by.x = "ensembl", by.y = "ensembl_gene_id", all.x = TRUE)
ggplot(final[brm2.padj < 0.05 & lm2.padj < 0.05,], aes(brm2.l2fc, lm2.l2fc)) + 
  geom_point()
#### Get filtered genes and STAR DEG ####

#### Get TSS of genes ####
# Remember to get hg19!!!!!!
bm <- useMart("ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
test <- listAttributes(bm)
bm.res <- getBM(c("ensembl_gene_id", "chromosome_name", "start_position", "transcription_start_site", "strand"), filters = "ensembl_gene_id", values = final$ensembl, bm)
bm.dt <- as.data.table(bm.res)
bm.dt <- bm.dt[,.(ensembl_gene_id, chr = chromosome_name, tss = transcription_start_site, strand = strand)]
# bm.dt <- bm.dt[,.(tss = transcription_start_site[transcription_start_site*strand == min(transcription_start_site*strand)],
#                   strand = strand[transcription_start_site*strand == min(transcription_start_site*strand)]), by = .(ensembl_gene_id)]
bm.dt[strand == 1, start:= as.numeric(tss)]
bm.dt[strand == 1, end:= as.numeric(tss+1)]
bm.dt[strand == -1, start:= as.numeric(tss-1)]
bm.dt[strand == -1, end:= as.numeric(tss)]
final.tss <- merge(final, bm.dt[,.(ensembl = ensembl_gene_id, chr, start, end, strand)], by = "ensembl", all.x = TRUE)
fwrite(final.tss, "par_brm2_lm2.expressed.txt", sep = "\t")
#### Get TSS of genes ####
