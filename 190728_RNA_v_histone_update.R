# Last Updated: 
# Author: Wesley Cai
# Purpose: Plot Histone vs RNA-seq

library(data.table)
library(ggplot2)

source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/plot_genes"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

rna <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/expressed_genes/par_brm2_lm2.expressed.txt")
rna <- rna[,.SD, .SDcols = -grep("chr|start|end|strand", colnames(rna))]
rna <- rna[!duplicated(rna),]

#### Get one TSS per gene based on H3K4me3 peak that changed the most between met and parental ####
h3k4me3 <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me3_merge_peaks.tss.bed")
h3k4me3[,ensembl := gsub("_.*", "", V8)]
h3k4me3 <- h3k4me3[V10 == 0,]
h3k4me3.res.tmp <- res.list[["h3k4me3.broad"]]
colnames(h3k4me3.res.tmp) <- paste0(colnames(h3k4me3.res.tmp), ".h3k4me3")
h3k4me3.merge <- merge(h3k4me3, h3k4me3.res.tmp, by.x = "V4", by.y = "peakid.h3k4me3", all.x = TRUE)
h3k4me3.rna.merge <- merge(h3k4me3.merge, rna, by = "ensembl", all.x = TRUE)
tss.final <- h3k4me3.rna.merge[,.SD[which.min(get(paste0("met_padj.h3k4me3")))], by = .(ensembl)]
tss.final.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/tss.final.bed"
tss.final.unsort.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/tss.final.unsort.bed"
fwrite(tss.final[,.(V5, V6, V7, V8, V9)], tss.final.unsort.file, sep = "\t", col.names = FALSE)
system(paste0("source ~/.bashrc; bedtools sort -i ", tss.final.unsort.file,
              " > ", tss.final.file))
#### Get one TSS per gene based on H3K4me3 peak that changed the most between met and parental ####

#### Promoter ####
h3k4me3 <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me3_merge_peaks.tss.bed")
h3k4me3[,ensembl := gsub("_.*", "", V8)]
h3k4me3 <- h3k4me3[V10 == 0,]
h3k4me3.res.tmp <- res.list[["h3k4me3.broad"]]
colnames(h3k4me3.res.tmp) <- paste0(colnames(h3k4me3.res.tmp), ".h3k4me3")
h3k4me3.merge <- merge(h3k4me3, h3k4me3.res.tmp, by.x = "V4", by.y = "peakid.h3k4me3", all.x = TRUE)
h3k4me3.rna.merge <- merge(h3k4me3.merge, rna, by = "ensembl", all.x = TRUE)

line <- "lm"
for(line in c("lm", "brm")){
  # Take 1 H3K4me3 peak within promoter region, pick one with lowest p.value
  prom <- h3k4me3.rna.merge[,.SD[which.min(get(paste0(line, "_padj.h3k4me3")))], by = .(ensembl)]
  prom[get(paste0(line, "_l2fc.h3k4me3")) < 0, h3k4me3.binary := "Decreased"]
  prom[get(paste0(line, "_l2fc.h3k4me3")) > 0, h3k4me3.binary := "Increased"]
  
  tmp <- prom[get(paste0(line, "_padj.h3k4me3")) < 0.05]
  p.value <- t.test(tmp[get("h3k4me3.binary") == "Increased", get(paste0(line, "2.l2fc"))], tmp[h3k4me3.binary == "Decreased", get(paste0(line, "2.l2fc"))])$p.value
  myplot <- ggplot(tmp, aes_string("h3k4me3.binary", paste0(line, "2.l2fc"), fill = "h3k4me3.binary")) +
    stat_boxplot(geom ='errorbar', width = 0.25) + 
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ylim(c(-2, 2)) + 
    xlab(label = "Promoter H3K4me3 direction") +
    ylab(label = "RNA log2(FC)") +
    #guides(fill=guide_legend(title="Enhancer\nDirection")) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    scale_x_discrete(labels=c("Decreased" = "Down", "Increased" = "Up"))
  myplot
  ggsave(paste0(line, "2.h3k4me3.upsig.downsig.unique.png"), myplot, height = 4, width = 3.2, units = "in")
  #### Promoter ####
}
#### Promoter ####

#### Enhancers ####
h3k27ac.broad.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K27ac_merge_peaks.bed"
tss.final.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/tss.final.bed"
h3k27ac.broad.nearestKtss.file <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K27ac_merge_peaks.nearestKtss.bed"

# Nearest 100 genes
system(paste0("source ~/.bashrc; bedtools closest -k 100 -D a -a ", h3k27ac.broad.file,
              " -b ", tss.final.file, " > ", h3k27ac.broad.nearestKtss.file))
h3k27ac.broad.nearestKtss <- fread(h3k27ac.broad.nearestKtss.file)
system(paste0("source ~/.bashrc; rm ", h3k27ac.broad.nearestKtss.file))
head(h3k27ac.broad.nearestKtss)
h3k27ac.broad.nearestKtss <- h3k27ac.broad.nearestKtss[order(abs(V10)),]
h3k27ac.broad.nearestKtss <- h3k27ac.broad.nearestKtss[order(V4),]
h3k27ac.broad.nearestKtss[, seq := seq(.N), by = .(V4)]

nearestnth <- c(1, 2, 3, 100)
nearestKtss.subset <- h3k27ac.broad.nearestKtss[seq %in% nearestnth,]
nearestKtss.subset[,ensembl := gsub("_.*", "", V8)]
nearestKtss.subset <- merge(nearestKtss.subset[,.(V4, V8, seq, ensembl)], rna, by = "ensembl", all.x = TRUE)
enhancers <- res.list[["h3k27ac.broad"]][peakid %in% chip.annotations[["h3k27ac.broad.enh.active"]],]
colnames(enhancers) <- paste0(colnames(enhancers), ".histone")

line <- "lm"
for(line in c("lm", "brm")){
  p.thr <- 5e-2
  line.padj.histone <- paste0(line, "_padj.histone")
  line.l2fc.histone <- paste0(line, "_l2fc.histone")
  line.padj.rna <- paste0(line, "2.padj")
  line.l2fc.rna <- paste0(line, "2.l2fc")
  subset.merged <- merge(nearestKtss.subset[V4 %in% chip.annotations[["h3k27ac.broad.enh.active"]],], enhancers, by.x = "V4", by.y = "peakid.histone", all.x = TRUE)
  subset.merged[get(line.padj.histone) < p.thr & get(line.l2fc.histone) > 0, Direction := "Up"]
  subset.merged[get(line.padj.histone) < p.thr & get(line.l2fc.histone) < 0, Direction := "Down"]
  subset.merged$seq <- factor(subset.merged$seq, nearestnth)
  subset.merged <- subset.merged[!is.na(line.padj.histone) & !is.na(line.l2fc.histone) & !is.na(Direction),]
  
  p.value <- list()
  for(i in unique(subset.merged$seq)){
    p.value[[i]] <- t.test(subset.merged[seq %in% i & Direction == "Up",get(line.l2fc.rna)], subset.merged[seq %in% i & Direction == "Down",get(line.l2fc.rna)])$p.value
  }
  fwrite(as.data.table(do.call(rbind, p.value), keep.rownames = "nearest"), paste0(line, ".pvalues.txt"), sep = "\t")
  
  pd = position_dodge(width = 0.75)
  myplot <- ggplot(subset.merged, aes(seq, lm2.l2fc, fill = Direction)) +
    stat_boxplot(geom ='errorbar', position=pd, width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.shape = NA, position=pd, width = 0.6) +
    ylim(c(-2, 2)) + 
    xlab(label = "Nth gene from enhancer") +
    ylab(label = "RNA log2(FC)") +
    guides(fill=guide_legend(title="Enhancer\nDirection")) +
    theme_bw()
  myplot
  ggsave(paste0(line,"2.h3k27ac.upsig.downsig.nthgene.unique.pdf"), myplot, height = 4, width = 4, units = "in")
}
#### Enhancers ####
