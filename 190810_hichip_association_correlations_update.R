# Last Updated: 
# Author: Wesley Cai
# Purpose: Updated on 190910 to include all RNA

library(data.table)
library(ggplot2)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/"
inputfolder <- "input/"
outputfolder <- "output/association_correlation_update"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

load("../inter.final.dt.RData")
rna <- res.list[["rna"]]
inter.final.lin <- rbind(inter.final.dt[,.(peakid = atac1.peakid, annote = annote1, target = gene2)],
                         inter.final.dt[,.(peakid = atac2.peakid, annote = annote2, target = gene1)])

inter.final.lin <- merge(inter.final.lin, res.list[["atac"]][,.(peakid, atac.brm_l2fc = brm_l2fc,
                                                                atac.brm_padj = brm_padj, atac.lm_l2fc = lm_l2fc,
                                                                atac.lm_padj = lm_padj)],
                         by = "peakid")
inter.final.lin <- merge(inter.final.lin, rna[,.(target = ensembl, rna.brm_l2fc = brm_l2fc,
                                                 rna.brm_padj = brm_padj, rna.lm_l2fc = lm_l2fc,
                                                 rna.lm_padj = lm_padj, hgnc_symbol = hgnc)],
                         by = "target")
inter.final.lin.enh <- inter.final.lin[annote == "enhancer",]
inter.final.lin.enh.gene <- inter.final.lin.enh
#inter.final.lin.enh.gene <- inter.final.lin.enh[hgnc_symbol != "",]
inter.final.lin.enh.gene <- inter.final.lin.enh.gene[!duplicated(inter.final.lin.enh.gene[,.(peakid, target, annote)])]
inter.final.lin.enh.gene <- merge(inter.final.lin.enh.gene, chip.annotations[["atac.h3k27ac.broad.overlap"]][,.(peakid = atac.peakid, h3k27ac.peakid)],
                                  by = "peakid")
inter.final.lin.enh.gene <- merge(inter.final.lin.enh.gene, res.list[["h3k27ac.broad"]][,.(h3k27ac.peakid = peakid,
                                                                                           h3k27ac.brm_l2fc = brm_l2fc,
                                                                                           h3k27ac.brm_padj = brm_padj,
                                                                                           h3k27ac.lm_l2fc = lm_l2fc,
                                                                                           h3k27ac.lm_padj = lm_padj)], by = "h3k27ac.peakid")
inter.final.lin.enh.gene[, h3k27ac.brm_dir := cut(h3k27ac.brm_l2fc, breaks = c(min(h3k27ac.brm_l2fc), 0, max(h3k27ac.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
inter.final.lin.enh.gene[, h3k27ac.lm_dir := cut(h3k27ac.lm_l2fc, breaks = c(min(h3k27ac.lm_l2fc), 0, max(h3k27ac.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
inter.final.lin.enh.gene[, rna.brm_dir := cut(rna.brm_l2fc, breaks = c(min(rna.brm_l2fc), 0, max(rna.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
inter.final.lin.enh.gene[, rna.lm_dir := cut(rna.lm_l2fc, breaks = c(min(rna.lm_l2fc), 0, max(rna.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
inter.final.lin.enh.gene[, atac.brm_dir := cut(atac.brm_l2fc, breaks = c(min(atac.brm_l2fc), 0, max(atac.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
inter.final.lin.enh.gene[, atac.lm_dir := cut(atac.lm_l2fc, breaks = c(min(atac.lm_l2fc), 0, max(atac.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]

output <- inter.final.lin.enh.gene
#output[,peakid := sub("-", ":", peakid)]
#output[,h3k27ac.peakid := sub("-", ":", h3k27ac.peakid)]
output2 <- cbind(output[,.SD,.SDcols = -grep("l2fc|padj", colnames(output))], apply(output[,.SD,.SDcols = grep("l2fc|padj", colnames(output))], c(1,2), function(x){signif(x,3)}))
expressed.tss <- sub("_.*", "", chip.annotations[["atac.tss.promoters"]][,get("gene")])
expressed.tss <- unique(expressed.tss)
length(unique(output2$target))/length(unique(expressed.tss))
# 81.8% of expressed genes have interactions!!!
# fwrite(output2, "inter.final.lin.enh.gene.txt", sep = "\t")

tmp <- inter.final.lin.enh.gene#[h3k27ac.brm_padj < 0.05 & rna.brm_padj < 0.05,]
pval <- wilcox.test(tmp[h3k27ac.brm_dir == "down", rna.brm_l2fc],tmp[h3k27ac.brm_dir == "up", rna.brm_l2fc])
ggplot(tmp[h3k27ac.brm_padj < 0.05 & rna.brm_padj < 0.05], aes(h3k27ac.brm_l2fc, rna.brm_l2fc)) +
  geom_point()

plotAssociation <- function(tmp, x, y, fill = NA, xlab, ylab){
  tmp2 <- tmp
  atac.rna.p <- wilcox.test(tmp2[get(x) == "down", get(y)],
                            tmp2[get(x) == "up", get(y)])$p.value
  mp <- ggplot(tmp2, aes_string(x, y, fill = fill)) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
    xlab(label = xlab) + ylab(label = ylab) +
    ylim(c(-5,5))
  return(mp)
}
testAssociations <- function(df, name, histone.padjCutoff){
  #histone.padjCutoff <- 0.05
  res <- list()
  tmp <- df
  line <- "lm"
  tmp[, h3k27ac.brm_dir := cut(h3k27ac.brm_l2fc, breaks = c(min(h3k27ac.brm_l2fc), 0, max(h3k27ac.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  tmp[, h3k27ac.lm_dir := cut(h3k27ac.lm_l2fc, breaks = c(min(h3k27ac.lm_l2fc), 0, max(h3k27ac.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  tmp[, rna.brm_dir := cut(rna.brm_l2fc, breaks = c(min(rna.brm_l2fc), 0, max(rna.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  tmp[, rna.lm_dir := cut(rna.lm_l2fc, breaks = c(min(rna.lm_l2fc), 0, max(rna.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  tmp[, atac.brm_dir := cut(atac.brm_l2fc, breaks = c(min(atac.brm_l2fc), 0, max(atac.brm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  tmp[, atac.lm_dir := cut(atac.lm_l2fc, breaks = c(min(atac.lm_l2fc), 0, max(atac.lm_l2fc)), include.lowest = TRUE, labels = c("down", "up"))]
  for(line in c("lm", "brm")){
    line.fixed <- lines.dict[line.lower == line,get("line.fixed")]
    
    tmp2 <- tmp[get(paste0("atac.", line, "_padj")) < 0.05 & get(paste0("rna.", line, "_padj")) < 0.05]
    mp <- plotAssociation(tmp2, paste0("atac.", line, "_dir"), paste0("rna.", line, "_l2fc"), fill = paste0("atac.", line, "_dir"), "ATAC", "RNA log2(FC)")
    ggsave(paste0(line, ".atac.rna.", name, ".png"), mp, width = 2, height = 4)
    atac.rna.p <- t.test(tmp2[get(paste0("atac.", line, "_dir")) == "down", get(paste0("rna.", line, "_l2fc"))],
                         tmp2[get(paste0("atac.", line, "_dir")) == "up", get(paste0("rna.", line, "_l2fc"))])$p.value
    
    tmp2 <- tmp[get(paste0("h3k27ac.", line, "_padj")) < histone.padjCutoff & get(paste0("rna.", line, "_padj")) < 0.05]
    histone.rna.p <- wilcox.test(tmp2[get(paste0("h3k27ac.", line, "_dir")) == "down", get(paste0("rna.", line, "_l2fc"))],
                            tmp2[get(paste0("h3k27ac.", line, "_dir")) == "up", get(paste0("rna.", line, "_l2fc"))])$p.value
    mp <- ggplot(tmp2, aes_string(paste0("h3k27ac.", line, "_dir"), paste0("rna.", line, "_l2fc"), fill = paste0("h3k27ac.", line, "_dir"))) +
      geom_hline(yintercept = 0, lty = "dotted") +
      stat_boxplot(geom = "errorbar", width = 0.25) + 
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
      xlab(label = paste0("Enhancer H3K27ac Direction")) + ylab(label = paste0("RNA log2(FC)")) +
      ylim(c(-5,5)) +
      theme(axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 16)) +
      scale_x_discrete(labels=c("down" = "Down", "up" = "Up"))
    ggsave(paste0(line, ".h3k27ac.rna.", name, ".png"), mp, width = 2.6, height = 4)
    
    tmp2 <- tmp[get(paste0("h3k27ac.", line, "_padj")) < histone.padjCutoff]
    histone.rna.p <- wilcox.test(tmp2[get(paste0("h3k27ac.", line, "_dir")) == "down", get(paste0("rna.", line, "_l2fc"))],
                                 tmp2[get(paste0("h3k27ac.", line, "_dir")) == "up", get(paste0("rna.", line, "_l2fc"))])$p.value
    mp <- ggplot(tmp2, aes_string(paste0("h3k27ac.", line, "_dir"), paste0("rna.", line, "_l2fc"), fill = paste0("h3k27ac.", line, "_dir"))) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
      xlab(label = paste0("H3K27ac")) + ylab(label = paste0("RNA log2(FC)")) +
      ylim(c(-5,5))
    ggsave(paste0(line, ".h3k27ac.nonsigrna.", name, ".png"), mp, width = 2, height = 4)
    
    tmp2 <- tmp[get(paste0("atac.", line, "_padj")) < 0.05 & get(paste0("h3k27ac.", line, "_padj")) < histone.padjCutoff]
    atac.histone.p <- wilcox.test(tmp2[get(paste0("atac.", line, "_dir")) == "down", get(paste0("h3k27ac.", line, "_l2fc"))],
                             tmp2[get(paste0("atac.", line, "_dir")) == "up", get(paste0("h3k27ac.", line, "_l2fc"))])$p.value
    mp <- ggplot(tmp2, aes_string(paste0("atac.", line, "_dir"), paste0("h3k27ac.", line, "_l2fc"), fill = paste0("atac.", line, "_dir"))) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
      xlab(label = paste0("ATAC")) + ylab(label = paste0("H3K27ac log2(FC)")) +
      ylim(c(-5,5))
    ggsave(paste0(line, ".atac.h3k27ac.", name, ".png"), mp, width = 2, height = 4)
    
    tmp2 <- tmp[get(paste0("rna.", line, "_padj")) < 0.05 & get(paste0("h3k27ac.", line, "_padj")) < histone.padjCutoff]
    rna.histone.p <- wilcox.test(tmp2[get(paste0("rna.", line, "_dir")) == "down", get(paste0("h3k27ac.", line, "_l2fc"))],
                            tmp2[get(paste0("rna.", line, "_dir")) == "up", get(paste0("h3k27ac.", line, "_l2fc"))])$p.value
    mp <- ggplot(tmp2, aes_string(paste0("rna.", line, "_dir"), paste0("h3k27ac.", line, "_l2fc"), fill = paste0("rna.", line, "_dir"))) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
      xlab(label = paste0("RNA")) + ylab(label = paste0(" H3K27ac log2(FC)")) +
      ylim(c(-5,5))
    ggsave(paste0(line, ".rna.h3k27ac.", name, ".png"), mp, width = 2, height = 4)
    
    tmp2 <- tmp[get(paste0("atac.", line, "_padj")) < 0.05 & get(paste0("rna.", line, "_padj")) < 0.05]
    rna.atac.p <- wilcox.test(tmp2[get(paste0("rna.", line, "_dir")) == "down", get(paste0("atac.", line, "_l2fc"))],
                         tmp2[get(paste0("rna.", line, "_dir")) == "up", get(paste0("atac.", line, "_l2fc"))])$p.value
    mp <- ggplot(tmp2, aes_string(paste0("rna.", line, "_dir"), paste0("atac.", line, "_l2fc"), fill = paste0("rna.", line, "_dir"))) +
      geom_boxplot(outlier.shape = NA, show.legend = FALSE) + theme_bw() +
      xlab(label = paste0("RNA")) + ylab(label = paste0("ATAC log2(FC)")) +
      ylim(c(-5,5))
    ggsave(paste0(line, ".rna.atac.", name, ".png"), mp, width = 2, height = 4)
    
    res.inner <- c(atac.rna.p, histone.rna.p, atac.histone.p, rna.histone.p, rna.atac.p)
    names(res.inner) <- c("atac.rna.p", "histone.rna.p", "atac.histone.p", "rna.histone.p", "rna.atac.p")
    res[[line]] <- res.inner
  }
  return(as.data.table(do.call(rbind, res), keep.rownames = "comp"))
}


#### Randomly shuffle gene-peak interactions to test for significance ####
tmp <- inter.final.lin.enh.gene[,.(peakid, target)]
length(unique((tmp$peakid)))

peak.shuffle <- sample(tmp$peakid, nrow(tmp), replace = FALSE)
identical(sort(tmp$peakid), sort(peak.shuffle))
#sample(tmp$target, nrow(tmp), replace = TRUE)
tmp[, peakid := peak.shuffle]
tmp <- merge(tmp, res.list[["atac"]][,.(peakid, atac.brm_l2fc = brm_l2fc,
                                        atac.brm_padj = brm_padj, atac.lm_l2fc = lm_l2fc,
                                        atac.lm_padj = lm_padj)], by = "peakid")
tmp <- merge(tmp, rna[,.(target = ensembl, rna.brm_l2fc = brm_l2fc,
                         rna.brm_padj = brm_padj, rna.lm_l2fc = lm_l2fc,
                         rna.lm_padj = lm_padj, hgnc_symbol = hgnc)], by = "target")
tmp <- tmp[!duplicated(tmp[,.(peakid, target)])]
tmp <- merge(tmp, chip.annotations[["atac.h3k27ac.broad.overlap"]][,.(peakid = atac.peakid, h3k27ac.peakid)],
             by = "peakid")
tmp <- merge(tmp, res.list[["h3k27ac.broad"]][,.(h3k27ac.peakid = peakid,
                                                 h3k27ac.brm_l2fc = brm_l2fc,
                                                 h3k27ac.brm_padj = brm_padj,
                                                 h3k27ac.lm_l2fc = lm_l2fc,
                                                 h3k27ac.lm_padj = lm_padj)], by = "h3k27ac.peakid")
#fwrite(testAssociations(tmp, "shuffle.histone-nonsig", 1), "shuffle.peaks.pval.histone-nonsig.txt", sep = "\t")
fwrite(testAssociations(tmp, "shuffle.histone-sig", 0.05), "shuffle.peaks.pval.histone-sig.txt", sep = "\t")
#fwrite(testAssociations(inter.final.lin.enh.gene, "actual.histone-nonsig", 1), "actual.peaks.pval.histone-nonsig.txt", sep = "\t")
fwrite(testAssociations(inter.final.lin.enh.gene, "actual.histone-sig", 0.05), "actual.peaks.pval.histone-sig.txt", sep = "\t")
#### Randomly shuffle gene-peak interactions to test for significance ####