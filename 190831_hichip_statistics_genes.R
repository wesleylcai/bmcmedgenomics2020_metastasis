# Last Updated: 
# Author: Wesley Cai
# Purpose: Assign genes to changed enhancers and promoter histone mod by significance and direction

library(data.table)
library(ggplot2)
library(reshape2)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/"
inputfolder <- "input/"
outputfolder <- "annotation/"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.all.RData")

# Remove prom-prom interactions of same genes
tmp.filter <- inter.final.all.dt
tmp.filter[is.na(annote1), annote1 := "unknown"]
tmp.filter[is.na(annote2), annote2 := "unknown"]
tmp.filter[, gene1 := gsub("_.*", "", gene1)]
tmp.filter[, gene2 := gsub("_.*", "", gene2)]
tmp.filter <- tmp.filter[!(gene1 == gene2 & annote1 == "promoter" & annote2 == "promoter"),]
tmp.filter.lin <- rbind(tmp.filter[annote1 == "promoter",.(gene = gene1, target.atac = atac2.peakid, target.annote = annote2, intercount.sum)],
                        tmp.filter[annote2 == "promoter",.(gene = gene2, target.atac = atac1.peakid, target.annote = annote1, intercount.sum)])
tmp.filter.lin <- tmp.filter.lin[target.annote == "enhancer",]

# Get H3K27ac enhancers
atac.h3k27ac <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/atac.h3k27ac.broad.bed")
atac.h3k27ac <- atac.h3k27ac[V9 == 0,]
tmp.h3k27ac <- res.list[["h3k27ac.broad"]]
colnames(tmp.h3k27ac) <- paste0("h3k27ac.", colnames(tmp.h3k27ac))
tmp.atac <- res.list[["atac"]]
colnames(tmp.atac) <- paste0("atac.", colnames(tmp.atac))
atac.h3k27ac <- merge(atac.h3k27ac[,.(atac.peakid = V4, h3k27ac.peakid = V8)], tmp.h3k27ac, by = "h3k27ac.peakid")
tmp.filter.lin <- merge(tmp.filter.lin, atac.h3k27ac, by.x = "target.atac", by.y = "atac.peakid")
tmp.filter.lin <- merge(tmp.filter.lin, tmp.atac, by.x = "target.atac", by.y = "atac.peakid")

# Get H3K4me3 promoters
atac.h3k4me3.tss <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/atac.h3k4me3.broad.tss.bed")
atac.h3k4me3.tss <- atac.h3k4me3.tss[V9 == 0 & V15 == 0]
atac.h3k4me3.tss[, V13 := gsub("_.*", "", V13)]
genes.prom <- merge(res.list[["rna"]][,.(ensembl)], atac.h3k4me3.tss[,.(ensembl = V13, atac.peakid = V4, h3k4me3.peakid = V8)], by = "ensembl")
tmp.h3k4me3 <- res.list[["h3k4me3.broad"]]
colnames(tmp.h3k4me3) <- paste0("h3k4me3.", colnames(tmp.h3k4me3))
genes.prom <- merge(genes.prom, tmp.h3k4me3, by = "h3k4me3.peakid")
genes.prom <- merge(genes.prom, tmp.atac, by = "atac.peakid")

genes.changed <- list()
genes.changed[["lm"]][["up"]] <- res.list[["rna"]][lm_padj < 0.05 & lm_l2fc > 0, .(ensembl)]
genes.changed[["lm"]][["down"]] <- res.list[["rna"]][lm_padj < 0.05 & lm_l2fc < 0, .(ensembl)]
genes.changed[["brm"]][["up"]] <- res.list[["rna"]][brm_padj < 0.05 & brm_l2fc > 0, .(ensembl)]
genes.changed[["brm"]][["down"]] <- res.list[["rna"]][brm_padj < 0.05 & brm_l2fc < 0, .(ensembl)]

line <- "lm"
direction <- "up"
# This takes significantly changed ATAC peaks and directionally changed H3K27ac and H3K4me3 peaks
for(line in c("lm", "brm")){
  for(direction in c("up", "down")){
    tmp <- genes.changed[[line]][[direction]]
    tmp[,enhancer.change := 0]
    tmp[,promoter.change := 0]
    if(direction == "up"){
      # Enhancers
      enhancer.link <- tmp.filter.lin
      # goi <- enhancer.link[get(paste0("h3k27ac.", line, "_l2fc")) > 0 & get(paste0("atac.", line, "_l2fc")) > 0 &
      #                        get(paste0("atac.", line, "_padj")) < 0.05, get("gene")]
      goi <- enhancer.link[get(paste0("h3k27ac.", line, "_l2fc")) > 0 & get(paste0("h3k27ac.", line, "_padj")) < 0.05, get("gene")]
      tmp[ensembl %in% goi, enhancer.change := 1]
      # Promoters
      promoter.link <- genes.prom
      # goi <- promoter.link[get(paste0("h3k4me3.", line, "_l2fc")) > 0 & get(paste0("atac.", line, "_l2fc")) > 0 &
      #                        get(paste0("atac.", line, "_padj")) < 0.05, get("ensembl")]
      goi <- promoter.link[get(paste0("h3k4me3.", line, "_l2fc")) > 0 & get(paste0("h3k4me3.", line, "_padj")) < 0.05, get("ensembl")]
      tmp[ensembl %in% goi, promoter.change := 1]
    } else if(direction == "down"){
      # Enhancers
      enhancer.link <- tmp.filter.lin
      # goi <- enhancer.link[get(paste0("h3k27ac.", line, "_l2fc")) < 0 & get(paste0("atac.", line, "_l2fc")) < 0 &
      #                        get(paste0("atac.", line, "_padj")) < 0.05, get("gene")]
      goi <- enhancer.link[get(paste0("h3k27ac.", line, "_l2fc")) < 0 & get(paste0("h3k27ac.", line, "_padj")) < 0.05, get("gene")]
      tmp[ensembl %in% goi, enhancer.change := 1]
      # Promoters
      promoter.link <- genes.prom
      # goi <- promoter.link[get(paste0("h3k4me3.", line, "_l2fc")) < 0 & get(paste0("atac.", line, "_l2fc")) < 0 &
      #                        get(paste0("atac.", line, "_padj")) < 0.05, get("ensembl")]
      goi <- promoter.link[get(paste0("h3k4me3.", line, "_l2fc")) < 0 & get(paste0("h3k4me3.", line, "_padj")) < 0.05, get("ensembl")]
      tmp[ensembl %in% goi, promoter.change := 1]
    }
    genes.changed[[line]][[direction]] <- tmp
  }
}

gct <- res.list[["rna"]][,.(ensembl)]
line <- "lm"
direction <- "up"
for(line in c("lm", "brm")){
  for(direction in c("up", "down")){
    gct[, eval(paste0(line, ".", direction)) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][, get("ensembl")], eval(paste0(line, ".", direction)) := 1]
    gct[, eval(paste0(line, ".enh.", direction)) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][enhancer.change == 1, get("ensembl")], eval(paste0(line, ".enh.", direction)) := 1]
    gct[, eval(paste0(line, ".prom.", direction)) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][promoter.change == 1, get("ensembl")], eval(paste0(line, ".prom.", direction)) := 1]
    
    gct[, eval(paste0(line)) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][, get("ensembl")], eval(paste0(line)) := 1]
    gct[, eval(paste0(line, ".enh")) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][enhancer.change == 1, get("ensembl")], eval(paste0(line, ".enh")) := 1]
    gct[, eval(paste0(line, ".prom")) := 0]
    gct[ensembl %in% genes.changed[[line]][[direction]][promoter.change == 1, get("ensembl")], eval(paste0(line, ".prom")) := 1]
  }
}

gct.freq <- data.table(type = c("lm", "shared", "brm"))
gct.freq[type == "lm", enhancer := nrow(gct[lm.enh == 1,])]
gct.freq[type == "lm", enhancer.fract := nrow(gct[lm.enh == 1,])/nrow(gct[lm == 1,])]
gct.freq[type == "shared", enhancer := nrow(gct[lm.enh == 1 & brm.enh == 1,])]
gct.freq[type == "shared", enhancer.fract := nrow(gct[lm.enh == 1 & brm.enh == 1,])/nrow(gct[lm == 1 & brm == 1,])]
gct.freq[type == "brm", enhancer := nrow(gct[brm.enh == 1,])]
gct.freq[type == "brm", enhancer.fract := nrow(gct[brm.enh == 1,])/nrow(gct[brm == 1,])]
gct.freq[type == "lm", promoter := nrow(gct[lm.prom == 1,])]
gct.freq[type == "lm", promoter.fract := nrow(gct[lm.prom == 1,])/nrow(gct[lm == 1,])]
gct.freq[type == "shared", promoter := nrow(gct[lm.prom == 1 & brm.prom == 1,])]
gct.freq[type == "shared", promoter.fract := nrow(gct[lm.prom == 1 & brm.prom == 1,])/nrow(gct[lm == 1 & brm == 1,])]
gct.freq[type == "brm", promoter := nrow(gct[brm.prom == 1,])]
gct.freq[type == "brm", promoter.fract := nrow(gct[brm.prom == 1,])/nrow(gct[brm == 1,])]
gct.freq$type <- factor(c("LM2", "Shared", "BrM2"), c("BrM2", "Shared", "LM2"))
gct.freq[, enhancer.text := paste0(enhancer, " (", signif(enhancer.fract, 3)*100, "%)")]
gct.freq[, promoter.text := paste0(promoter, " (", signif(promoter.fract, 3)*100, "%)")]

ymax <- 850
mp <- ggplot(gct.freq, aes(type, enhancer)) +
  geom_bar(stat = "identity", fill = "darkorange", color = "black") +
  geom_text(aes(label = enhancer.text), nudge_y = ymax/8, color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Genes (n)") +
  ylim(c(0, ymax))
mp
ggsave("enhancer.genes.png", mp, width = 5, height = 3)

ymax <- 600
mp <- ggplot(gct.freq, aes(type, promoter)) +
  geom_bar(stat = "identity", fill = "dodgerblue3", color = "black") +
  geom_text(aes(label = promoter.text), nudge_y = ymax/8, color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Genes (n)") +
  ylim(c(0, ymax))
mp
ggsave("promoter.genes.png", mp, width = 5, height = 3)
#### Plot all fractions combined ####
inter.count.select <- 2
fract.dt.melt.combined <- data.table()
for(line in names(inter.count.res)){
  for(direction in names(inter.count.res[[line]])){
    tmp <- inter.count.res[[line]][[direction]][inter.count == inter.count.select,]
    tmp[, group := paste(line, direction, sep = ".")]
    fract.dt.melt.combined <- rbind(fract.dt.melt.combined, tmp)
  }
}

mp <- ggplot(fract.dt.melt.combined, aes(group, value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d()
mp
ggsave(paste0(inter.count.select, ".intercount.all_fractions.png"), mp, width = 6, height = 5)
#### Plot all fractions combined ####

