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
outputfolder <- "annotation_revision/"
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
# atac.h3k4me3.tss <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/atac.h3k4me3.broad.tss.bed")
# atac.h3k4me3.tss <- atac.h3k4me3.tss[V9 == 0 & V15 == 0]
# atac.h3k4me3.tss[, V13 := gsub("_.*", "", V13)]
# genes.prom <- merge(res.list[["rna"]][,.(ensembl)], atac.h3k4me3.tss[,.(ensembl = V13, atac.peakid = V4, h3k4me3.peakid = V8)], by = "ensembl")
# tmp.h3k4me3 <- res.list[["h3k4me3.broad"]]
# colnames(tmp.h3k4me3) <- paste0("h3k4me3.", colnames(tmp.h3k4me3))
# genes.prom <- merge(genes.prom, tmp.h3k4me3, by = "h3k4me3.peakid")
# genes.prom <- merge(genes.prom, tmp.atac, by = "atac.peakid")

# Get H3K4me3 promoters w/o ATAC
h3k4me3.tss <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me3_merge_peaks.tss.bed")
h3k4me3.tss <- h3k4me3.tss[V10 == 0,]
h3k4me3.tss[, ensembl := gsub("_.*", "", V8)]
h3k4me3.tss.prom <- merge(h3k4me3.tss[,.(h3k4me3.peakid = V4, ensembl)], res.list[["h3k4me3.broad"]], by.x = "h3k4me3.peakid", by.y = "peakid")
h3k4me3.tss.prom <- h3k4me3.tss.prom[,.SD[which.min(met_padj)], by = .(h3k4me3.peakid, ensembl)]

genes.changed <- res.list[["rna"]][,.(ensembl)]
hist.pthr <- 0.05
for(line in c("lm", "brm")){
  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) > 0, get("ensembl")], paste0(line, ".rna") := "up"]
  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) < 0, get("ensembl")], paste0(line, ".rna") := "down"]

  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) > 0, get("ensembl")] &
                  ensembl %in% tmp.filter.lin[get(paste0("h3k27ac.", line, "_padj")) < hist.pthr & get(paste0("h3k27ac.", line, "_l2fc")) > 0, get("gene")], paste0(line, ".enh") := "up"]
  message(nrow(genes.changed[get(paste0(line, ".enh")) == "up"]))
  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) < 0, get("ensembl")] &
                  ensembl %in% tmp.filter.lin[get(paste0("h3k27ac.", line, "_padj")) < hist.pthr & get(paste0("h3k27ac.", line, "_l2fc")) < 0, get("gene")], paste0(line, ".enh") := "down"]
  message(nrow(genes.changed[get(paste0(line, ".enh")) == "up"]))
  
  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) > 0, get("ensembl")] &
                  ensembl %in% h3k4me3.tss.prom[get(paste0(line, "_padj")) < hist.pthr & get(paste0(line, "_l2fc")) > 0, get("ensembl")], paste0(line, ".prom") := "up"]
  genes.changed[ensembl %in% res.list[["rna"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) < 0, get("ensembl")] &
                  ensembl %in% h3k4me3.tss.prom[get(paste0(line, "_padj")) < hist.pthr & get(paste0(line, "_l2fc")) < 0, get("ensembl")], paste0(line, ".prom") := "down"]
}

gct.freq <- data.table(direction = rep(c("up", "down"), 3), type = c("lm", "lm", "shared", "shared", "brm", "brm"))

nrow(genes.changed[brm.rna == "down" | brm.rna == "up",])
nrow(genes.changed[lm.rna == "down" | lm.rna == "up",])

for(direct in c("up", "down")){
  lm.rna <- genes.changed[lm.rna == direct,get("ensembl")]
  brm.rna <- genes.changed[brm.rna == direct,get("ensembl")]
  shared.rna <- intersect(lm.rna, brm.rna)
  lm.enh <- genes.changed[lm.enh == direct & lm.rna == direct,get("ensembl")]
  brm.enh <- genes.changed[brm.enh == direct & brm.rna == direct,get("ensembl")]
  shared.enh <- intersect(lm.enh, brm.enh)
  lm.prom <- genes.changed[lm.prom == direct & lm.rna == direct,get("ensembl")]
  brm.prom <- genes.changed[brm.prom == direct & brm.rna == direct,get("ensembl")]
  shared.prom <- intersect(lm.prom, brm.prom)
  
  gct.freq[type == "lm" & direction == direct, enh := length(which(!(lm.enh %in% shared.enh)))]
  gct.freq[type == "brm" & direction == direct, enh := length(which(!(brm.enh %in% shared.enh)))]
  gct.freq[type == "shared" & direction == direct, enh := length(shared.enh)]
  # gct.freq[type == "lm" & direction == direct, enh.fract := length(which(!(lm.enh %in% shared.enh)))/length(which(!(lm.rna %in% shared.rna)))]
  # gct.freq[type == "brm" & direction == direct, enh.fract := length(which(!(brm.enh %in% shared.enh)))/length(which(!(brm.rna %in% shared.rna)))]
  # gct.freq[type == "shared" & direction == direct, enh.fract := length(shared.enh)/length(shared.rna)]

  gct.freq[type == "lm" & direction == direct, prom := length(which(!(lm.prom %in% shared.prom)))]
  gct.freq[type == "brm" & direction == direct, prom := length(which(!(brm.prom %in% shared.prom)))]
  gct.freq[type == "shared" & direction == direct, prom := length(shared.prom)]
  # gct.freq[type == "lm" & direction == direct, prom.fract := length(which(!(lm.prom %in% shared.prom)))/length(which(!(lm.rna %in% shared.rna)))]
  # gct.freq[type == "brm" & direction == direct, prom.fract := length(which(!(brm.prom %in% shared.prom)))/length(which(!(brm.rna %in% shared.rna)))]
  # gct.freq[type == "shared" & direction == direct, prom.fract := length(shared.prom)/length(shared.rna)]
  # 
}

gct.freq[type == "lm",type.fix := "LM2"]
gct.freq[type == "shared",type.fix := "Shared"]
gct.freq[type == "brm",type.fix := "BrM2"]
gct.freq[direction == "up",Direction := "Up"]
gct.freq[direction == "down",Direction := "Down"]

gct.freq$Direction <- factor(gct.freq$Direction, c("Up", "Down"))
gct.freq$type.fix <- factor(gct.freq$type.fix, c("BrM2", "Shared", "LM2"))

gct.fract <- gct.freq[,.(enh.total = sum(enh), prom.total = sum(prom)), by = .(type)]

lm.rna <- genes.changed[lm.rna == "up" | lm.rna == "down",get("ensembl")]
brm.rna <- genes.changed[brm.rna == "up" | brm.rna == "down",get("ensembl")]
shared.rna <- intersect(lm.rna, brm.rna)
lm.rna.len <- length(which(!(lm.rna %in% shared.rna)))
brm.rna.len <- length(which(!(brm.rna %in% shared.rna)))
shared.rna.len <- length(shared.rna)

gct.fract[type == "lm", enh.fract := enh.total/lm.rna.len]
gct.fract[type == "brm", enh.fract := enh.total/brm.rna.len]
gct.fract[type == "shared", enh.fract := enh.total/shared.rna.len]
gct.fract[type == "lm", prom.fract := prom.total/lm.rna.len]
gct.fract[type == "brm", prom.fract := prom.total/brm.rna.len]
gct.fract[type == "shared", prom.fract := prom.total/shared.rna.len]

# Total genes linked to changes
prom.genes <- sum(gct.fract$prom.total)
prom.genes.fract <- sum(gct.fract$prom.total)/(lm.rna.len + brm.rna.len + shared.rna.len)
enh.genes <- sum(gct.fract$enh.total)
enh.genes.fract <- sum(gct.fract$enh.total)/(lm.rna.len + brm.rna.len + shared.rna.len)


gct.freq <- merge(gct.freq, gct.fract, by = "type")
gct.freq[direction == "up", enh.text := paste0(enh.total, " (", signif(enh.fract, 3)*100, "%)")]
gct.freq[direction == "up", prom.text := paste0(prom.total, " (", signif(prom.fract, 3)*100, "%)")]
gct.freq[direction == "up", enh.text.pos := enh.total + 100]
gct.freq[direction == "up", prom.text.pos := prom.total + 50]

gct.freq[type == "lm",type.fix := "LM2"]
gct.freq[type == "shared",type.fix := "Shared"]
gct.freq[type == "brm",type.fix := "BrM2"]
gct.freq[direction == "up",Direction := "Up"]
gct.freq[direction == "down",Direction := "Down"]


ymax <- 600
mp <- ggplot(gct.freq, aes(type.fix, enh, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = enh.text, y = enh.text.pos), color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Genes (n)") +
  ylim(c(0, ymax)) +
  scale_fill_manual(values = c("royalblue1", "midnightblue")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
mp
ggsave("enhancer.genes_update.png", mp, width = 5, height = 3)

ymax <- 250
mp <- ggplot(gct.freq, aes(type.fix, prom, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = prom.text, y = prom.text.pos), color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Genes (n)") +
  ylim(c(0, ymax)) +
  scale_fill_manual(values = c("lightgreen", "forestgreen")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
mp
ggsave("promoter.genes_update.png", mp, width = 5, height = 3)

