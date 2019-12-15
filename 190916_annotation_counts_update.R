# Last Updated: 
# Author: Wesley Cai
# Purpose: Plot histone changes counts with annotations

library(data.table)
library(ggplot2)

# load gsea
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/annotation"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))


#### Enhancer ####
peaks.changed <- res.list[["h3k27ac.broad"]][peakid %in% chip.annotations[["h3k27ac.broad.enh.active"]],.(peakid)]
hist.pthr <- 0.05

for(line in c("lm", "brm")){
  peaks.changed[peakid %in% res.list[["h3k27ac.broad"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) > 0, get("peakid")], paste0(line, ".enh") := "up"]
  peaks.changed[peakid %in% res.list[["h3k27ac.broad"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) < 0, get("peakid")], paste0(line, ".enh") := "down"]
}

gct.freq <- data.table(direction = rep(c("up", "down"), 3), type = c("lm", "lm", "shared", "shared", "brm", "brm"))

for(direct in c("up", "down")){
  lm.enh <- peaks.changed[lm.enh == direct,get("peakid")]
  brm.enh <- peaks.changed[brm.enh == direct,get("peakid")]
  shared.enh <- intersect(lm.enh, brm.enh)
  
  gct.freq[type == "lm" & direction == direct, enh := length(which(!(lm.enh %in% shared.enh)))]
  gct.freq[type == "brm" & direction == direct, enh := length(which(!(brm.enh %in% shared.enh)))]
  gct.freq[type == "shared" & direction == direct, enh := length(shared.enh)]
}

gct.freq[type == "lm",type.fix := "LM2"]
gct.freq[type == "shared",type.fix := "Shared"]
gct.freq[type == "brm",type.fix := "BrM2"]
gct.freq[direction == "up",Direction := "Up"]
gct.freq[direction == "down",Direction := "Down"]

gct.freq$Direction <- factor(gct.freq$Direction, c("Up", "Down"))
gct.freq$type.fix <- factor(gct.freq$type.fix, c("BrM2", "Shared", "LM2"))

gct.fract <- gct.freq[,.(enh.total = sum(enh)), by = .(type)]

gct.freq <- merge(gct.freq, gct.fract, by = "type")
gct.freq[direction == "up", enh.text := paste0(enh.total)]
gct.freq[direction == "up", enh.text.pos := enh.total + 150]

ymax <- 1800
mp <- ggplot(gct.freq, aes(type.fix, enh, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = enh.text, y = enh.text.pos), color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Peaks (n)") +
  ylim(c(0, ymax)) +
  scale_fill_manual(values = c("royalblue1", "midnightblue")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
mp
ggsave("enhancer_update.png", mp, width = 5, height = 3)
#### Enhancer ####

#### Promoter ####
peaks.changed <- res.list[["h3k4me3.broad"]][peakid %in% chip.annotations[["h3k4me3.broad.prom"]],.(peakid)]
hist.pthr <- 0.05

for(line in c("lm", "brm")){
  peaks.changed[peakid %in% res.list[["h3k4me3.broad"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) > 0, get("peakid")], paste0(line, ".prom") := "up"]
  peaks.changed[peakid %in% res.list[["h3k4me3.broad"]][get(paste0(line, "_padj")) < 0.05 & get(paste0(line, "_l2fc")) < 0, get("peakid")], paste0(line, ".prom") := "down"]
}

gct.freq <- data.table(direction = rep(c("up", "down"), 3), type = c("lm", "lm", "shared", "shared", "brm", "brm"))

for(direct in c("up", "down")){
  lm.prom <- peaks.changed[lm.prom == direct,get("peakid")]
  brm.prom <- peaks.changed[brm.prom == direct,get("peakid")]
  shared.prom <- intersect(lm.prom, brm.prom)
  
  gct.freq[type == "lm" & direction == direct, prom := length(which(!(lm.prom %in% shared.prom)))]
  gct.freq[type == "brm" & direction == direct, prom := length(which(!(brm.prom %in% shared.prom)))]
  gct.freq[type == "shared" & direction == direct, prom := length(shared.prom)]
}

gct.freq[type == "lm",type.fix := "LM2"]
gct.freq[type == "shared",type.fix := "Shared"]
gct.freq[type == "brm",type.fix := "BrM2"]
gct.freq[direction == "up",Direction := "Up"]
gct.freq[direction == "down",Direction := "Down"]

gct.freq$Direction <- factor(gct.freq$Direction, c("Up", "Down"))
gct.freq$type.fix <- factor(gct.freq$type.fix, c("BrM2", "Shared", "LM2"))

gct.fract <- gct.freq[,.(prom.total = sum(prom)), by = .(type)]

gct.freq <- merge(gct.freq, gct.fract, by = "type")
gct.freq[direction == "up", prom.text := paste0(prom.total)]
gct.freq[direction == "up", prom.text.pos := prom.total + 20]

ymax <- 250
mp <- ggplot(gct.freq, aes(type.fix, prom, fill = Direction)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = prom.text, y = prom.text.pos), color = "black") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(label = "") +
  ylab(label = "Peaks (n)") +
  ylim(c(0, ymax)) +
  scale_fill_manual(values = c("lightgreen", "forestgreen")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
mp
ggsave("promoter_update.png", mp, width = 5, height = 3)
#### Promoter ####
