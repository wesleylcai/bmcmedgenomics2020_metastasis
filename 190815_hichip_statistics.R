# Last Updated: 
# Author: Wesley Cai
# Purpose:

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

load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.all.RData")

## Total number of unique, replicated interactions
length(which(duplicated(c(paste0(inter.final.all.dt$atac1.peakid, inter.final.all.dt$atac2.peakid), paste0(inter.final.all.dt$atac2.peakid, inter.final.all.dt$atac1.peakid)))))
# Since ^ is 0 that means all interactions are unique and total inter is defined as nrow
nrow(inter.final.all.dt)
# 1659221 interactions

tmp <- inter.final.all.dt#[intercount.sum > 20]
# Percent of cross-chromosome interactions (only 0.652%)
tmp2 <- tmp[,.(atac1.peakid, atac2.peakid)]
tmp2[,atac1.chr := gsub("-.*", "", atac1.peakid)]
tmp2[,atac2.chr := gsub("-.*", "", atac2.peakid)]
nrow(tmp2[atac1.chr != atac2.chr])/nrow(tmp2)

# Remove prom-prom interactions of same genes
tmp.filter <- inter.final.all.dt
tmp.filter[is.na(annote1), annote1 := "unknown"]
tmp.filter[is.na(annote2), annote2 := "unknown"]
tmp.filter[, gene1 := gsub("_.*", "", gene1)]
tmp.filter[, gene2 := gsub("_.*", "", gene2)]
tmp.filter <- tmp.filter[!(gene1 == gene2 & annote1 == "promoter" & annote2 == "promoter"),]

#### Plot fractions with respect to connection strength ####
max.inter.count <- 50
inter.count.res <- list()
for(line in c("lm", "brm", "all")){
  for(direction in c("up", "down")){

    if(line == "all"){
      direction <- "all"
      tmp.filter.line <- tmp.filter
    } else {
      # Only keep changed peaks
      if(direction == "up"){
        changed.peaks <- res.list[["atac"]][get(paste0(line, "_l2fc")) > 0 & get(paste0(line, "_padj")) < 0.05,get("peakid")]
      } else {
        changed.peaks <- res.list[["atac"]][get(paste0(line, "_l2fc")) < 0 & get(paste0(line, "_padj")) < 0.05,get("peakid")]
      }
      tmp.filter.line <- tmp.filter[atac1.peakid %in% changed.peaks | atac2.peakid %in% changed.peaks,]
    }
    
    fract.dt <- data.table(inter.count = 0:max.inter.count)
    for(i in fract.dt$inter.count){
      message(i)
      tmp <- tmp.filter.line[intercount.sum > i]
      total.n <- nrow(tmp)
      prom.prom.n <- nrow(tmp[annote1 == "promoter" & annote2 == "promoter",])
      prom.enh.n <- nrow(tmp[(annote1 == "enhancer" & annote2 == "promoter") | (annote2 == "enhancer" & annote1 == "promoter"),])
      enh.enh.n <- nrow(tmp[annote1 == "enhancer" & annote2 == "enhancer",])
      unkown.unknown.n <- nrow(tmp[annote1 == "unknown" & annote2 == "unknown",])
      enh.unknown.n <- nrow(tmp[(annote1 == "enhancer" & annote2 == "unknown") | (annote2 == "enhancer" & annote1 == "unknown"),])
      prom.unknown.n <- nrow(tmp[(annote1 == "promoter" & annote2 == "unknown") | (annote2 == "promoter" & annote1 == "unknown"),])
      
      fract.dt[inter.count == i, prom.prom := prom.prom.n/total.n]
      fract.dt[inter.count == i, prom.enh := prom.enh.n/total.n]
      fract.dt[inter.count == i, enh.enh := enh.enh.n/total.n]
      fract.dt[inter.count == i, unkown.unknown := unkown.unknown.n/total.n]
      fract.dt[inter.count == i, enh.unknown := enh.unknown.n/total.n]
      fract.dt[inter.count == i, prom.unknown := prom.unknown.n/total.n]
    }
    
    fract.dt.melt <- melt(fract.dt, id.vars = "inter.count", measure.vars = colnames(fract.dt)[2:ncol(fract.dt)])
    mp <- ggplot(fract.dt.melt, aes(inter.count, value, fill = variable)) +
      geom_area() +
      xlim(c(0,max.inter.count)) +
      scale_fill_viridis_d()
    mp
    ggsave(paste0(max.inter.count, ".intercount.v.annote.fract.", line, ".", direction, ".png"), mp, width = 3, height = 5)
    
    inter.count.res[[line]][[direction]] <- fract.dt.melt
    if(line == "all") break
  }
}
#### Plot fractions with respect to connection strength ####

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

## Total number of promoter-promoter interactions: 156763 (9.45%)
tmp.count <- nrow(tmp[annote1 == "promoter" & annote2 == "promoter",])
tmp.count/nrow(tmp)

## Total number of promoter-enhancer interactions: 230804 (13.9%)
tmp.count <- nrow(tmp[(annote1 == "enhancer" & annote2 == "promoter") | (annote2 == "enhancer" & annote1 == "promoter"),])
tmp.count/nrow(tmp)

## Total number of enhancer-enhancer interactions: 364870 (22.0%)
tmp.count <- nrow(tmp[annote1 == "enhancer" & annote2 == "enhancer",])
tmp.count/nrow(tmp)

## Total number of NA-NA interactions: 173308 (10.4%)
tmp.count <- nrow(tmp[is.na(annote1) & is.na(annote2),])
tmp.count/nrow(tmp)

## Total number of enhancer-NA interactions: 442582 (27%)
tmp.count <- nrow(tmp[(annote1 == "enhancer" & is.na(annote2)) | (annote2 == "enhancer" & is.na(annote1)),])
tmp.count/nrow(tmp)

## Total number of promoter-NA interactions: 290894 (17.5%)
tmp.count <- nrow(tmp[(annote1 == "promoter" & is.na(annote2)) | (annote2 == "promoter" & is.na(annote1)),])
tmp.count/nrow(tmp)

## Gene interactions (all): 11961 (92.1%)
length(unique(c(tmp.filter$gene1, tmp.filter$gene2)))/
  length(unique(res.list[["rna"]][,get("ensembl")]))

## Gene interactions (enh/prom): 11399 (87.7%)
tmp.filter.annote <- tmp.filter[annote1 %in% c("promoter", "enhancer") & annote2 %in% c("promoter", "enhancer")]
length(unique(c(tmp.filter.annote$gene1, tmp.filter.annote$gene2)))/
  length(unique(res.list[["rna"]][,get("ensembl")]))
