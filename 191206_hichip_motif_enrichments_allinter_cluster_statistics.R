# Last Updated: 
# Author: Wesley Cai
# Purpose: 

library(data.table)
library(ggplot2)
library(openxlsx)
library(ggrepel)
library(VennDiagram)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/"
inputfolder <- "input/"
outputfolder <- "output/motif_allinter_cluster_statistics"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

load(file = "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_cluster/hichip.motif.res.list.filter.RData")
# Filter by expressed genes and significance
cluster.dict <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/motif/RSAT/output/rsat.melt.dict.txt")
cluster.dict <- cluster.dict[,.(cluster, motif = tf, final.id)]
cluster.dict <- cluster.dict[!duplicated(cluster.dict)]
cluster.dict[, motif := toupper(motif)]
all.genes <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/genomes/genes/output/ensembl_hgnc_map.txt")
cluster.dict <- cluster.dict[motif %in% all.genes$hgnc_symbol]

#### Functions ####
writeworkbook <- function(reslist, name){
  tmp.motif <- reslist
  for(line in c("lm", "brm")){
    wb <- createWorkbook(paste0(line, "_", name))
    databar.cols <- grep("neglogpadj|l2fc", colnames(tmp.motif[[line]][[comp]]$up))
    for(comp in names(tmp.motif[[line]])){
      comp.name <- gsub("atac", "a", gsub("histones", "h", gsub("changed_", "", gsub("_inter", "", comp))))
      comp.name <- gsub("changed_", "", gsub("_inter", "", comp))
      
      addWorksheet(wb, sheetName = paste0(comp.name, "_up"), zoom = 150)
      writeData(wb, paste0(comp.name, "_up"), tmp.motif[[line]][[comp]]$up)
      conditionalFormatting(wb,paste0(comp.name, "_up"),cols = databar.cols, rows =2:100,type = 'databar', gradient = FALSE)
      addWorksheet(wb, sheetName = paste0(comp.name, "_down"), zoom = 150)
      writeData(wb, paste0(comp.name, "_down"), tmp.motif[[line]][[comp]]$down)
      conditionalFormatting(wb,paste0(comp.name, "_down"),cols = databar.cols, rows =2:100,type = 'databar', gradient = FALSE)
    }
    saveWorkbook(wb, paste0(line, "_", name, ".xlsx"), overwrite = TRUE)
  }
}
#reslist <- hichip.motif.genes
#name <- "test"
writevenndiagram <- function(reslist, name){
  hichip.motif <- reslist
  for(direction in c("up", "down")){
    brm <- unique(hichip.motif[["brm"]][["changed_atac_inter_distal_prom"]][[direction]][,get("motif")])
    lm <- unique(hichip.motif[["lm"]][["changed_atac_inter_distal_prom"]][[direction]][,get("motif")])
    express <- res.list[["rna"]][hgnc %in% c(brm, lm),]
    express <- merge(express, cluster.dict[,.(hgnc = motif, final.id)], by = "hgnc")
    brm.specific <- brm[!(brm %in% intersect(brm, lm))]
    lm.specific <- lm[!(lm %in% intersect(brm,lm))]
    shared <- intersect(brm,lm)
    fwrite(express[hgnc %in% brm.specific], paste0(name, "_", direction, ".brm.specific.txt"), sep = "\t", col.names = TRUE)
    fwrite(express[hgnc %in% lm.specific], paste0(name, "_", direction, ".lm.specific.txt"), sep = "\t", col.names = TRUE)
    fwrite(express[hgnc %in% shared], paste0(name, "_", direction, ".shared.txt"), sep = "\t", col.names = TRUE)
    pdf(paste0(name, "_", direction, ".venn.pdf"), width = 2, height = 2)
    if(length(lm) < length(brm)){
      manual.rotate <- 180
    } else {
      manual.rotate <- 0
    }
    draw.pairwise.venn(area1 = length(lm), area2 = length(brm), cross.area = length(shared),
                       category = c("", ""), rotation.degree = manual.rotate)
    dev.off()
    
    pdf(paste0(name, "_", direction, ".venn.nolabel.pdf"), width = 2, height = 2)
    if(length(lm) < length(brm)){
      manual.rotate <- 180
    } else {
      manual.rotate <- 0
    }
    draw.pairwise.venn(area1 = length(lm), area2 = length(brm), cross.area = length(shared),
                       category = c("", ""), rotation.degree = manual.rotate, cex = 0, ext.text = FALSE)
    dev.off()
  }
}
#### Functions ####


# Map motif back to genes, all motifs regardless of expression
# hichip.motif.genes <- hichip.cluster.res.list.sig
hichip.motif.genes <- hichip.motif.res.list.filter

# Only take motifs that occur in either sig up or sig down peaks (take the more significant one)
line <- "lm"
comp <- "changed_atac_inter_distal_prom"
# test <- list()
# test$lm$changed_atac_inter_distal_prom$up <- data.table(motif = c("TEAD", "FOX", "AP2"),
#                                                         p.hyper = c(0.01, 0.001, 0.1))
# test$lm$changed_atac_inter_distal_prom$down <- data.table(motif = c("TEAD", "FOX", "AP2"),
#                                                         p.hyper = c(0.02, 0.0001, 0.2))
# hichip.motif.genes <- test

for(line in names(hichip.motif.genes)){
  for(comp in names(hichip.motif.genes[[line]])){
    common.motif <- intersect(hichip.motif.genes[[line]][[comp]]$up$motif,
                              hichip.motif.genes[[line]][[comp]]$down$motif)
    up.dt <- hichip.motif.genes[[line]][[comp]]$up[motif %in% common.motif]
    up.unique.dt <- hichip.motif.genes[[line]][[comp]]$up[!(motif %in% common.motif)]
    down.dt <- hichip.motif.genes[[line]][[comp]]$down[motif %in% common.motif]
    down.unique.dt <- hichip.motif.genes[[line]][[comp]]$down[!(motif %in% common.motif)]
    setkey(up.dt, "motif")
    setkey(down.dt, "motif")
    up.dt.update <- up.dt[up.dt[common.motif, get("p.hyper")] < down.dt[common.motif, get("p.hyper")]]
    down.dt.update <- down.dt[!(up.dt[common.motif, get("p.hyper")] < down.dt[common.motif, get("p.hyper")])]
    up.dt.update <- rbind(up.dt.update, up.unique.dt)
    up.dt.update <- up.dt.update[order(p.hyper),]
    down.dt.update <- rbind(down.dt.update, down.unique.dt)
    down.dt.update <- down.dt.update[order(p.hyper),]
    
    hichip.motif.genes[[line]][[comp]]$up <- up.dt.update
    hichip.motif.genes[[line]][[comp]]$down <- down.dt.update
  }
}


writeworkbook(reslist = hichip.motif.genes, name = "all_cluster")
writevenndiagram(reslist = hichip.motif.genes, name = "all_cluster")

hichip.genes <- hichip.motif.genes[["lm"]][["changed_atac_inter_distal_prom"]][["up"]]

distal.prom.motif.res <- data.table()
for(line in c("lm", "brm")){
  for(direction in c("up", "down")){
    tmp <- hichip.motif.genes[[line]][["changed_atac_inter_distal_prom"]][[direction]][,.(hgnc = motif, class = final.id, cluster, padj.hyper, neglogpadj.hyper)]
    tmp$comp <- paste0(line, ".", direction)
    distal.prom.motif.res <- rbind(distal.prom.motif.res, tmp)
  }
}
distal.prom.motif.res <- as.data.table(reshape2::dcast(distal.prom.motif.res, formula = hgnc + class ~ comp, value.var = "padj.hyper"))
length(unique(distal.prom.motif.res$hgnc))
distal.prom.motif.res <- merge(distal.prom.motif.res, res.list[["rna"]], by = "hgnc")
# Only take one motif per gene (use lower cluster number first)
distal.prom.motif.res[, cluster_no := as.numeric(gsub(".*_", "", class))]
distal.prom.motif.res <- distal.prom.motif.res[order(cluster_no)]
distal.prom.motif.res <- distal.prom.motif.res[!duplicated(hgnc)]
save(distal.prom.motif.res, file = "distal.prom.motif.res.RData")


