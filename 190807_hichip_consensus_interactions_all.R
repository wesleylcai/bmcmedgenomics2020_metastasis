# Last Updated: 
# Author: Wesley Cai
# Purpose: I can use this for LM2 because it is technically a subpopulation of 231

library(data.table)
library(ggplot2)
library(enrichR)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/"
inputfolder <- "input/"
outputfolder <- "output/"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#### Proof of concept for finding unique interactions ####
# library(data.table)
# dt1 <- data.table(V1 = c("A", "B", "A", "A", "C", "A"),
#                   V2 = c("B", "A", "B", "B", "B", "C"))
# dt1[, id := paste0(V1, "_", V2)]
# dt2 <- dt1[,.(id, V1 = V2, V2 = V1)]
# dt3 <- rbind(dt1[,.(id, V1, V2)], dt2[,.(id, V1, V2)])
# dt4 <- dt3[,.(id, sum = .N), by = .(V1, V2)]
# dt5 <- dt4[!(duplicated(id)) & !(duplicated(dt4[,.(V1, V2)]))]
#### Proof of concept for finding unique interactions ####

#### Get region-region interactions that are associated with H3K27ac or H3K4me3 ####
hc.reps <- list()
hc.reps[["rep1"]] <- fread("../associations/rep1.unique.intercount.txt.gz")
hc.reps[["rep2"]] <- fread("../associations/rep2.unique.intercount.txt.gz")
rep <- "rep1"
tmp.atac.tss <- chip.annotations[["atac.tss.promoters"]][,.(peakid, gene)]
# tmp.atac.tss[, gene := gsub("_.*", "", gene)]
# for(rep in names(hc.reps)){
#   # Count up all rows where the interactions are between two of the same peaks
#   tmp <- hc.reps[[rep]][,.(id = paste0(atac1.peakid, "_", atac2.peakid),
#                            atac1.peakid, atac2.peakid, intercount.sum)]
#   tmp2 <- rbind(tmp[,.(id, peak1 = atac1.peakid, peak2 = atac2.peakid, intercount.sum)],
#                 tmp[,.(id, peak1 = atac2.peakid, peak2 = atac1.peakid, intercount.sum)])
#   tmp2 <- tmp2[,.(id, sum.final = sum(intercount.sum)), by = .(peak1, peak2)]
#   hc.reps[[rep]] <- tmp2[!duplicated(id) & !duplicated(tmp2[,.(peak1, peak2)]),
#                          .(atac1.peakid = peak1, atac2.peakid = peak2, intercount.sum = sum.final)]
#   hc.reps[[rep]]
# 
#   hc.reps[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote1 := "enhancer"]
#   hc.reps[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote2 := "enhancer"]
#   hc.reps[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote1 := "promoter"]
#   hc.reps[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote2 := "promoter"]
#   hc.reps[[rep]] <- merge(hc.reps[[rep]], tmp.atac.tss[,.(peakid, gene1 = gene)], by.x = "atac1.peakid", by.y = "peakid", all.x = TRUE)
#   hc.reps[[rep]] <- merge(hc.reps[[rep]], tmp.atac.tss[,.(peakid, gene2 = gene)], by.x = "atac2.peakid", by.y = "peakid", all.x = TRUE)
#   hc.reps[[rep]][gene1 == gene2, equal := TRUE]
# }
# save(hc.reps, file = "hc.reps.RData")
load(file = "hc.reps.RData")
#### Get region-region interactions that are associated with H3K27ac or H3K4me3 ####

#### Get region-region interactions that exist in both replicates ####
# inter.rep1 <- paste0(hc.reps[["rep1"]]$atac1.peakid, "_", hc.reps[["rep1"]]$atac2.peakid)
# inter.rep2 <- paste0(hc.reps[["rep2"]]$atac1.peakid, "_", hc.reps[["rep2"]]$atac2.peakid)
# inter.rep2.swap <- paste0(hc.reps[["rep2"]]$atac2.peakid, "_", hc.reps[["rep2"]]$atac1.peakid)
# inter.rep1.rep2 <- intersect(inter.rep1, inter.rep2)
# inter.rep1.rep2.swap <- intersect(inter.rep1, inter.rep2.swap)
# inter.final <- unique(inter.rep1.rep2, inter.rep1.rep2.swap)
# 
# tmp <- hc.reps[["rep1"]]
# tmp[, id := paste0(atac1.peakid, "_", atac2.peakid)]
# inter.final.all.dt <- tmp[id %in% inter.final,.(atac1.peakid, atac2.peakid, annote1, annote2, gene1, gene2, intercount.sum, equal)]
# save(inter.final.all.dt, file = "inter.final.dt.all.RData")
load("inter.final.dt.all.RData")
load("inter.final.dt.RData")
#### Get region-region interactions that exist in both replicates ####

inter.count.freq <- as.data.table(table(inter.final.all.dt$intercount.sum))

plot(hist(inter.final.all.dt$intercount.sum), xlim = c(0, 67))
# Most interactions occur only once
hist(inter.final.all.dt$intercount.sum, breaks = 5000, xlim = c(0,67))
