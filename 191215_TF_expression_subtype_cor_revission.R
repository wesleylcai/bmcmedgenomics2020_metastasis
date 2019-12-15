# Last Updated: 
# Author: Wesley Cai
# Purpose: TF expression correlation with subtype

library(data.table)
library(ggplot2)

# load data
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/tf_subtype_correlation"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

# Setup files and data and variables
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_cluster_statistics/distal.prom.motif.res.RData")
distal.prom.motif.res <- as.data.table(distal.prom.motif.res)

#### Get subtype correlation from TCGA ####
tfs <- distal.prom.motif.res[,.(symbol = hgnc)]
load(file = "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/cgdsr/new.total.genes.df.log2.RData")
pam50 <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/Motherlode/dataset/UNC_PAM50.txt", header = TRUE, sep = "\t", check.names = TRUE)
pam50[,Sample.ID := gsub("-", ".", Sample.ID)]
tcga.all <- new.total.genes.df.log2[,.SD,.SDcols = c("patient_id", intersect(tfs$symbol, colnames(new.total.genes.df.log2)))]
tcga.all <- merge(tcga.all, pam50[,.(patient_id = Sample.ID, PAM50)])
tcga.all <- tcga.all[PAM50 %in% c("LumA", "LumB", "Her2", "Basal")]

rown <- 1
subtype.rown <- 1
tfs.subtype.res <- list()
dir.create("boxplots", showWarnings = FALSE)
#tfs <- tfs[symbol == "SIX1"]
toi <- "SIX1"
# for(rown in 1:nrow(tfs)){
#   toi <- tfs[rown,get("symbol")]
#   if(toi %in% colnames(tcga.all)){
#     subtypes <- as.data.table(t(combn(c("LumA", "LumB", "Her2", "Basal"), 2)))
#     for(subtype.rown in 1:nrow(subtypes)){
#       subtype1 <- subtypes[subtype.rown,get("V1")]
#       subtype2 <- subtypes[subtype.rown,get("V2")]
#       t.test.res <- t.test(tcga.all[PAM50 == subtype1,get(toi)], tcga.all[PAM50 == subtype2,get(toi)])
#       subtypes[subtype.rown,statistic := t.test.res$statistic]
#       subtypes[subtype.rown,p.value := t.test.res$p.value]
#     }
#     tfs.subtype.res[[toi]] <- subtypes
#     message(toi)
#     mp <- ggplot(tcga.all[,.(PAM50, gene = get(toi))], aes(PAM50, gene, fill = PAM50)) + geom_boxplot() + theme_bw()
#     ggsave(paste0(toi, ".pdf"), mp, width = 3, height = 3)
#   }
# }
# save(tfs.subtype.res, file = "tfs.subtype.res.RData")
load(file = "tfs.subtype.res.RData")

# Get best subtype per toi
toi <- "EHF"
tfs.subtype.specificity <- data.table(symbol = names(tfs.subtype.res), subtype = "",
                                      best.p.value = 1, best.statistic = 0, worst.subtype = "",
                                      d.subtype = "", d.best.p.value = 1, d.best.statistic = 0,
                                      d.worst.subtype = "")
# Most enriched
for(toi in names(tfs.subtype.res)){
  tmp <- tfs.subtype.res[[toi]]
  # All have to be significant first
  if(tmp[which.min(p.value),get("p.value")] < 0.05){
    best.statistic.value <- tmp[which.min(p.value),get("statistic")]
    best.p.value.value <- tmp[which.min(p.value),get("p.value")]
    if(best.statistic.value < 0){
      best.subtype.value <- tmp[which.min(p.value),get("V2")]
      worst.subtype.value <- tmp[which.min(p.value),get("V1")]
    } else {
      best.subtype.value <- tmp[which.min(p.value),get("V1")]
      worst.subtype.value <- tmp[which.min(p.value),get("V2")]
    }
    rown <- 1
    for(rown in 1:nrow(tmp)){
      if(is.na(best.subtype.value)){
        best.subtype.value <- NA
      } else {
        if(tmp[rown,get("V1")] == best.subtype.value){
          if(sign(tmp[rown,get("statistic")]) == 1 & tmp[rown,get("p.value")] < 0.05){
            best.subtype.value <- best.subtype.value
          } else {
            best.subtype.value <- NA
          }
        } else if(tmp[rown,get("V2")] == best.subtype.value){
          if(sign(tmp[rown,get("statistic")]) == -1 & tmp[rown,get("p.value")] < 0.05){
            best.subtype.value <- best.subtype.value
          } else {
            best.subtype.value <- NA
          }
        }
      }
    }
    if(is.na(best.subtype.value)){
      tfs.subtype.specificity[symbol == toi,subtype := NA]
      tfs.subtype.specificity[symbol == toi,best.p.value := NA]
      tfs.subtype.specificity[symbol == toi,best.statistic := NA]
      tfs.subtype.specificity[symbol == toi,worst.subtype := NA]
    } else {
      tfs.subtype.specificity[symbol == toi,subtype := best.subtype.value]
      tfs.subtype.specificity[symbol == toi,best.p.value := best.p.value.value]
      tfs.subtype.specificity[symbol == toi,best.statistic := abs(best.statistic.value)]
      tfs.subtype.specificity[symbol == toi,worst.subtype := worst.subtype.value]
    }
  }
}
# Most depleted
tmp <- tfs.subtype.res[["GATA3"]]
toi <- "GATA3"
for(toi in names(tfs.subtype.res)){
  tmp <- tfs.subtype.res[[toi]]
  # All have to be significant first
  if(tmp[which.min(p.value),get("p.value")] < 0.05){
    best.statistic.value <- tmp[which.min(p.value),get("statistic")]
    best.p.value.value <- tmp[which.min(p.value),get("p.value")]
    if(best.statistic.value > 0){
      best.subtype.value <- tmp[which.min(p.value),get("V2")]
      worst.subtype.value <- tmp[which.min(p.value),get("V1")]
    } else {
      best.subtype.value <- tmp[which.min(p.value),get("V1")]
      worst.subtype.value <- tmp[which.min(p.value),get("V2")]
    }
    rown <- 1
    for(rown in 1:nrow(tmp)){
      if(is.na(best.subtype.value)){
        best.subtype.value <- NA
      } else {
        if(tmp[rown,get("V1")] == best.subtype.value){
          if(sign(tmp[rown,get("statistic")]) == -1 & tmp[rown,get("p.value")] < 0.05){
            best.subtype.value <- best.subtype.value
          } else {
            best.subtype.value <- NA
          }
        } else if(tmp[rown,get("V2")] == best.subtype.value){
          if(sign(tmp[rown,get("statistic")]) == 1 & tmp[rown,get("p.value")] < 0.05){
            best.subtype.value <- best.subtype.value
          } else {
            best.subtype.value <- NA
          }
        }
      }
    }
    if(is.na(best.subtype.value)){
      tfs.subtype.specificity[symbol == toi,d.subtype := NA]
      tfs.subtype.specificity[symbol == toi,d.best.p.value := NA]
      tfs.subtype.specificity[symbol == toi,d.best.statistic := NA]
      tfs.subtype.specificity[symbol == toi,d.worst.subtype := NA]
    } else {
      tfs.subtype.specificity[symbol == toi,d.subtype := best.subtype.value]
      tfs.subtype.specificity[symbol == toi,d.best.p.value := best.p.value.value]
      tfs.subtype.specificity[symbol == toi,d.best.statistic := abs(best.statistic.value)]
      tfs.subtype.specificity[symbol == toi,d.worst.subtype := worst.subtype.value]
    }
  }
}
save(tfs.subtype.specificity, file = "tfs.subtype.specificity.RData")
fwrite(tfs.subtype.specificity, "tfs.subtype.specificity.txt", sep = "\t")

