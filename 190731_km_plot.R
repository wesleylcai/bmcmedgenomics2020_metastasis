# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(jetset)
library(reshape2)
library(survminer)
library(survival)
library(RColorBrewer)
library(openxlsx)

# Get the jetset best probe per gene
jetset = function(symbol){
  return(as.character(jmap(chip = "hgu133a", symbol = symbol)))
}

drawSurvPlot <- function(time, event, name, df, title = "", ylab = ""){
  quant.length <- length(na.omit(unique(df$quant)))
  if(quant.length == 2){
    pal = c("#2166AC", "#B2182B")
  } else if(quant.length == 3){
    pal = c("#2166AC", "#A9A9A9", "#B2182B")
  } else if(quant.length == 4){
    pal = c("#2166AC", "#92C5DE", "#F4A582", "#B2182B")
  }
  dir.create(dirname(name), showWarnings = FALSE, recursive = TRUE)
  form <- as.formula(paste0("Surv(", time, ",", event, ") ~ quant"))
  fit <- surv_fit(form, data = df)
  myplot <- ggsurvplot(fit, data = df, pval = TRUE, pval.method = TRUE,
                       title = title,
                       ylab = ylab,
                       xlab = "Time (Months)",
                       palette = pal,
                       ggtheme = theme_survminer() + theme(plot.title = element_text(hjust = 0.5))) +
    guides(colour = guide_legend(nrow = quant.length))

  pdf(paste0(name, ".pdf"), width = 4, height = 3.75+quant.length/8, onefile = FALSE)
  print(myplot)
  dev.off()
  message(name)
}
i <- "LMO4"
gene <- "TFAP2C"
relapse.type <- "Lung.relapse"
manualPlot <- function(gene, relapse.type){
  if(is.na(jetset(gene))){
    message("no probe")
  } else {
    i <- gene
    dat <- merge(clinical.final, metDataset[,.SD, .SDcols = c("patient", jetset(i))], by = "patient", all.x = TRUE)
    colnames(dat)[ncol(dat)] <- i
    dat[, quant := cut(get(i), quantile(get(i), c(0,1/2,1), na.rm = TRUE), labels = c("Low", "High"), include.lowest = TRUE)]
    dat$quant <- factor(dat$quant, c("Low", "High"))
    drawSurvPlot(time = "MFS", event = relapse.type, name = file.path(paste0(i, "_", relapse.type)), df = dat, title = paste0(i), ylab = paste0(gsub("\\.relapse", "", relapse.type), " Metastasis Free (%)"))
  }
}
gene1 <- "TFAP2C"
gene2 <- "GATA3"
relapse.type <- "Lung.relapse"
manualMultiPlot <- function(gene1, gene2, relapse.type){
  gene1.low <- paste0(gene1, " Low")
  gene1.high <- paste0(gene1, " High")
  gene2.low <- paste0(gene2, " Low")
  gene2.high <- paste0(gene2, " High")
  dat <- merge(clinical.final, metDataset[,.SD, .SDcols = c("patient", jetset(gene1), jetset(gene2))], by = "patient", all.x = TRUE)
  colnames(dat)[(ncol(dat)-1):ncol(dat)] <- c(gene1, gene2)
  dat[, quant1 := cut(get(gene1), quantile(get(gene1), c(0,1/2,1)), labels = c(gene1.low, gene1.high), include.lowest = TRUE)]
  dat[, quant2 := cut(get(gene2), quantile(get(gene2), c(0,1/2,1)), labels = c(gene2.low, gene2.high), include.lowest = TRUE)]
  dat[, quant := paste0(quant1, "; ", quant2)]
  dat$quant <- factor(dat$quant, c(paste0(gene1.low, "; ", gene2.high),
                                   paste0(gene1.low, "; ", gene2.low),
                                   paste0(gene1.high, "; ", gene2.high),
                                   paste0(gene1.high, "; ", gene2.low)))
  drawSurvPlot(time = "MFS", event = relapse.type, name = file.path(paste0(gene1, "_", gene2, "_", relapse.type)), df = dat, title = "", ylab = paste0(gsub("\\.relapse", "", relapse.type), " Metastasis Free (%)"))
}

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/manual_kmplots"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#load(file = "some.multivariate.res.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/GEO datasets/metastasis_patient_datasets/morales/metDataset.Rdata")
#load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/cgdsr/new.total.genes.df.log2.RData")
# Transpose metDataset
#test <- t(metDataset[1:10,1:10])
metDataset <- t(metDataset)#[1:10,1:10])
metDataset <- as.data.table(metDataset, keep.rownames = "patient")

clinical <- fread("/Users/wcai/Google_Drive/_Lab/Data/GEO datasets/metastasis_patient_datasets/morales/annotations/perou_annotation.txt", sep = "\t", header = TRUE, check.names = TRUE)
clinical <- clinical[Cohort != "NKI295",c(.(patient = GEO.array), .SD), .SDcols = c(grep("", colnames(clinical)))]
#clinical <- clinical[,c(.(patient = GEO.array), .SD), .SDcols = c(grep("", colnames(clinical)))]
covariates <- c("ER.clinical", "PR.clinical", "Her2.clinical",
               "Age", "LN.status", "T.Stage",
               "Differention", "Chemo")
#clinical.final <- clinical[complete.cases(clinical[,.SD,.SDcols = covariates])]

#covariates <- c("Subtype", "Age", "T.Stage", "Differention", "Chemo")
clinical.final <- clinical[complete.cases(clinical[,.SD,.SDcols = covariates]),.SD,
                           .SDcols = c("patient", "MFS", grep("relapse", colnames(clinical), value = TRUE),
                                       covariates)]
nrow(clinical.final)

gois <- c("SMAD1", "RFX7", "ELF4", "FOXM1", "CEBPB", "ATF3", "TEAD4", "NFE2L3", "BATF3", "NRBF2", "TEAD3", "FOXK2", "RUNX2", "JUN", "ELK3")
gois <- c("RARA", "TFAP2C", "ELF4", "APC2", "SSRP1", "RUNX1", "HSF2", "HMGN3", "TFAM", "RUNX2", "CEBPB", "CREB3L1", "MAFG", "GATA3", "ATF6", "TEAD3", "ATF6B", "TEAD4", "THRA", "HMGN1", "TEF", "HSF1", "HMGB3", "THRB")
gois <- c("NRBF2")

gois <- unique(read.xlsx("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update/additional_file_1.xlsx", sheet = 2)[,1],
               read.xlsx("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update/additional_file_1.xlsx", sheet = 2)[,1])
gois <- c("RUNX2")
for(goi in gois){
  for(relapse.type in c("Brain.relapse", "Lung.relapse", "Any.relapse", "Bone.relapse")){
    manualPlot(gene = goi, relapse.type = relapse.type)
  }
}

# manualPlot(gene = "NRBF2", relapse.type = "Brain.relapse")
# manualPlot(gene = "LMO2", relapse.type = "Any.relapse")
# manualPlot(gene = "LMO2", relapse.type = "Lung.relapse")
# manualPlot(gene = "LMO2", relapse.type = "Brain.relapse")
# manualPlot(gene = "ETS1", relapse.type = "Any.relapse")
# manualPlot(gene = "ETS1", relapse.type = "Lung.relapse")
# manualPlot(gene = "ETS1", relapse.type = "Brain.relapse")
# manualPlot(gene = "RARA", relapse.type = "Lung.relapse")
# manualPlot(gene = "TFAP2C", relapse.type = "Lung.relapse")
# manualPlot(gene = "RUNX2", relapse.type = "Lung.relapse")
# manualPlot(gene = "GATA3", relapse.type = "Lung.relapse")
# manualPlot(gene = "JUN", relapse.type = "Brain.relapse")
# manualPlot(gene = "RUNX2", relapse.type = "Brain.relapse")
# manualPlot(gene = "HIF1A", relapse.type = "Brain.relapse")
# manualPlot(gene = "TEAD4", relapse.type = "Brain.relapse")
# manualPlot(gene = "SIX1", relapse.type = "Brain.relapse")
# manualMultiPlot(gene1 = "TFAP2C", gene2 = "GATA3", relapse.type = "Lung.relapse")
# manualMultiPlot(gene1 = "JUN", gene2 = "RUNX2", relapse.type = "Brain.relapse")
# 
# goi.list <- list()
# goi.list[["lm.down.specific"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_DBD_statistics/expfilter_DBD_down.lm.specific.txt", header = FALSE)
# goi.list[["brm.down.specific"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_DBD_statistics/expfilter_DBD_down.brm.specific.txt", header = FALSE)
# goi.list[["lm.up.specific"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_DBD_statistics/expfilter_DBD_up.lm.specific.txt", header = FALSE)
# goi.list[["brm.up.specific"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_DBD_statistics/expfilter_DBD_up.brm.specific.txt", header = FALSE)
# goi.list[["up.shared"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_DBD_statistics/expfilter_DBD_up.shared.txt", header = FALSE)
# for(i in names(goi.list)){
#   for(goi in goi.list[[i]]$V1){
#     for(relapse.type in c("Brain.relapse", "Lung.relapse", "Any.relapse")){
#       manualPlot(gene = goi, relapse.type = relapse.type)
#     }
#   }
# }
# 
# 
# load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update/natmeth.hichip.some.multivariate.res.RData")
# 
# #### Get Overlapping and Specific Survival, sig by either cox or log-rank ####
# test <- natmeth.hichip.some.multivariate.res$Brain.relapse
# lung.concordant.sig <- natmeth.hichip.some.multivariate.res[["Lung.relapse"]][lm.concordant == TRUE & (quantHigh.pval < 0.05 | pval.logrank < 0.05),get("symbol")]
# brain.concordant.sig <- natmeth.hichip.some.multivariate.res[["Brain.relapse"]][brm.concordant == TRUE & (quantHigh.pval < 0.05 | pval.logrank < 0.05),get("symbol")]
# intersect(lung.concordant.sig, brain.concordant.sig)
# gois <- c(lung.concordant.sig, brain.concordant.sig)
# for(goi in gois){
#   for(relapse.type in c("Brain.relapse", "Lung.relapse", "Any.relapse")){
#     manualPlot(gene = goi, relapse.type = relapse.type)
#   }
# }
# #### Get Overlapping and Specific Survival, sig by either cox or log-rank  ####
