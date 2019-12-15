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
gene <- "RARA"
relapse.type <- "Lung.relapse"
manualPlot <- function(gene, relapse.type, subtype = "All"){
  if(is.na(jetset(gene))){
    message("no probe")
  } else {
    i <- gene
    dat <- merge(clinical.final, metDataset[,.SD, .SDcols = c("patient", jetset(i))], by = "patient", all.x = TRUE)
    colnames(dat)[ncol(dat)] <- i
    dat[, quant := cut(get(i), quantile(get(i), c(0,1/2,1), na.rm = TRUE), labels = c("Low", "High"), include.lowest = TRUE)]
    dat$quant <- factor(dat$quant, c("Low", "High"))
    
    if(subtype == "All"){
      dat <- dat
    } else {
      dat <- dat[Subtype == subtype,]
    }
    drawSurvPlot(time = "MFS", event = relapse.type, name = file.path(paste0(i, "_", subtype, "_", relapse.type)), df = dat, title = paste(i, subtype, sep = " "), ylab = paste0(gsub("\\.relapse", "", relapse.type), " Metastasis Free (%)"))
  }
}

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/manual_kmplots_update"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#load(file = "some.multivariate.res.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/GEO_datasets/metastasis_patient_datasets/morales/metDataset.Rdata")
#load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/cgdsr/new.total.genes.df.log2.RData")
# Transpose metDataset
#test <- t(metDataset[1:10,1:10])
metDataset <- t(metDataset)#[1:10,1:10])
metDataset <- as.data.table(metDataset, keep.rownames = "patient")
head(metDataset)
clinical <- fread("/Users/wcai/Google_Drive/_Lab/Data/GEO_datasets/metastasis_patient_datasets/morales/annotations/perou_annotation.txt", sep = "\t", header = TRUE, check.names = TRUE)
clinical <- clinical[Cohort != "NKI295",c(.(patient = GEO.array), .SD), .SDcols = c(grep("", colnames(clinical)))]
#clinical <- clinical[,c(.(patient = GEO.array), .SD), .SDcols = c(grep("", colnames(clinical)))]
clinical.final <- clinical
covariates <- c("Subtype", "ER.clinical", "PR.clinical", "Her2.clinical",
               "Age", "LN.status", "T.Stage",
               "Differention", "Chemo")

# covariates <- c("Subtype",
#                 "Age", "LN.status", "T.Stage",
#                 "Differention", "Chemo")

covariates <- c("Subtype", "LN.status", "Age")
#clinical.final <- clinical[complete.cases(clinical[,.SD,.SDcols = covariates])]

#covariates <- c("Subtype", "Age", "T.Stage", "Differention", "Chemo")
clinical.final <- clinical[complete.cases(clinical[,.SD,.SDcols = covariates]),.SD,
                           .SDcols = c("patient", "MFS", grep("relapse", colnames(clinical), value = TRUE),
                                       covariates)]
nrow(clinical.final)
cohorts.remaining <- clinical[patient %in% clinical.final$patient]
unique(cohorts.remaining$Cohort)


gois <- c("SMAD1", "RFX7", "ELF4", "FOXM1", "CEBPB", "ATF3", "TEAD4", "NFE2L3", "BATF3", "NRBF2", "TEAD3", "FOXK2", "RUNX2", "JUN", "ELK3")
gois <- c("RARA", "TFAP2C", "ELF4", "APC2", "SSRP1", "RUNX1", "HSF2", "HMGN3", "TFAM", "RUNX2", "CEBPB", "CREB3L1", "MAFG", "GATA3", "ATF6", "TEAD3", "ATF6B", "TEAD4", "THRA", "HMGN1", "TEF", "HSF1", "HMGB3", "THRB")
gois <- c("TFAP2C", "RUNX2", "ELF4", "SMAD1", "NRBF2", "RARA", "JUN", "RFX7")

gois <- unique(read.xlsx("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update/additional_file_1.xlsx", sheet = 2)[,1],
               read.xlsx("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update/additional_file_1.xlsx", sheet = 2)[,1])
gois <- c("RUNX2")
for(goi in gois){
  for(relapse.type in c("Brain.relapse", "Lung.relapse", "Any.relapse", "Bone.relapse")){
    manualPlot(gene = goi, relapse.type = relapse.type)
  }
}

gois <- c("POU2F1", "POU2F2", "POU3F2", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "ATF4",
          "CEBPB", "CEBPD", "CEBPG", "DBP", "NFIL3", "TEF")
gois <- c("TFAP2C", "JUN", "RARA", "RFX7", "ELF4", "RUNX2", "ETS2")
motifs <- as.data.table(read.xlsx("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_cluster_statistics/compiled_all_formatted.xlsx", sheet = 1))
gois <- unique(motifs[direction %in% c("brm-up", "down-brm"), get("hgnc")])
for(goi in gois){
  for(relapse.type in c("Brain.relapse")){
    manualPlot(gene = goi, relapse.type = relapse.type)
    manualPlot(gene = goi, relapse.type = relapse.type, subtype = "Basal")
  }
}
gois <- unique(motifs[direction %in% c("down-shared", "up-shared"), get("hgnc")])
for(goi in gois){
  for(relapse.type in c("Lung.relapse", "Brain.relapse")){
    manualPlot(gene = goi, relapse.type = relapse.type)
    manualPlot(gene = goi, relapse.type = relapse.type, subtype = "Basal")
  }
}
gois <- unique(motifs[direction %in% c("lm-up", "down-lm"), get("hgnc")])
for(goi in gois){
  for(relapse.type in c("Lung.relapse")){
    manualPlot(gene = goi, relapse.type = relapse.type)
    manualPlot(gene = goi, relapse.type = relapse.type, subtype = "Basal")
  }
}

manualPlot(gene = "TFAP2C", relapse.type = "Lung.relapse", subtype = "Basal")
manualPlot(gene = "JUN", relapse.type = "Brain.relapse", subtype = "Basal")
manualPlot(gene = "RARA", relapse.type = "Lung.relapse", subtype = "Basal")
manualPlot(gene = "RFX7", relapse.type = "Brain.relapse", subtype = "Basal")
manualPlot(gene = "ELF4", relapse.type = "Lung.relapse", subtype = "Basal")
manualPlot(gene = "ELF4", relapse.type = "Brain.relapse", subtype = "Basal")
manualPlot(gene = "RUNX2", relapse.type = "Lung.relapse", subtype = "Basal")
manualPlot(gene = "RUNX2", relapse.type = "Brain.relapse", subtype = "Basal")
manualPlot(gene = "ETS2", relapse.type = "Lung.relapse", subtype = "Basal")

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
