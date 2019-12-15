# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(jetset)
library(reshape2)
library(survminer)
library(survival)
library(openxlsx)
library(VennDiagram)

# Get the jetset best probe per gene
jetset = function(symbol){
  return(as.character(jmap(chip = "hgu133a", symbol = symbol)))
}

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/tf_multivariate_update_revision"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

lines.dict <- data.table(line.lower = c("bom", "brm", "lm", "common"),
                         line.fixed = c("BoM", "BrM2", "LM2", "Met"),
                         tissue = c("Bone", "Brain", "Lung", "Bone|Brain|Lung"))

#load(file = "some.multivariate.res.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/GEO_datasets/metastasis_patient_datasets/morales/metDataset.Rdata")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/cgdsr/new.total.genes.df.log2.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/motif_allinter_cluster_statistics/distal.prom.motif.res.RData")

metDataset <- t(metDataset)#[1:10,1:10])
metDataset <- as.data.table(metDataset, keep.rownames = "patient")

clinical <- fread("/Users/wcai/Google_Drive/_Lab/Data/GEO_datasets/metastasis_patient_datasets/morales/annotations/perou_annotation.txt", sep = "\t", header = TRUE, check.names = TRUE)
clinical <- clinical[Cohort != "NKI295",c(.(patient = GEO.array), .SD), .SDcols = c(grep("", colnames(clinical)))]

# Keep some variables
clinical$Subtype
clinical$ER.clinical <- factor(clinical$ER.clinical, c("0", "1"))
clinical$PR.clinical <- factor(clinical$PR.clinical, c("0", "1"))
clinical$Her2.clinical <- factor(clinical$Her2.clinical, c("0", "1"))
clinical$Age
clinical$LN.status <- factor(clinical$LN.status, c("0", "1"))
#clinical$Dscore # remove because custom score
#clinical$Proliferation # remove because custom score
clinical$T.Stage <- factor(clinical$T.Stage, c("1", "2", "3", "4"))
clinical$Differention <- factor(clinical$Differention, c("1", "2", "3", "4"))
clinical$Chemo <- factor(clinical$Chemo, c("0", "1"))
#clinical$Hormone <- factor(clinical$Hormone, c("0", "1")) # remove because already have hormone info

# Only keep some variables
covariates <- c("Subtype", "ER.clinical", "PR.clinical", "Her2.clinical",
                "Age", "T.Stage", "Differention", "Chemo")
#covariates <- c("Subtype", "Age", "T.Stage", "Differention", "Chemo")
clinical.final <- clinical[complete.cases(clinical[,.SD,.SDcols = covariates]),.SD,
                           .SDcols = c("patient", "MFS", grep("relapse", colnames(clinical), value = TRUE),
                                       covariates)]
nrow(clinical.final)
colnames(clinical.final)[1] <- "patient"
i <- "TFAP2C"
relapse.type <- "Lung.relapse"
final.list <- list()
no_probe <- c()

distal.prom.motif.res <- as.data.table(distal.prom.motif.res)
goi <- distal.prom.motif.res$hgnc

nrow(clinical.final[Subtype == "Basal"])

# For subtype expression specificity
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_subtype_correlation/tfs.subtype.specificity.RData")

# Perform all multivariate cox
relapse.type <- "Lung.relapse"

#### Test RUNX2 and RARA ####
met.subset <- metDataset[,.(patient, RUNX2 = get(jetset("RUNX2")),
                            RARA = get(jetset("RARA")),
                            ETS2 = get(jetset("ETS2")))]
cor.test(met.subset[,get("RUNX2")], met.subset[,get("RARA")])
cor.test(met.subset[,get("RUNX2")], met.subset[,get("ETS2")])
met.subset <- merge(met.subset, clinical.final, by = "patient")
ggplot(met.subset[Subtype == "Basal"], aes(RUNX2, RARA, color = Subtype)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

#metDataset[,.SD,.SDcols = 1]

#### Test RUNX2 and RARA ####

for(relapse.type in c("Lung.relapse", "Brain.relapse", "Any.relapse")){
  relapse.list <- list()
  message(relapse.type)
  i <- "TFAP2C"
  for(i in goi){
    if(is.na(jetset(i))){
      no_probe <- c(no_probe, i)
    } else {
      dat <- merge(clinical.final, metDataset[,.SD, .SDcols = c("patient", jetset(i))], by = "patient", all.x = TRUE)
      colnames(dat)[ncol(dat)] <- i
      dat[, quant := cut(get(i), quantile(get(i), c(0,1/2,1)), labels = c("Low", "High"), include.lowest = TRUE)]
      res.cox <- coxph(as.formula(paste0("Surv(MFS,", relapse.type, ") ~ ","quant + ", paste0(covariates, collapse = "+"))), data = dat)
      res.cox.sum <- summary(res.cox)
      res.cox.sum.df <- signif(res.cox.sum$coefficients[,c(2,5)],3)
      colnames(res.cox.sum.df) <- c("HR", "pval")
      HR <- t(res.cox.sum.df)["HR",]
      names(HR) <- paste0(names(HR), ".HR")
      pval <- t(res.cox.sum.df)["pval",]
      names(pval) <- paste0(names(pval), ".pval")
      
      logrank <- survdiff(as.formula(paste0("Surv(MFS,", relapse.type, ") ~ ","quant")), data = dat)
      logrank.p <- 1 - pchisq(logrank$chisq, 1)
      logrank.HR = (logrank$obs[2]/logrank$exp[2])/(logrank$obs[1]/logrank$exp[1])
      names(logrank.p) <- "logrank.pval"
      names(logrank.HR) <- "logrank.HR"
      idx <- order(c(seq_along(pval), seq_along(HR)))
      relapse.list[[i]] <- c(signif(logrank.p,3), signif(logrank.HR,3), (c(pval,HR))[idx])
      #relapse.list[[i]] <- c(signif(logrank.p,3), signif(logrank.HR,3), pval, HR, signif(res.cox.sum$waldtest[3],3))
    }
  }
  final.list[[relapse.type]] <- data.table(do.call(rbind, relapse.list), keep.rownames = "symbol")
}

# Adjust p.value and add column for Motif+HR concordance+Subtype specificity
final.list.adjust <- final.list
for(relapse.type in c("Lung.relapse", "Brain.relapse", "Any.relapse")){
  tmp <- final.list[[relapse.type]]
  # for(pval in grep("\\.pval", colnames(tmp), value = TRUE)){
  #   tmp[, eval(paste0(pval, ".adj")) := p.adjust(get(pval), method = "BH")]
  # }
  tmp <- merge(tmp, distal.prom.motif.res, by.x = "symbol", by.y = "hgnc", all.x = TRUE)
  tmp <- merge(tmp, tfs.subtype.specificity[,.(symbol, enriched.subtype = subtype, depleted.subtype = d.subtype)], by = "symbol", all.x = TRUE)
  if(relapse.type == "Lung.relapse"){
    tmp[quantHigh.HR < 1 & logrank.HR < 1 & !is.na(lm.down), lm.concordant := TRUE]
    tmp[quantHigh.HR > 1 & logrank.HR > 1 & !is.na(lm.up), lm.concordant := TRUE]
  } else if(relapse.type == "Brain.relapse"){
    tmp[quantHigh.HR < 1 & logrank.HR < 1 & !is.na(brm.down), brm.concordant := TRUE]
    tmp[quantHigh.HR > 1 & logrank.HR > 1 & !is.na(brm.up), brm.concordant := TRUE]
  }
  final.list.adjust[[relapse.type]] <- tmp
}
test <- final.list.adjust[["Lung.relapse"]]
fwrite(final.list.adjust[["Lung.relapse"]], "RSAT.hichip.some.lung.multivariate.padj.txt", sep = "\t")
fwrite(final.list.adjust[["Brain.relapse"]], "RSAT.hichip.some.brain.multivariate.padj.txt", sep = "\t")
fwrite(final.list.adjust[["Any.relapse"]], "RSAT.hichip.some.any.multivariate.padj.txt", sep = "\t")

wb <- createWorkbook("RSAT.res")
for(type in names(final.list.adjust)){
  tmp <- final.list.adjust[[type]]
  addWorksheet(wb, sheetName = type, zoom = 150)
  writeData(wb, type, tmp)
}
saveWorkbook(wb, paste0("RSAT.hichip.some.multivariate.xlsx"), overwrite = TRUE)

RSAT.hichip.some.multivariate.res <- final.list.adjust
save(RSAT.hichip.some.multivariate.res, file = "RSAT.hichip.some.multivariate.res.RData")


load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/tf_multivariate_update_revision/RSAT.hichip.some.multivariate.res.RData")

#### Get Overlapping and Specific Survival, sig by either cox or log-rank ####
lung.concordant.sig <- RSAT.hichip.some.multivariate.res[["Lung.relapse"]][lm.concordant == TRUE & (quantHigh.pval < 0.05 | logrank.pval < 0.05),.(symbol, logrank.pval, logrank.HR, quantHigh.pval, quantHigh.HR, class, brm.down, brm.up, lm.down, lm.up, ensembl, brm_l2fc, brm_padj, lm_l2fc, lm_padj, lm.concordant)]
brain.concordant.sig <- RSAT.hichip.some.multivariate.res[["Brain.relapse"]][brm.concordant == TRUE & (quantHigh.pval < 0.05 | logrank.pval < 0.05),.(symbol, logrank.pval, logrank.HR, quantHigh.pval, quantHigh.HR, class, brm.down, brm.up, lm.down, lm.up, ensembl, brm_l2fc, brm_padj, lm_l2fc, lm_padj, brm.concordant)]
both.concordant.sig <- lung.concordant.sig[symbol %in% intersect(lung.concordant.sig$symbol, brain.concordant.sig$symbol),.(symbol, class, brm.down, brm.up, lm.down, lm.up, ensembl, brm_l2fc, brm_padj, lm_l2fc, lm_padj)]
wb <- createWorkbook("RSAT.concordant")
addWorksheet(wb, sheetName = "lung.concordant", zoom = 150)
writeData(wb, "lung.concordant", lung.concordant.sig)
addWorksheet(wb, sheetName = "brain.concordant", zoom = 150)
writeData(wb, "brain.concordant", brain.concordant.sig)
addWorksheet(wb, sheetName = "both.concordant", zoom = 150)
writeData(wb, "both.concordant", both.concordant.sig)
saveWorkbook(wb, paste0("RSAT.concordant.xlsx"), overwrite = TRUE)

fwrite(lung.concordant.sig, "lung.concordant.sig.txt", sep = "\t")
fwrite(brain.concordant.sig, "brain.concordant.sig.txt", sep = "\t")
fwrite(both.concordant.sig, "both.concordant.sig.txt", sep = "\t")
#### Get Overlapping and Specific Survival, sig by either cox or log-rank  ####

pdf(paste0("survival.venn.pdf"), width = 2, height = 2)
draw.pairwise.venn(area1 = length(lung.concordant.sig$symbol), area2 = length(brain.concordant.sig$symbol), cross.area = length(both.concordant.sig$symbol),
                   category = c("", ""))
dev.off()
