# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)

source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/cgdsr/new.total.genes.df.log2.RData")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/greenleaf_update"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))
# Remove Bone line from DICT

lines.dict <- lines.dict[line.lower != "bom",]
# Last Updated: 190110
# Author: Wesley Cai
# Purpose: To stratify TCGA patients by ATAC-seq and see if it correlates with survival
library(data.table)
library(ggplot2)
library(ggrepel)
#### Functions ####
survivalDict <- data.table(event.abbrev = c("OS", "DSS", "DFI", "PFI"),
                           event.full = c("Overall Survival", "Disease-specific Survival",
                                          "Disease-Free Interval", "Progression-Free Interval"))
lines.dict <- data.table(line.lower = c("bom", "brm", "lm"),
                         line.fixed = c("BoM", "BrM2", "LM2"))
# time <- "OS.time"
# event <- "OS"
# df <- tmp
# name <- paste0(mda.sig.name, "/", newdir, "/", type, "_tertile")
# time <- as.character(paste0(type, ".time"))
# event <- type
# name <- paste0(mda.sig.name, "/", newdir, "/", type, "_tertile")
# df <- tmp
# title <- plottitle
drawSurvPlot <- function(time, event, name, df, title = "", quant = "quant", groups = 5){
  df <- as.data.table(df)
  if(groups == 2){
    pal = c("#2166AC", "#B2182B")
  } else if(groups == 3){
    pal = c("#2166AC", "#A9A9A9", "#B2182B")
  } else if(groups == 5){
    pal = c("#CA0020", "#F4A582", "#F7F7F7", "#92C5DE", "#0571B0")
  }

  dir.create(dirname(name), showWarnings = FALSE, recursive = TRUE)
  form <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", quant))
  fit <- surv_fit(form, data = df)

  myplot <- ggsurvplot(fit, data = df, pval = TRUE, pval.method = TRUE,
                       title = title,
                       ylab = survivalDict[event.abbrev == event,get("event.full")],
                       xlab = "Time (Months)",
                       palette = pal,
                       ggtheme = theme_survminer() + theme(plot.title = element_text(hjust = 0.5)))
  myplot
  pdf(paste0(name, ".pdf"), width = 4.5, height = 4, onefile = FALSE)
  print(myplot)
  dev.off()
  
  cox <- coxph(form, data = df)
  cox.sum <- summary(cox)
  print(cox.sum)
  p.value <- signif(cox.sum$coefficients["quantHigh", "Pr(>|z|)"],4)
  hr1 <- signif(exp(cox.sum$coefficients["quantHigh", "coef"]),4)
  hr2 <- signif(exp(-cox.sum$coefficients["quantHigh", "coef"]),4)
  fwrite(data.table(comp = "quantHigh", p.value = p.value, hr1 = hr1, hr2 = hr2), paste0(name, ".coxph.txt"), sep = "\t")

  message(name)
}
drawPCAPlot <- function(df, name){
  pc <- prcomp(t(df))
  pc_out <- as.data.frame(pc$x)
  percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
  percentage <- paste(colnames(pc_out), paste0("(", as.character(percentage), "%", ")") )
  #pc_out$samples <- row.names(pc_out)#sub("_rep.*", "", row.names(pc_out))
  pc_out$samples <- sub("\\..*", "", row.names(pc_out))
  mp <- ggplot(pc_out, aes(PC1, PC2, color = samples)) +
    geom_point() +
    xlab(percentage[1]) +
    ylab(percentage[2]) +
    theme_bw() +
    ggtitle(label = "Principal Components with Greenleaf BRCA Peaks")
  ggsave(file.path("pca", paste0(gsub("\\.pdf", "", name), ".pdf")), mp, width = 5.5, height = 4)
}
# df <- tmp
# meta <- final.meta
drawColoredPCAPlot <- function(df, name, meta){
  
  pc <- prcomp(t(df))
  pc_out <- as.data.frame(pc$x)
  percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
  percentage <- paste(colnames(pc_out), paste0("(", as.character(percentage), "%", ")") )
  #pc_out$samples <- row.names(pc_out)#sub("_rep.*", "", row.names(pc_out))
  pc_out$sample <- sub("\\..*", "", row.names(pc_out))
  pc_out <- merge(pc_out, meta, by = "sample")
  mp <- ggplot(pc_out, aes(PC1, PC2)) +
    geom_point(aes(fill = cor), pch = 21, color = "black", size = 2) +
    xlab(percentage[1]) +
    ylab(percentage[2]) +
    theme_bw() +
    scale_fill_viridis_c() +
    geom_text_repel(aes(label = label), force = 10) +
    ggtitle(label = "Principal Components with Greenleaf BRCA Peaks")
  mp
  ggsave(name, mp, width = 5.5, height = 4)
}
#### Functions ####

library(data.table)
library(edgeR)
library(ggplot2)
library(preprocessCore)
library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)
library(survival)
library(survminer)
library(cgdsr)

#### Read in mda counts ####
mda.list <- list()
mda.files <- gsub("\\.txt", "", grep(".txt", dir(path = "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190110_greenleaf_survival/mda_counts/"), value = TRUE))
for(i in mda.files){
  mda.list[[i]] <- fread(file = paste0("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190110_greenleaf_survival/mda_counts/", i, ".txt"), data.table = FALSE, header = FALSE)[,1]
}
mda.matrix <- as.data.table(do.call(cbind, mda.list))
mda.matrix.cpm <- log2(cpm(mda.matrix+1))
peaks <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190110_greenleaf_survival/mda_counts/peaks.bed")
peaks[, peakid := paste(V1, V2, V3, sep = "-")]
mda.matrix.peaks <- cbind(peaks[,.(peakid)], mda.matrix)
mda.matrix.cpm <- as.data.frame(mda.matrix.cpm)
row.names(mda.matrix.cpm) <- peaks$peakid
#### Read in mda counts ####

#### Get differential accessibility for top diffacc peaks ####
# res <- list()
# mda.matrix.df <- as.data.frame(mda.matrix)
# row.names(mda.matrix.df) <- peaks$peakid
# for(comp.name in c("met", "lm", "brm", "bom")){
#   comp <- comp.name
#   if(comp.name == "met"){
#     comp <- c("lm", "brm", "bom")
#   }
#   comp.cols <- grep(paste(comp, collapse = "|"), names(mda.matrix.df))
#   par.cols <- grep("par", names(mda.matrix.df))
#   type <- rep("other", ncol(mda.matrix.df))
#   type[comp.cols] <- "met"
#   type[par.cols] <- "par"
# 
#   colData <- data.frame(row.names = names(mda.matrix.df), type = type)
#   dds <- DESeqDataSetFromMatrix(mda.matrix.df, colData = colData, design = ~ type)
#   dds <- DESeq(dds)
#   result.df <- results(dds, contrast = c("type", "met", "par"), tidy = TRUE)
#   result.df[,c("chr", "start", "end")] <- do.call(rbind, strsplit(result.df$row, split = "-"))
#   #fwrite(result.df[,c("chr", "start", "end", "row", "log2FoldChange", "padj")], "test.bed", sep = "\t", col.names = FALSE)
#   res[[comp.name]] <- as.data.table(result.df[,c("chr", "start", "end", "row", "log2FoldChange", "padj")])
# }
# for(i in names(res)){
#   res[[i]] <- res[[i]][,.(row, log2FoldChange, padj)]
# }
# save(res, file = "mda_deseq2.RData")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190110_greenleaf_survival/mda_deseq2.RData")
#### Get differential accessibility for top diffacc peaks ####

#### Loop through all ####
# Get gene expression profile
pam50 <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/Motherlode/dataset/UNC_PAM50.txt", header = TRUE, sep = "\t", check.names = TRUE)
pam50$Sample.ID <- gsub("-", ".", substr(pam50$Sample.ID, 1, 12))
pam50 <- pam50[,c("Sample.ID", "PAM50")]
pam50 <- pam50[!duplicated(Sample.ID)]
m <- as.data.table(pam50)
m[,bcr_patient_barcode := gsub("\\.", "-", substr(Sample.ID, 1, 12))]

mycgds <- CGDS("http://www.cbioportal.org/")
clinDat <- as.data.table(getClinicalData(mycgds, "brca_tcga_all"), keep.rownames = TRUE)
clinDat[, bcr_patient_barcode := substring(gsub("\\.", "-", rn), 1, 12)]

mda.diffacc.res <- as.data.table(do.call(cbind,res.greenleaf))
mda.diffacc.res <- mda.diffacc.res[,.(peakid = met.row,
                                      met.l2fc = met.log2FoldChange, met.padj,
                                      brm.l2fc = brm.log2FoldChange, brm.padj,
                                      lm.l2fc = lm.log2FoldChange, lm.padj)]
numPadj <- 5e-5


mda.sig.list <- list()
tmp.common <- mda.diffacc.res[which(brm.padj < numPadj & lm.padj < numPadj),]
tmp.common <- tmp.common[which(sign(brm.l2fc) == sign(lm.l2fc)),]
mda.sig.list[["common"]] <- tmp.common[,get("peakid")]
mda.sig.list[["lm_all"]] <- mda.diffacc.res[which(lm.padj < numPadj),get("peakid")]
mda.sig.list[["brm_all"]] <- mda.diffacc.res[which(brm.padj < numPadj),get("peakid")]

# Write out the signature peaks
dir.create("signature_peaks", showWarnings = FALSE)
# First write out the background peakset
tmp <- mda.diffacc.res
tmp[,c("chr", "start", "end") := as.data.table(do.call(rbind, strsplit(peakid, "-")))]
fwrite(tmp[,.(chr, start, end)], paste0("signature_peaks/all_peaks.bed"), sep = "\t", col.names = FALSE)
# Next write out each peakset

mda.sig.all.list <- list()
mda.sig.all.list[[paste0("padj", numPadj)]] <- list()
for(comp in names(mda.sig.list)){
  tmp <- mda.diffacc.res[peakid %in% mda.sig.list[[comp]]]
  tmp[,c("chr", "start", "end") := as.data.table(do.call(rbind, strsplit(peakid, "-")))]
  # Write out all
  fwrite(tmp, paste0("signature_peaks/", comp, ".padj", numPadj, ".txt"), sep = "\t")
  # Write out for GREAT
  if(comp == "common"){
    comp.line <- "lm"
  } else {
    comp.line <- gsub("_.*", "", comp)
  }
  # All peaks
  out.file <- paste0("signature_peaks/", comp, ".padj", numPadj, ".great_all")
  fwrite(tmp[,.(chr, start, end)], paste0(out.file, ".bed"), sep = "\t", col.names = FALSE)
  mp <- ggplot(tmp, aes(lm.l2fc, brm.l2fc)) +
    geom_point(size = 0.5) + theme_bw()
  ggsave(paste0(out.file, ".correlation.png"), mp, width = 4, height = 4)
  
  mda.sig.all.list[[paste0("padj", numPadj)]][[comp]] <- tmp
  # Up peaks
  out.file <- paste0("signature_peaks/", comp, ".padj", numPadj, ".great_up")
  fwrite(tmp[get(paste0(comp.line, ".l2fc")) > 0,.(chr, start, end)], paste0(out.file, ".bed"), sep = "\t", col.names = FALSE)
  mp <- ggplot(tmp[get(paste0(comp.line, ".l2fc")) > 0,], aes(lm.l2fc, brm.l2fc)) +
    geom_point(size = 0.5) + theme_bw()
  ggsave(paste0(out.file, ".correlation.png"), mp, width = 4, height = 4)
  # Down peaks
  out.file <- paste0("signature_peaks/", comp, ".padj", numPadj, ".great_down")
  fwrite(tmp[get(paste0(comp.line, ".l2fc")) < 0,.(chr, start, end)], paste0(out.file, ".bed"), sep = "\t", col.names = FALSE)
  mp <- ggplot(tmp[get(paste0(comp.line, ".l2fc")) < 0,], aes(lm.l2fc, brm.l2fc)) +
    geom_point(size = 0.5) + theme_bw()
  ggsave(paste0(out.file, ".correlation.png"), mp, width = 4, height = 4)
}
save(mda.sig.all.list, file = "mda.sig.all.list.RData")

#Prepare brca reads
brca <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/greenleaf_atac_tcga/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts/BRCA_raw_counts.txt")
brca[, peakid := paste(seqnames, start, end, sep = "-")]
brca.mda <- merge(brca, mda.matrix.peaks, by = "peakid")

#Prepare clinical annotation
brca.id <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/greenleaf_atac_tcga/TCGA_identifier_mapping.txt", data.table = FALSE)
row.names(brca.id) <- gsub("-", "_", brca.id$bam_prefix)

mda.sig.name <- "lm_all"
mda.sig.name <- "common"
for(mda.sig.name in names(mda.sig.list)){
  # Select data
  mda.sig <- mda.sig.list[[mda.sig.name]]
  brca.mda.subset <- brca.mda[peakid %in% mda.sig,]
  brca.mda.subset <- brca.mda.subset[,.SD,.SDcols = !grep("mda_bom", colnames(brca.mda.subset))]
  brca.mda.subset.ma <- brca.mda.subset[,.SD, .SDcols = 7:ncol(brca.mda.subset)]
  
  # Normalize
  brca.mda.cpm <- cpm(brca.mda.subset.ma, log = TRUE, prior.count = 5)
  brca.mda.cpm.norm <- as.data.frame(normalize.quantiles(brca.mda.cpm))
  colnames(brca.mda.cpm.norm) <- colnames(brca.mda.cpm)
  
  # Case IDs, take only first 15 with is the patient sample including tumor type
  colnames(brca.mda.cpm.norm)[1:141] <- substr(brca.id[colnames(brca.mda.cpm.norm)[1:141], "Case_ID"], 1, 15)
  colnames(brca.mda.cpm.norm)[142:147] <- gsub("_rep.*", "", colnames(brca.mda.cpm.norm)[142:147])
  
  # Average replicates
  brca.mda.cpm.norm.avg <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
    sapply(unique(names(brca.mda.cpm.norm)), # for each unique column name
           function(col) rowMeans(brca.mda.cpm.norm[names(brca.mda.cpm.norm) == col]) # calculate row means
    )
  )
  brca.mda.cpm.norm.avg <- as.data.table(cbind(peakid = brca.mda.subset$peakid, brca.mda.cpm.norm.avg))
  #fwrite(brca.mda.cpm.norm.avg, "tmp.txt", sep = "\t")
  # Average signal for mets
  brca.mda.cpm.norm.avg$mda_met <- rowMeans(brca.mda.cpm.norm.avg[,.(mda_brm, mda_lm)])
  #save(brca.mda.cpm.norm.avg, file = "brca.mda.cpm.norm.avg.RData")
  
  # Draw PCA
  tmp <- brca.mda.cpm.norm.avg[peakid %in% mda.sig, .SD, .SDcols = c(2:ncol(brca.mda.cpm.norm.avg))]
  colnames(tmp) <- sub("_rep.*", "", colnames(tmp))
  colnames(tmp) <- sub("TCGA.*", "TCGA", colnames(tmp))
  drawPCAPlot(tmp, paste0(mda.sig.name, "_pr_analysis_", numPadj))
  
  # Classify samples
  brca.sam <- brca.mda.cpm.norm.avg[,.SD, .SDcols = grep("TCGA", names(brca.mda.cpm.norm.avg))]
  mda.sam <- brca.mda.cpm.norm.avg[,.SD, .SDcols = grep("mda_", names(brca.mda.cpm.norm.avg))]
  
  df.surv <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/greenleaf_atac_tcga/brca.greenleaf.survival.txt")
  # Remove those with PFI.time == 0
  df.surv <- df.surv[PFI.time != 0,]
  # Make time into months
  df.surv[, PFI.time := signif(PFI.time/30,3)]
  
  correlations <- list()
  comp <- "mda_met"
  comp.list <- c("mda_lm", "mda_brm", "mda_met")
  #comp.list <- c("mda_met")
  for(comp in comp.list){
    newdir <- paste0(comp, "_padj_", numPadj)
    cor.par <- do.call(rbind, lapply(brca.sam, function(x) cor.test(x, mda.sam$mda_par)$estimate))
    cor.met <- do.call(rbind, lapply(brca.sam, function(x) cor.test(x, mda.sam[,get(comp)])$estimate))
    
    # Bigger more like met
    cor.par.met <- as.data.frame(cor.met - cor.par + 1)
    #hist(cor.par.met$cor)
    cor.par.met$bcr_patient_barcode <- substr(row.names(cor.par.met), 1, 12)
    
    final <- as.data.table(merge(cor.par.met, df.surv, by = "bcr_patient_barcode"))
    #remove duplicated patient samples
    final <- final[!duplicated(bcr_patient_barcode),]
    
    # Draw PCA and color by correlation
    tmp <- brca.mda.cpm.norm.avg[, .SD, .SDcols = c(2:ncol(brca.mda.cpm.norm.avg))]
    colnames(tmp) <- substr(colnames(tmp),1, 12)
    final.meta <- data.table(sample = colnames(tmp))
    final.meta <- merge(final.meta, final, by.x = "sample", by.y = "bcr_patient_barcode", all.x = TRUE)
    final.meta[is.na(cor),label := sample]
    dir.create(file.path(mda.sig.name, "pca"), showWarnings = FALSE, recursive = TRUE)
    drawColoredPCAPlot(tmp, file.path(mda.sig.name, "pca", paste0(comp,".colorPCA.pdf")), final.meta)
    
    fwrite(final, paste0(mda.sig.name, "_", comp, "_", numPadj, "_datatable.txt"), sep = "\t")
    correlations[[comp]] <- final
    type <- "PFI"
    for(type in c("PFI")){
      # for(type in c("OS", "DSS", "PFI", "DFI")){
      plottitle <- paste0(lines.dict[line.lower == gsub("mda_", "", comp), get("line.fixed")], " Score")

      tmp <- as.data.table(merge(cor.par.met, df.surv, by = "bcr_patient_barcode"))
      q <- quantile(tmp$cor, na.rm = TRUE, probs = c(0.25, 0.75))
      tmp[cor < q[[1]],quant := "Low"]
      tmp[cor >= q[[2]],quant := "High"]
      tmp$quant <- factor(tmp$quant, c("Low", "High"))
      drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, "_quartile_upper_lower"), tmp, plottitle)
      # 
      tmp <- as.data.table(merge(cor.par.met, df.surv, by = "bcr_patient_barcode"))
      q <- quantile(tmp$cor, na.rm = TRUE, probs = c(0.5, 1))
      tmp[cor < q[[1]],quant := "Low"]
      tmp[cor >= q[[1]],quant := "High"]
      tmp$quant <- factor(tmp$quant, c("Low", "High"))
      drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, "_twotile"), tmp, plottitle)
      
      # ER positivity
      tmp <- as.data.table(merge(cor.par.met, df.surv, by = "bcr_patient_barcode"))
      tmp.clin <- merge(tmp, clinDat[,.(bcr_patient_barcode, ER_STATUS_BY_IHC)], by = "bcr_patient_barcode")
      # Only keep those with "negative" and "positive" IHC
      tmp.clin <- merge(tmp.clin, new.total.genes.df.log2[,.(bcr_patient_barcode = gsub("\\.", "-", substr(patient_id,1,12)), ESR1),], by = "bcr_patient_barcode")
      q <- quantile(tmp.clin$cor, na.rm = TRUE, probs = c(0.5, 1))
      tmp.clin[cor < q[[1]],quant := "Low"]
      tmp.clin[cor >= q[[1]],quant := "High"]
      tmp.clin$quant <- factor(tmp.clin$quant, c("Low", "High"))
      tmp.clin <- tmp.clin[ER_STATUS_BY_IHC %in% c("Negative", "Positive"),]
      tmp.clin[PFI == 1, label := PFI.time]
      tmp.clin[,bcr_patient_barcode := gsub("-", "\\.", substr(bcr_patient_barcode, 1, 12))]

      tmp.clin.pam <- merge(tmp.clin, pam50, by.x = "bcr_patient_barcode", by.y = "Sample.ID")
      
      # Subtype
      for(subtype in c("LumA", "LumB", "Basal", "Her2")){
        drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, ".", subtype, "_twotile_ER_IHC_positive"), tmp.clin.pam[PAM50 %in% subtype], subtype)
      }
      
      drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, "_PAM50"), tmp.clin.pam, "PAM50", "PAM50")
      
      drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, "_twotile_ER_IHC_positive"), tmp.clin[ER_STATUS_BY_IHC == "Positive"], "IHC ER positive")
      drawSurvPlot(as.character(paste0(type, ".time")), type, paste0(mda.sig.name, "/", newdir, "/", type, "_twotile_ER_IHC_negative"), tmp.clin[ER_STATUS_BY_IHC == "Negative"], "IHC ER negative")
      cor.pval <- signif(cor.test(tmp.clin$ESR1, tmp.clin$cor, method = "spearman")$p.value,3)
      cor.estimate <- signif(cor.test(tmp.clin$ESR1, tmp.clin$cor, method = "spearman")$estimate,2)
      
      mp <- ggplot(tmp.clin, aes(ER_STATUS_BY_IHC, cor, fill = ER_STATUS_BY_IHC)) +
        geom_violin() +
        geom_jitter() +
        #annotate("text", x = 12, y = -0.05, label = paste0("r = ", cor.estimate, "\np = ", cor.pval)) +
        theme_bw()
      mp
      er.cor <- t.test(tmp.clin[ER_STATUS_BY_IHC == "Positive",get("cor")],
                       tmp.clin[ER_STATUS_BY_IHC == "Negative",get("cor")])
      er.cor$p.value
      er.cor <- wilcox.test(tmp.clin[ER_STATUS_BY_IHC == "Positive",get("cor")],
                            tmp.clin[ER_STATUS_BY_IHC == "Negative",get("cor")])
      er.cor$p.value
      
      ggsave(paste0(mda.sig.name, "/", newdir, "/ER_vs_cor_violin.png"), mp)
      
      mp <- ggplot(tmp.clin, aes(ESR1, cor, color = ER_STATUS_BY_IHC)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        annotate("text", x = 12, y = -0.05, label = paste0("r = ", cor.estimate, "\np = ", cor.pval)) +
        theme_bw()
      mp
      ggsave(paste0(mda.sig.name, "/", newdir, "/ER_vs_cor.png"), mp)
      mp <- ggplot(tmp.clin, aes(ER_STATUS_BY_IHC, ESR1, fill = ER_STATUS_BY_IHC)) +
        geom_boxplot() +
        geom_jitter() +
        theme_bw()
      mp
      ggsave(paste0(mda.sig.name, "/", newdir, "/ER_IHC.png"), mp)
      pval <- tryCatch(signif(t.test(tmp.clin[ER_STATUS_BY_IHC == "Positive" & quant == "Low",get("ESR1")],tmp.clin[ER_STATUS_BY_IHC == "Positive" & quant == "High",get("ESR1")])$p.value,3),
                       error = function(e) return(NA))
      mp <- ggplot(tmp.clin[ER_STATUS_BY_IHC == "Positive",], aes(quant, ESR1, fill = quant)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) +
        geom_point(data = subset(tmp.clin[ER_STATUS_BY_IHC == "Positive",], PFI == 0), color = "grey") +
        geom_point(data = subset(tmp.clin[ER_STATUS_BY_IHC == "Positive",], PFI == 1), color = "firebrick") +
        geom_text_repel(aes(label = label), force = 50) +
        labs(fill = "metATAC") +
        xlab(label = "") +
        annotate("text", x = 1.5, y = 5, label = paste0("p = ", pval)) +
        theme_bw()
      mp
      ggsave(paste0(mda.sig.name, "/", newdir, "/ER_IHC_Positive.png"), mp, height = 4, width = 5)
      pval <- tryCatch(signif(t.test(tmp.clin[ER_STATUS_BY_IHC == "Negative" & quant == "Low",get("ESR1")],tmp.clin[ER_STATUS_BY_IHC == "Negative" & quant == "High",get("ESR1")])$p.value,3),
                       error = function(e) return(NA))
      mp <- ggplot(tmp.clin[ER_STATUS_BY_IHC == "Negative",], aes(quant, ESR1, fill = quant)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5) +
        geom_point(data = subset(tmp.clin[ER_STATUS_BY_IHC == "Negative",], PFI == 0), color = "grey") +
        geom_point(data = subset(tmp.clin[ER_STATUS_BY_IHC == "Negative",], PFI == 1), color = "firebrick") +
        geom_text_repel(aes(label = label), force = 50) +
        labs(fill = "metATAC") +
        xlab(label = "") +
        annotate("text", x = 1.5, y = 5, label = paste0("p = ", pval)) +
        theme_bw()
      mp
      ggsave(paste0(mda.sig.name, "/", newdir, "/ER_IHC_Negative.png"), mp, height = 4, width = 5)

    }
  }
}

#### Loop through all ####
