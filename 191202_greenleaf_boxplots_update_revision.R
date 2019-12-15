# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/greenleaf_subtype_boxplots_update"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

lines.dict <- data.table(line.lower = c("bom", "brm", "lm", "met"),
                         line.fixed = c("BoM", "BrM2", "LM2", "Common"))

score.list <- list()
score.list[["lm_5e-05"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/greenleaf_update/lm_all_mda_lm_5e-05_datatable.txt")
score.list[["brm_5e-05"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/greenleaf_update/brm_all_mda_brm_5e-05_datatable.txt")
score.list[["met_5e-05"]] <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/greenleaf_update/common_mda_met_5e-05_datatable.txt")

# Get gene expression profile
library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/")
test(mycgds)
study <- cgdsr::getCancerStudies(mycgds)
mycaselist <- getCaseLists(mycgds, "brca_tcga")
cgdsr.clin <- as.data.table(getClinicalData(mycgds, caseList = "brca_tcga_all"), keep.rownames = "Sample.ID")
cgdsr.clin.er <- cgdsr.clin[,.(Sample.ID, ER_STATUS_BY_IHC)]
length(intersect(cgdsr.clin.er$Sample.ID, pam50$Sample.ID))

pam50 <- read.table("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/Motherlode/dataset/UNC_PAM50.txt", header = TRUE, sep = "\t")
pam50$Sample.ID <- gsub("-", ".", pam50$Sample.ID)
pam50 <- pam50[,c("Sample.ID", "PAM50")]
pam50 <- merge(pam50, cgdsr.clin.er, by = "Sample.ID", all.x = TRUE)
clinDat <- as.data.table(pam50)
clinDat[,bcr_patient_barcode := gsub("\\.", "-", substr(Sample.ID, 1, 12))]

# Boxplots for subtype
line <- "lm_5e-05"
line <- "met_5e-05"
for(line in names(score.list)){
  scores <- score.list[[line]]
  scores <- merge(scores, clinDat[,.(bcr_patient_barcode, PAM50, ER_STATUS_BY_IHC)])
  scores <- scores[PFI.time != 0,]
  xlabel <- paste0(lines.dict[line.lower == gsub("_.*", "", line), get("line.fixed")], " metATAC score")
  scores <- scores[PAM50 %in% c("Basal", "Her2", "LumA", "LumB"),]
  mp <- ggplot(scores, aes_string("PAM50", "cor", fill = "PAM50")) +
    geom_violin(show.legend = FALSE) +
    geom_jitter(show.legend = FALSE) +
    #geom_smooth(method = "lm", se = FALSE) +
    #annotate("text", label = paste0("r2 = ", rsquared, "\n", "p = ", pval),
    #          x = annote.x, y = annote.y, hjust = 0, size = 5) +
    xlab(label = "") +
    ylab(label = xlabel) +
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 14, color = "black"))
  mp
  #ggsave("test.pdf", mp, width = 5.5, height = 5)
  ggsave(paste0(gsub("_.*", "", line), ".subtype.metATAC.png"), mp, width = 4.5, height = 4)
  
  er.scores <- scores[ER_STATUS_BY_IHC != ""]
  mp <- ggplot(er.scores, aes_string("ER_STATUS_BY_IHC", "cor", fill = "ER_STATUS_BY_IHC")) +
    geom_violin() +
    geom_jitter(show.legend = FALSE) +
    #geom_smooth(method = "lm", se = FALSE) +
    #annotate("text", label = paste0("r2 = ", rsquared, "\n", "p = ", pval),
    #          x = annote.x, y = annote.y, hjust = 0, size = 5) +
    xlab(label = "") +
    ylab(label = xlabel) +
    theme_bw() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 14, color = "black"))
  mp
  ggsave(paste0(gsub("_.*", "", line), ".ERstatus.metATAC.png"), mp, width = 4.5, height = 4)
  
  wilcox.test(er.scores[ER_STATUS_BY_IHC == "Positive",get("cor")],
              er.scores[ER_STATUS_BY_IHC == "Negative",get("cor")], paired = FALSE)
  
  res <- as.data.table(t(combn(unique(scores$PAM50), 2)))
  colnames(res) <- c("subtype1", "subtype2")
  for(i in 1:nrow(res)){
    subtype1.comp <- res[i, get("subtype1")]
    subtype2.comp <- res[i, get("subtype2")]
    res[i, t.test := tryCatch(signif(t.test(scores[PAM50 == subtype1.comp, get("cor")], scores[PAM50 == subtype2.comp, get("cor")])$p.value, 3), error = function(e){return(NA)})]
    res[i, wilcox := tryCatch(signif(wilcox.test(scores[PAM50 == subtype1.comp, get("cor")], scores[PAM50 == subtype2.comp, get("cor")])$p.value, 3), error = function(e){return(NA)})]
  }
  write.table(res, paste0(gsub("_.*", "", line), ".subtype.pvals.txt"), sep = "\t")
  
}
