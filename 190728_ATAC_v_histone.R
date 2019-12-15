# Last Updated: 
# Author: Wesley Cai
# Purpose: Plot histone vs ATAC

library(data.table)
library(ggplot2)

# load data
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/joint"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

annotations.overlap <- grep("overlap", names(chip.annotations), value = TRUE)
i <- annotations.overlap[[1]]
for(i in annotations.overlap){
  overlap.names <- strsplit(i, "\\.")[[1]][1:2]
  histone <- paste0(overlap.names[2], ".broad")
  overlap <- chip.annotations[[i]]
  overlap.l2fc <- merge(overlap, res.list[["atac"]][,.(peakid, atac.l2fc = met_l2fc, atac.padj = met_padj)], by.x = "atac.peakid", by.y = "peakid", all.x = TRUE)
  overlap.l2fc <- merge(overlap.l2fc, res.list[[histone]][,.(peakid, histone.l2fc = met_l2fc, histone.padj = met_padj)], by.x = paste0(overlap.names[2], ".peakid"), by.y = "peakid", all.x = TRUE)
  overlap.l2fc <- overlap.l2fc[complete.cases(overlap.l2fc),]
  overlap.l2fc.sig <- overlap.l2fc[atac.padj < 0.05 & histone.padj < 0.05,]
  ggplot(overlap.l2fc.sig, aes(atac.l2fc, histone.l2fc)) +
    geom_point()
}


#### Plot correlations ####
lines <- c("brm", "lm", "met")

i <- "atac.h3k27ac.broad.overlap"
annote.int <- data.table(annote = c("atac.h3k27ac.broad.overlap", "atac.h3k4me3.broad.overlap", "atac.h3k4me1.broad.overlap"),
                         histone = c("h3k27ac", "h3k4me3", "h3k4me1"),
                         title = c("H3K27ac", "H3K4me3", "H3K4me1"))
i <- 1
line <- "met"
for(i in 1:nrow(annote.int)){
  for(line in lines){
    dir.create(file.path(getwd(), line), showWarnings = FALSE, recursive = TRUE)
    annote <- annote.int[i,get("annote")]
    histone <- annote.int[i,get("histone")]
    title <- annote.int[i,get("title")]
    annote.peaks <- chip.annotations[[annote]]
    sig.atac <- res.list[["atac"]][get(paste0(line, "_padj")) < 0.05,]
    sig.histone <- res.list[[paste0(histone, ".broad")]][get(paste0(line, "_padj")) < 0.05,]
    tmp.cor <- merge(annote.peaks, sig.atac[,.(peakid, brm.atac = brm_l2fc, lm.atac = lm_l2fc, met.atac = met_l2fc)], by.x = "atac.peakid", by.y = "peakid")
    tmp.cor <- merge(tmp.cor, sig.histone[,.(peakid, brm.histone = brm_l2fc, lm.histone = lm_l2fc, met.histone = met_l2fc)], by.x = paste0(histone, ".peakid"), by.y = "peakid")
    #tmp.cor$Density <- get_density(tmp.cor[,get(paste0(line, ".atac"))], tmp.cor[,get(paste0(line, ".histone"))])
    cor <- cor.test(tmp.cor[,get(paste0(line, ".atac"))], tmp.cor[,get(paste0(line, ".histone"))])
    r2 <- signif(cor$estimate, 3)
    pval <- signif(cor$p.value, 3)
    if(pval < 2.2e-16){
      pval <- paste0("p < 2.2e-16")
    } else {
      pval <- paste0("p = ", pval)
    }
    mp <- ggplot(tmp.cor, aes_string(paste0(line, ".atac"), paste0(line, ".histone"))) +#, color = "Density")) +
      geom_point(size = 0.5) + 
      #scale_color_viridis() +
      geom_smooth(method = "lm", col = "darkgray", linetype = "dashed") +
      theme_bw() +
      ylab(label = bquote("ATAC "*log[2]*" Fold Change")) +
      xlab(label = bquote(log[2]*" Fold Change")) +
      ggtitle(label = title) +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("       r2 = ", r2, "\n", "       ", pval),
               hjust = 0, vjust = 2, size = 6) +
      theme(axis.title = element_text(size = 20),
            axis.text = element_text(size = 14),
            plot.title = element_text(size = 30, hjust = 0.5),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"))
    mp
    ggsave(file.path(line, paste0(line, ".atac-v-", histone, ".annote-", annote, ".png")), units = "in", dpi = 300, width = 4.4, height = 6)#, width = 5, height = 6)
  }
}

