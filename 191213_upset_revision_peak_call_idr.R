library(ChIPpeakAnno)
library(GenomicRanges)

#### IDR Version ####
for(mark in c("H3K4me3", "H3K27ac")){
  BrM1 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_BrM_", mark, "_R1_peaks.broadPeak"), format = "broadPeak")
  BrM2 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_BrM_", mark, "_R2_peaks.broadPeak"), format = "broadPeak")
  BrM.idr <- IDRfilter(BrM1, BrM2,
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_BrM_", mark, "_R1.mqsd.bam"),
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_BrM_", mark, "_R2.mqsd.bam"))
  
  LM1 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_LM_", mark, "_R1_peaks.broadPeak"), format = "broadPeak")
  LM2 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_LM_", mark, "_R2_peaks.broadPeak"), format = "broadPeak")
  LM.idr <- IDRfilter(LM1, LM2,
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_LM_", mark, "_R1.mqsd.bam"),
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_LM_", mark, "_R2.mqsd.bam"))
  
  Par1 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_Par_", mark, "_R1_peaks.broadPeak"), format = "broadPeak")
  Par2 <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/single_new/MDAMB231_Par_", mark, "_R2_peaks.broadPeak"), format = "broadPeak")
  Par.idr <- IDRfilter(Par1, Par2,
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_Par_", mark, "_R1.mqsd.bam"),
                       paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/bam/bamfiles/MDAMB231_Par_", mark, "_R2.mqsd.bam"))
  
  consens <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/merge/MDAMB231_Par_BrM_LM_", mark, "_merge_peaks.broadPeak"), format = "broadPeak")
  
  hits <- findOverlaps(consens, BrM.idr)
  overlaps <- pintersect(consens[queryHits(hits)], BrM.idr[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  BrM.idr.consens <- consens[queryHits(hits)[percentOverlap > 0.5]]
  BrM.idr.consens <- as.data.table(BrM.idr.consens)[,.(seqnames, start, end)]
  BrM.idr.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  BrM.idr.consens[, brm := 1]
  
  hits <- findOverlaps(consens, LM.idr)
  overlaps <- pintersect(consens[queryHits(hits)], LM.idr[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  LM.idr.consens <- consens[queryHits(hits)[percentOverlap > 0.5]]
  LM.idr.consens <- as.data.table(LM.idr.consens)[,.(seqnames, start, end)]
  LM.idr.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  LM.idr.consens[, lm := 1]
  
  hits <- findOverlaps(consens, Par.idr)
  overlaps <- pintersect(consens[queryHits(hits)], Par.idr[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  Par.idr.consens <- consens[queryHits(hits)[percentOverlap > 0.5]]
  Par.idr.consens <- as.data.table(Par.idr.consens)[,.(seqnames, start, end)]
  Par.idr.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  Par.idr.consens[, par := 1]
  
  consens.dt <- as.data.table(consens)[,.(seqnames, start, end)]
  consens.dt[, peakid := paste(seqnames, start, end, sep = "-")]
  
  consens.dt <- merge(consens.dt[,.(peakid)], BrM.idr.consens[,.(peakid, brm)], by = "peakid", all = TRUE)
  consens.dt <- merge(consens.dt, LM.idr.consens[,.(peakid, lm)], by = "peakid", all = TRUE)
  consens.dt <- merge(consens.dt, Par.idr.consens[,.(peakid, par)], by = "peakid", all = TRUE)
  
  consens.dt <- cbind(consens.dt[,.(peakid)],
        apply(consens.dt[,.(par,brm,lm)], c(1,2), function(x){if(is.na(x)){return(0)}else{return(x)}}))
  fwrite(consens.dt, paste0(mark, ".peaks.txt"), sep = "\t")
}
#### IDR Version ####

#### Regular Version ####
for(mark in c("H3K4me3", "H3K27ac")){
  BrM <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/merge_new/MDAMB231_BrM_", mark, "_merge_peaks.broadPeak"), format = "broadPeak")
  LM <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/merge_new/MDAMB231_LM_", mark, "_merge_peaks.broadPeak"), format = "broadPeak")
  Par <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/merge_new/MDAMB231_Par_", mark, "_merge_peaks.broadPeak"), format = "broadPeak")
 
  consens <- toGRanges(paste0("/gpfs/ysm/project/yan_qin/wc376/chip/MDAMB231/peaks/macs2/merge/MDAMB231_Par_BrM_LM_", mark, "_merge_peaks.broadPeak"), format = "broadPeak")
  
  hits <- findOverlaps(consens, BrM)
  overlaps <- pintersect(consens[queryHits(hits)], BrM[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  BrM.consens <- consens[queryHits(hits)[percentOverlap > 0]]
  BrM.consens <- as.data.table(BrM.consens)[,.(seqnames, start, end)]
  BrM.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  BrM.consens[, brm := 1]
  
  hits <- findOverlaps(consens, LM)
  overlaps <- pintersect(consens[queryHits(hits)], LM[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  LM.consens <- consens[queryHits(hits)[percentOverlap > 0]]
  LM.consens <- as.data.table(LM.consens)[,.(seqnames, start, end)]
  LM.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  LM.consens[, lm := 1]
  
  hits <- findOverlaps(consens, Par)
  overlaps <- pintersect(consens[queryHits(hits)], Par[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
  Par.consens <- consens[queryHits(hits)[percentOverlap > 0]]
  Par.consens <- as.data.table(Par.consens)[,.(seqnames, start, end)]
  Par.consens[, peakid := paste(seqnames, start, end, sep = "-")]
  Par.consens[, par := 1]
  
  consens.dt <- as.data.table(consens)[,.(seqnames, start, end)]
  consens.dt[, peakid := paste(seqnames, start, end, sep = "-")]
  
  consens.dt <- merge(consens.dt[!(duplicated(peakid)),.(peakid)], BrM.consens[!(duplicated(peakid)),.(peakid, brm)], by = "peakid", all = TRUE)
  consens.dt <- merge(consens.dt, LM.consens[!(duplicated(peakid)),.(peakid, lm)], by = "peakid", all = TRUE)
  consens.dt <- merge(consens.dt, Par.consens[!(duplicated(peakid)),.(peakid, par)], by = "peakid", all = TRUE)
  
  consens.dt <- cbind(consens.dt[,.(peakid)],
                      apply(consens.dt[,.(par,brm,lm)], c(1,2), function(x){if(is.na(x)){return(0)}else{return(x)}}))
  fwrite(consens.dt, paste0(mark, ".peaks.txt"), sep = "\t")
}

BrM <- toGRanges("/gpfs/ysm/project/yan_qin/wc376/atac/peaks/macs2/merge/daughtery/sharp/MDAMB231_BrM_atac.merge.sharp.daughtery_peaks.narrowPeak", format = "narrowPeak")
LM <- toGRanges("/gpfs/ysm/project/yan_qin/wc376/atac/peaks/macs2/merge/daughtery/sharp/MDAMB231_LM_atac.merge.sharp.daughtery_peaks.narrowPeak", format = "narrowPeak")
Par <- toGRanges("/gpfs/ysm/project/yan_qin/wc376/atac/peaks/macs2/merge/daughtery/sharp/MDAMB231_Par_atac.merge.sharp.daughtery_peaks.narrowPeak", format = "narrowPeak")

consens <- toGRanges("/gpfs/ysm/project/yan_qin/wc376/atac/peaks/macs2/multimerge/daughtery/sharp/MDAMB231_Par_BrM_LM_atac.merge.nomt.pt.mqsd.cutsite.sharp.daughtery_peaks.narrowPeak", format = "narrowPeak")

hits <- findOverlaps(consens, BrM)
overlaps <- pintersect(consens[queryHits(hits)], BrM[subjectHits(hits)])
percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
BrM.consens <- consens[queryHits(hits)[percentOverlap > 0]]
BrM.consens <- as.data.table(BrM.consens)[,.(seqnames, start, end)]
BrM.consens[, peakid := paste(seqnames, start, end, sep = "-")]
BrM.consens[, brm := 1]

hits <- findOverlaps(consens, LM)
overlaps <- pintersect(consens[queryHits(hits)], LM[subjectHits(hits)])
percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
LM.consens <- consens[queryHits(hits)[percentOverlap > 0]]
LM.consens <- as.data.table(LM.consens)[,.(seqnames, start, end)]
LM.consens[, peakid := paste(seqnames, start, end, sep = "-")]
LM.consens[, lm := 1]

hits <- findOverlaps(consens, Par)
overlaps <- pintersect(consens[queryHits(hits)], Par[subjectHits(hits)])
percentOverlap <- width(overlaps)/width(consens[queryHits(hits)])
Par.consens <- consens[queryHits(hits)[percentOverlap > 0]]
Par.consens <- as.data.table(Par.consens)[,.(seqnames, start, end)]
Par.consens[, peakid := paste(seqnames, start, end, sep = "-")]
Par.consens[, par := 1]

consens.dt <- as.data.table(consens)[,.(seqnames, start, end)]
consens.dt[, peakid := paste(seqnames, start, end, sep = "-")]

consens.dt <- merge(consens.dt[!(duplicated(peakid)),.(peakid)], BrM.consens[!(duplicated(peakid)),.(peakid, brm)], by = "peakid", all = TRUE)
consens.dt <- merge(consens.dt, LM.consens[!(duplicated(peakid)),.(peakid, lm)], by = "peakid", all = TRUE)
consens.dt <- merge(consens.dt, Par.consens[!(duplicated(peakid)),.(peakid, par)], by = "peakid", all = TRUE)

consens.dt <- cbind(consens.dt[,.(peakid)],
                    apply(consens.dt[,.(par,brm,lm)], c(1,2), function(x){if(is.na(x)){return(0)}else{return(x)}}))
fwrite(consens.dt, "atac.peaks.txt", sep = "\t")

#### Regular Version ####

h3k4me3 <- fread("H3K4me3.peaks.txt")
colnames(h3k4me3) <- c("peakid", "Par", "BrM2", "LM2")
upset(h3k4me3, sets = c("BrM2", "LM2", "Par"),
      order.by = "degree", keep.order = TRUE, point.size = 3,
      text.scale = c(2, 2, 3, 3, 2, 1.75), number.angles = 0, set_size.show = FALSE)

h3k27ac <- fread("H3K27ac.peaks.txt")
colnames(h3k27ac) <- c("peakid", "Par", "BrM2", "LM2")
upset(h3k27ac, sets = c("BrM2", "LM2", "Par"),
      order.by = "degree", keep.order = TRUE, point.size = 3,
      text.scale = c(2, 2, 3, 3, 2, 1.75), number.angles = 0, set_size.show = FALSE)

atac <- fread("atac.peaks.txt")
colnames(atac) <- c("peakid", "Par", "BrM2", "LM2")
upset(atac, sets = c("BrM2", "LM2", "Par"),
      order.by = "degree", keep.order = TRUE, point.size = 3,
      text.scale = c(2, 2, 3, 3, 2, 1.75), number.angles = 0, set_size.show = FALSE)
