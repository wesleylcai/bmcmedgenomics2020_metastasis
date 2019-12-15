# Last Updated: 
# Author: Wesley Cai
# Purpose: I can use this for LM2 because it is technically a subpopulation of 231

library(data.table)
library(ggplot2)
library(enrichR)

#### Functions ####
writeWig <- function(goi, interactions){
  goi.chrom <- strsplit(chip.annotations[["atac.tss.promoters"]][grep(paste0("_", goi, "_"), gene),get("peakid")], "-")[[1]][1]
  goi.chrom <- paste0("variableStep chrom=chr", goi.chrom)
  interactions.goi <- interactions[unique(c(grep(goi, interactions$gene1), grep(goi, interactions$gene2))),]
  exclude <- intersect(grep(goi, interactions.goi$gene1), grep(goi, interactions.goi$gene2))
  if(length(exclude) > 0){
    interactions.goi <- interactions.goi[-exclude,]
  }
  #rep1.goi[grep("goi", gene1),][,.(atac2.peakid, paste0(intercount.sum, "_inter_", gsub("^.+?_", "", gene1)))]
  interactions.goi.lin <- rbind(interactions.goi[grep(goi, gene1),][,.(peakid = atac2.peakid, intercount.sum)], 
                                interactions.goi[grep(goi, gene2),][,.(peakid = atac1.peakid, intercount.sum)])
  interactions.goi.lin <- cbind(interactions.goi.lin, do.call(rbind, strsplit(interactions.goi.lin$peakid, "-")))
  interactions.goi.lin <- interactions.goi.lin[intercount.sum > 1,]
  if(nrow(interactions.goi.lin) > 0){
    res <- data.table()
    rown <- 1
    for(rown in 1:nrow(interactions.goi.lin)){
      res <- rbind(res, data.table(start = seq(interactions.goi.lin[rown,get("V2")],
                                               interactions.goi.lin[rown,get("V3")]),
                                   value = interactions.goi.lin[rown,get("intercount.sum")]))
    }
    res <- res[order(start, decreasing = FALSE),]
    fwrite(res, paste0(goi, ".intercount.tmp.wig"), sep = "\t", col.names = FALSE)
    system(paste0("echo ", goi.chrom, " | cat - ", paste0(goi, ".intercount.tmp.wig"), " > ", paste0(goi, ".intercount.wig"), " && rm ", paste0(goi, ".intercount.tmp.wig")))
  } else {
    message(paste0("No contact for ", goi))
  }
}
#### Functions ####



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
# hc.reps <- list()
# hc.reps[["rep1"]] <- fread("../associations/rep1.unique.intercount.txt")
# hc.reps[["rep2"]] <- fread("../associations/rep2.unique.intercount.txt")
# 
# hc.genes <- list()
# tmp.atac.tss <- chip.annotations[["atac.tss.promoters"]][,.(peakid, gene)]
# tmp.atac.tss[, gene := gsub("_.*", "", gene)]
# for(rep in names(hc.reps)){
#   # Determine percent of ATAC with a connection (looks like about 92% of peaks have a connection)
#   #message(length(c(unique(hc.reps[[rep]]$atac1.peakid, hc.reps[[rep]]$atac2.peakid)))/length(unique(chip.annotations[["atac.all"]])))
#   
#   # Need to annotate first to prune list of interactions
#   hc.reps[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote1 := "enhancer"]
#   hc.reps[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote2 := "enhancer"]
#   hc.reps[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote1 := "promoter"]
#   hc.reps[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote2 := "promoter"]
#   #table(rep1[,.(annote1, annote2)])
#   
#   hc.genes[[rep]] <- merge(hc.reps[[rep]], chip.annotations[["atac.tss.promoters"]][,.(peakid, gene1 = gene)], by.x = "atac1.peakid", by.y = "peakid", all.x = TRUE)
#   hc.genes[[rep]] <- merge(hc.genes[[rep]], chip.annotations[["atac.tss.promoters"]][,.(peakid, gene2 = gene)], by.x = "atac2.peakid", by.y = "peakid", all.x = TRUE)
#   
#   # Remove rows with NA in both gene columns
#   hc.genes[[rep]] <- hc.genes[[rep]][!(is.na(gene1) & is.na(gene2))]
# 
#   # Remove rows with no annotation in either annote columns
#   hc.genes[[rep]] <- hc.genes[[rep]][!(is.na(annote1) | is.na(annote2))]
#   
#   # Count up all rows where the interactions are between two of the same peaks
#   tmp <- hc.genes[[rep]][,.(id = paste0(atac1.peakid, "_", atac2.peakid),
#                            atac1.peakid, atac2.peakid, intercount.sum)]
#   tmp2 <- rbind(tmp[,.(id, peak1 = atac1.peakid, peak2 = atac2.peakid, intercount.sum)],
#                 tmp[,.(id, peak1 = atac2.peakid, peak2 = atac1.peakid, intercount.sum)])
#   tmp2 <- tmp2[,.(id, sum.final = sum(intercount.sum)), by = .(peak1, peak2)]
#   hc.genes[[rep]] <- tmp2[!duplicated(id) & !duplicated(tmp2[,.(peak1, peak2)]),
#                          .(atac1.peakid = peak1, atac2.peakid = peak2, intercount.sum = sum.final)]
#   
#   # Redo annotations
#   hc.genes[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote1 := "enhancer"]
#   hc.genes[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k27ac.broad.enh.active"]], annote2 := "enhancer"]
#   hc.genes[[rep]][atac1.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote1 := "promoter"]
#   hc.genes[[rep]][atac2.peakid %in% chip.annotations[["atac.h3k4me3.broad.prom"]], annote2 := "promoter"]
#   hc.genes[[rep]] <- merge(hc.genes[[rep]], tmp.atac.tss[,.(peakid, gene1 = gene)], by.x = "atac1.peakid", by.y = "peakid", all.x = TRUE)
#   hc.genes[[rep]] <- merge(hc.genes[[rep]], tmp.atac.tss[,.(peakid, gene2 = gene)], by.x = "atac2.peakid", by.y = "peakid", all.x = TRUE)
#   hc.genes[[rep]] <- hc.genes[[rep]][!(is.na(gene1) & is.na(gene2))]
#   hc.genes[[rep]] <- hc.genes[[rep]][!(is.na(annote1) | is.na(annote2))]
#   hc.genes[[rep]][gene1 == gene2, equal := TRUE]
#   hc.genes[[rep]] <- hc.genes[[rep]][is.na(equal),.(atac1.peakid, atac2.peakid, intercount.sum, annote1, annote2, gene1, gene2)]
#   #test <- hc.genes[[rep]]
# }
# rm(hc.reps)
# save(hc.genes, file = "hc.genes.RData")
# load(file = "hc.genes.RData")
#### Get region-region interactions that are associated with H3K27ac or H3K4me3 ####

#### Get region-region interactions that exist in both replicates ####
# inter.rep1 <- paste0(hc.genes[["rep1"]]$atac1.peakid, "_", hc.genes[["rep1"]]$atac2.peakid)
# inter.rep2 <- paste0(hc.genes[["rep2"]]$atac1.peakid, "_", hc.genes[["rep2"]]$atac2.peakid)
# inter.rep2.swap <- paste0(hc.genes[["rep2"]]$atac2.peakid, "_", hc.genes[["rep2"]]$atac1.peakid)
# inter.rep1.rep2 <- intersect(inter.rep1, inter.rep2)
# inter.rep1.rep2.swap <- intersect(inter.rep1, inter.rep2.swap)
# inter.final <- unique(inter.rep1.rep2, inter.rep1.rep2.swap)
# 
# tmp <- hc.genes[["rep1"]][,.(atac1.peakid, atac2.peakid, annote1, annote2, gene1, gene2)]
# tmp[, id := paste0(atac1.peakid, "_", atac2.peakid)]
# inter.final.dt <- tmp[id %in% inter.final,.(atac1.peakid, atac2.peakid, annote1, annote2, gene1, gene2)]
# save(inter.final.dt, file = "inter.final.dt.RData")
load("inter.final.dt.RData")
#### Get region-region interactions that exist in both replicates ####

length(which(duplicated(tmp[,.(atac1.peakid, atac2.peakid)])))
length(which(duplicated(tmp2[,.(peak1, peak2)])))
enh.genes <- hc.genes[["rep1"]]
enh.genes[, gene1 := gsub("_.*", "", gene1)]
enh.genes[, gene2 := gsub("_.*", "", gene2)]
enh.genes[gene1 == gene2, equal := TRUE]
enh.genes <- enh.genes[is.na(equal)]
#enh.genes <- enh.genes[annote1 == "enhancer" | annote2 == "enhancer"]
enh.genes.lin <- rbind(enh.genes[,.(peakid = atac2.peakid, intercount.sum, annote = annote2, gene = gene1)],
                   enh.genes[,.(peakid = atac1.peakid, intercount.sum, annote = annote1, gene = gene2)])
#enh.genes.lin <- enh.genes.lin[annote == "enhancer"]
enh.genes.lin <- enh.genes.lin[,.(interaction.sum = sum(intercount.sum)), by = .(peakid, gene)]
enh.genes.lin.merged <- merge(enh.genes.lin, res.list[["atac"]][,.(peakid, lm_l2fc, lm_padj)])
enh.genes.lin.merged.rna <- merge(enh.genes.lin.merged, rna[,.(gene = ensembl, hgnc_symbol, lm2.l2fc, lm2.padj)], by = "gene")
hist <- chip.annotations[["atac.h3k27ac.broad.overlap"]][,.(peakid = atac.peakid, histone = h3k27ac.peakid)]
hist <- merge(hist, res.list[["h3k27ac.broad"]][,.(histone = peakid, hist.lm_l2fc = lm_l2fc, hist.lm_padj = lm_padj)], all.x = TRUE)
enh.genes.lin.merged.rna.histone <- merge(enh.genes.lin.merged.rna, hist, by.x = "peakid", by.y = "peakid", all.x = TRUE)

enh.genes.lin.merged.rna.histone.sig <- enh.genes.lin.merged.rna.histone[hist.lm_padj < 0.05]
enh.genes.lin.merged.rna.histone.sig <- enh.genes.lin.merged.rna.histone.sig[interaction.sum > 1]
tmp <- enh.genes.lin.merged.rna.histone.sig[lm2.padj < 0.05,]
cor.test(tmp$hist.lm_l2fc, tmp$lm2.l2fc)
ggplot(tmp, aes(hist.lm_l2fc, lm2.l2fc)) +
  geom_point()

enh.genes.lin.merged.sig.rna.motif <- merge(enh.genes.lin.merged.sig.rna, motif.occur.list[["lm"]], by.x = "peakid", by.y = "peakid", all.x = TRUE)
motifs.sel <- c("motif1", "motif2", "motif3", "motif4", "motif5")
selection <- c(paste0("up.", motifs.sel), paste0("down.", motifs.sel))
enh.subset <- enh.genes.lin.merged.sig.rna.motif[,.SD,.SDcols = c(1:8, which(colnames(enh.genes.lin.merged.sig.rna.motif) %in% selection))]

matrix.genes <- as.matrix(enh.subset[lm2.padj < 1e-20,.SD,.SDcols = grep("motif", colnames(enh.subset))])
row.names(matrix.genes) <- enh.subset[lm2.padj < 1e-20,get("hgnc_symbol")]
matrix.genes <- apply(matrix.genes, c(1,2), function(x){if(is.na(x)){0}else{1}})
pdf("heatmap.inter3.pdf", height = 20, width = 10)
print(ComplexHeatmap::Heatmap(matrix.genes))
dev.off()



enh.gene.hgnc <- enh.genes.lin.merged.sig[!duplicated(gene)]
enh.gene.hgnc <- merge(enh.gene.hgnc[,.(gene = gene)], rna[,.(gene = ensembl, hgnc_symbol, lm2.l2fc, lm2.padj)], by = "gene")
fwrite(data.table(hgnc = enh.gene.hgnc[hgnc_symbol != "" & lm2.l2fc > 0 & lm2.padj < 0.05,get("hgnc_symbol")]), "enh.up.genes.txt")
#enh.genes.lin.merged <- merge(enh.genes.lin, res.list[["atac"]][,.(peakid, lm_l2fc, lm_padj)])

rm(tmp)
tmp <- hc.genes[[1]]
tmp[, gene1 := gsub("_.*", "", gene1)]
tmp[, gene2 := gsub("_.*", "", gene2)]
tmp.sample <- tmp[sample(1:nrow(tmp), size = 1e5),]

gene.id <- rna[hgnc_symbol == "SPARC", get("ensembl")]
tmp.sparc <- tmp[(gene1 == gene.id | gene2 == gene.id),]
tmp.sparc[gene1 == gene2, equal := TRUE]
tmp.sparc <- tmp.sparc[is.na(equal)]
#tmp.sparc <- tmp.sparc[gene1 != gene2,]
tmp.sparc.lin <- rbind(tmp.sparc[gene1 == gene.id,.(peakid = atac2.peakid, intercount.sum, annote = annote2)],
                   tmp.sparc[gene2 == gene.id,.(peakid = atac1.peakid, intercount.sum, annote = annote1)])
tmp.sparc.lin <- merge(tmp.sparc.lin, res.list[["atac"]][,.(peakid, lm_l2fc, lm_padj)])
tmp.sparc.lin <- tmp.sparc.lin[lm_padj < 0.05,]
tmp.sparc.lin <- tmp.sparc.lin[,.(intercount.sum.sum = sum(intercount.sum), lm_l2fc.unique = unique(lm_l2fc), lm_padj.unique = unique(lm_padj)), by = .(peakid)]
# Write out peaks to get TFs
tmp.sparc.lin.inter.gt1 <- tmp.sparc.lin[intercount.sum.sum > 1,]
tmp.sparc.lin.peaks <- as.data.table(cbind(do.call(rbind, strsplit(tmp.sparc.lin.inter.gt1$peakid, "-")), paste0(tmp.sparc.lin.inter.gt1$intercount.sum.sum,
                                                                                                                 "_", signif(tmp.sparc.lin.inter.gt1$lm_l2fc.unique, 3),
                                                                                                                 "_", signif(tmp.sparc.lin.inter.gt1$lm_padj.unique, 3))))
tmp.sparc.lin.peaks[,V1 := paste0("chr", V1)]
fwrite(tmp.sparc.lin.peaks, "sparc.inter.bed", sep = "\t", col.names = FALSE)



tmp.sparc.lin$intercount.sum <- factor(tmp.sparc.lin$intercount.sum, c(sort(unique(tmp.sparc.lin$intercount.sum))))
mp <- ggplot(tmp.sparc.lin, aes(intercount.sum, lm_l2fc)) +
  geom_boxplot()
ggsave("SPARC_interactions.png", mp)

# Looks like peaks with 6,7,10 interactions are decreased in LM2 (these may be insulators!!)

tmp.sample.1 <- tmp[gene1 == gene.id,.(peakid = atac2.peakid, intercount.sum, annote = annote2)]
tmp.sample.2 <- tmp[gene2 == gene.id,.(peakid = atac1.peakid, intercount.sum, annote = annote1)]
tmp.sample.total <- rbind(tmp.sample.1, tmp.sample.2)
tmp.sample.total <- merge(tmp.sample.total, res.list[["atac"]][,.(peakid, lm_l2fc, lm_padj)], by = "peakid", all.x = TRUE)
tmp.sample.total <- tmp.sample.total[lm_padj < 0.05,]
tmp.sample.total$intercount.sum <- factor(tmp.sample.total$intercount.sum, c(sort(unique(tmp.sample.total$intercount.sum))))
ggplot(tmp.sample.total, aes(intercount.sum, lm_l2fc)) +
  geom_boxplot()
tmp.sample <- tmp.sample[annote1 == "enhancer" | annote2 == "enhancer",]
tmp.sample <- merge(tmp.sample,
                    res.list[["atac"]][,.(atac1.peakid = peakid, atac1.lm_l2fc = lm_l2fc, atac2.lm_padj = lm_padj)],
                    by = "atac1.peakid", all.x = TRUE)
tmp.sample <- merge(tmp.sample,
                    res.list[["atac"]][,.(atac2.peakid = peakid, atac2.lm_l2fc = lm_l2fc, atac2.lm_padj = lm_padj)],
                    by = "atac2.peakid", all.x = TRUE)
tmp.sample <- merge(tmp.sample,
                    rna[,.(gene1 = ensembl, gene1.lm2.l2fc = lm2.l2fc, gene1.lm2.padj = lm2.padj)],
                    by = "gene1", all.x = TRUE)
tmp.sample <- merge(tmp.sample,
                    rna[,.(gene2 = ensembl, gene2.lm2.l2fc = lm2.l2fc, gene2.lm2.padj = lm2.padj)],
                    by = "gene2", all.x = TRUE)



# Test results using GOI
goi <- "SPARC"
goi <- "TRIM32"
goi <- "MMP1"
goi <- "GATA3"
goi <- "TFAP2C"
goi <- "IGSF1"
gois <- c("PLCB1", "MMP3", "MUC15", "CD70", "ST6GALNAC3", "CSF3", "FOXA1", "FGF13", "DCLK1")
goi <- "DCLK1"
test <- chip.annotations[["atac.tss.promoters"]]
writeWig("TNC")
