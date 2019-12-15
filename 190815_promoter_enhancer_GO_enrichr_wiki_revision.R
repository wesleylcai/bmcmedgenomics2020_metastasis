# Last Updated: 
# Author: Wesley Cai
# Purpose:

library(data.table)
library(ggplot2)
library(httr)
library(jsonlite)
library(reshape2)
library(openxlsx)

# load gsea
# source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/gsea/gmt/load_gmt.R")
# load for revamp
source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

databases <- c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "KEGG_2016",
              "WikiPathways_2016", "Reactome_2016", "Panther_2016")
#### Functions ####
# Main function
gene.list <- prom.res[[line]][[direction]][["hgnc"]]
databases <- "KEGG_2016"
padj.cutoff <- 0.05
# KEGG_2016_Human Panther_2016 Reactome_2016 WikiPathways_2016
enrichr <- function(gene.list, databases = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "KEGG_2016",
                                             "WikiPathways_2016", "Reactome_2016", "Panther_2016"), padj.cutoff = 0.05){
  # Send query using gene list
  req.body <- list(list = paste(gene.list, collapse = "\n"))
  post.req <- httr::POST("http://amp.pharm.mssm.edu/Enrichr/addList", encode = "multipart", body = I(req.body))
  query.id <- jsonlite::fromJSON(content(post.req, type = "text", encoding = "UTF-8"))$userListId
  # Retreive results in supplied databases
  enrichr.res <- list()
  #db <- "GO_Molecular_Function_2018"
  for(db in databases){
    get.req <- httr::GET(paste0("http://amp.pharm.mssm.edu/Enrichr/enrich?userListId=", query.id, "&backgroundType=", db))
    res <- jsonlite::fromJSON(content(get.req, type = "text", encoding = "UTF-8"))
    enrichr.res[[db]] <- as.data.table(do.call(rbind, res[[db]]))
    tmp <- enrichr.res[[db]]
    tmp <- tmp[,.(index = unlist(V1), name = unlist(V2),
                  pval = unlist(V3), odds_ratio = unlist(V4),
                  combined_score = unlist(V5), targets = V6,
                  padj = unlist(V7))]
    tmp <- tmp[padj < padj.cutoff,]
    tmp <- tmp[order(padj),]
    enrichr.res[[db]] <- tmp
  }
  if(length(databases) == 1){
    return(enrichr.res[[databases]])
  } else {
    return(enrichr.res)
  }
}
#### Functions ####

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/promoter_enhancerGO_wiki_revision"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

expressed <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/output/expressed_genes/par_brm2_lm2.expressed.txt")
inter <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.lin.enh.hgnc.txt")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.all.RData")

wb <- createWorkbook("enrichr_results")
for(database in databases){
  #### Promoter ####
  prom <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/MDAMB231_Par_BrM_LM_H3K4me3_merge_peaks.tss.bed")
  prom <- prom[V10 == 0,.(peakid = V4, ensembl = V8)]
  prom[, ensembl := gsub("_.*", "", ensembl)]
  prom.rna <- merge(prom, res.list[["rna"]][,.(ensembl, hgnc, rna.brm_l2fc = brm_l2fc,
                                               rna.brm_padj = brm_padj, rna.lm_l2fc = lm_l2fc,
                                               rna.lm_padj = lm_padj)], by = "ensembl")
  prom.rna.hist <- merge(prom.rna, res.list[["h3k4me3.broad"]][,.(peakid, hist.brm_l2fc = brm_l2fc,
                                                                  hist.brm_padj = brm_padj, hist.lm_l2fc = lm_l2fc,
                                                                  hist.lm_padj = lm_padj)], by = "peakid")
  prom.res <- list()
  padj.cutoff <- 0.05
  for(line in c("brm", "lm")){
    line.l2fc <- paste0(line, "_l2fc")
    line.padj <- paste0(line, "_padj")
    for(direction in c("down", "up")){
      if(direction == "down"){
        prom.res[[line]][[direction]] <-prom.rna.hist[hgnc != "" & get(paste0("rna.", line.l2fc)) < 0 & get(paste0("hist.", line.l2fc)) < 0 & get(paste0("rna.", line.padj)) < padj.cutoff & get(paste0("hist.", line.padj)) < padj.cutoff, .(hgnc)]
      } else {
        prom.res[[line]][[direction]] <-prom.rna.hist[hgnc != "" & get(paste0("rna.", line.l2fc)) > 0 & get(paste0("hist.", line.l2fc)) > 0 & get(paste0("rna.", line.padj)) < padj.cutoff & get(paste0("hist.", line.padj)) < padj.cutoff, .(hgnc)]
      }
    }
  }
  prom.res.go.dt <- data.table()
  for(line in names(prom.res)){
    for(direction in names(prom.res[[line]])){
      #message(paste(line, direction))
      #print(prom.res[[line]][[direction]][["hgnc"]])
      tmp <- enrichr(prom.res[[line]][[direction]][["hgnc"]], databases = database, padj.cutoff = 1.1)
      tmp$line <- line
      tmp$direction <- direction
      prom.res.go.dt <- rbind(prom.res.go.dt, tmp[1:20])
    }
  }
  prom.res.go.dt[, targets.paste := do.call(rbind, lapply(targets, function(x) paste(sort(x), collapse = ";")))]
  prom.res.go.dt[, targets.num := do.call(rbind, lapply(targets, length))]
  prom.res.go.dt <- prom.res.go.dt[,.SD[which.min(padj)], by = .(targets.paste, line, direction)]
  prom.res.go.dt <- prom.res.go.dt[targets.num > 1 & odds_ratio > 10]
  prom.res.go.dt[, type := "promoter"]
  #### Promoter ####
  
  #### Enhancer ####
  inter.annote <- inter.final.all.dt
  inter.annote[is.na(annote1), annote1 := "unknown"]
  inter.annote[is.na(annote2), annote2 := "unknown"]
  inter.annote <- inter.annote[!(annote1 == "unknown" & annote2 == "unknown")]
  inter.annote[, gene1 := gsub("_.*", "", gene1)]
  inter.annote[, gene2 := gsub("_.*", "", gene2)]
  inter.annote <- inter.annote[!(annote1 == "promoter" & annote2 == "promoter" & gene1 == gene2),]
  inter.annote.lin <- rbind(inter.annote[,.(peakid = atac1.peakid, target = gene2, annote = annote1, intercount.sum)],
                            inter.annote[,.(peakid = atac2.peakid, target = gene1, annote = annote2, intercount.sum)])
  inter.annote.lin <- inter.annote.lin[annote != "unknown"]
  inter.annote.lin <- inter.annote.lin[!duplicated(inter.annote.lin)]
  inter.annote.lin <- inter.annote.lin[!is.na(target)]
  atac.h3k27ac <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/peaks/atac.h3k27ac.broad.bed")
  
  ### Enhancer Only ###
  atac.h3k27ac.enh <- atac.h3k27ac[V8 %in% chip.annotations[["h3k27ac.broad.enh.active"]] & V9 == 0,]
  inter.annote.lin.h3k27ac <- merge(inter.annote.lin, atac.h3k27ac.enh[,.(peakid = V4, h3k27ac.peakid = V8)], by = "peakid", all.x = TRUE)
  inter.annote.lin.h3k27ac <- inter.annote.lin.h3k27ac[!is.na(h3k27ac.peakid)]
  rna.tmp <- res.list[["rna"]]
  colnames(rna.tmp)[grep("_", colnames(rna.tmp))] <- paste0("rna.", grep("_", colnames(rna.tmp), value = TRUE))
  inter.annote.lin.h3k27ac.res <- merge(inter.annote.lin.h3k27ac, rna.tmp, by.x = "target", by.y = "ensembl", all.x = TRUE)
  h3k27ac.tmp <- res.list[["h3k27ac.broad"]]
  colnames(h3k27ac.tmp)[grep("_", colnames(h3k27ac.tmp))] <- paste0("h3k27ac.", grep("_", colnames(h3k27ac.tmp), value = TRUE))
  inter.annote.lin.h3k27ac.res <- merge(inter.annote.lin.h3k27ac.res, h3k27ac.tmp, by.x = "h3k27ac.peakid", by.y = "peakid", all.x = TRUE)
  inter.annote.lin.h3k27ac.res <- inter.annote.lin.h3k27ac.res[(rna.brm_padj < 0.05 & h3k27ac.brm_padj < 0.05) |
                                                                 (rna.lm_padj < 0.05 & h3k27ac.lm_padj < 0.05)]
  inter.annote.lin.h3k27ac.res <- inter.annote.lin.h3k27ac.res[,.SD[which.max(intercount.sum)], by = .(h3k27ac.peakid)]
  
  enh.any.res <- list()
  for(line in c("brm", "lm")){
    line.l2fc <- paste0(line, "_l2fc")
    line.padj <- paste0(line, "_padj")
    for(direction in c("down", "up")){
      if(direction == "down"){
        enh.any.res[[line]][[direction]] <- unique(inter.annote.lin.h3k27ac.res[hgnc != "" & get(paste0("rna.", line.l2fc)) < 0 & get(paste0("h3k27ac.", line.l2fc)) < 0 &
                                                                                  get(paste0("rna.", line.padj)) < 0.05 & get(paste0("h3k27ac.", line.padj)) < 0.05, get("hgnc")])
      } else {
        enh.any.res[[line]][[direction]] <- unique(inter.annote.lin.h3k27ac.res[hgnc != "" & get(paste0("rna.", line.l2fc)) > 0 & get(paste0("h3k27ac.", line.l2fc)) > 0 &
                                                                                  get(paste0("rna.", line.padj)) < 0.05 & get(paste0("h3k27ac.", line.padj)) < 0.05, get("hgnc")])
      }
    }
  }
  enh.any.res.go.dt <- data.table()
  for(line in names(enh.any.res)){
    for(direction in names(enh.any.res[[line]])){
      #message(paste(line, direction))
      #print(prom.res[[line]][[direction]][["hgnc"]])
      tmp <- enrichr(enh.any.res[[line]][[direction]], databases = database, padj.cutoff = 1.1)
      tmp$line <- line
      tmp$direction <- direction
      enh.any.res.go.dt <- rbind(enh.any.res.go.dt, tmp[1:20])
    }
  }
  enh.any.res.go.dt[, targets.paste := do.call(rbind, lapply(targets, function(x) paste(sort(x), collapse = ";")))]
  enh.any.res.go.dt[, targets.num := do.call(rbind, lapply(targets, length))]
  enh.any.res.go.dt <- enh.any.res.go.dt[,.SD[which.min(padj)], by = .(targets.paste, line, direction)]
  enh.any.res.go.dt <- enh.any.res.go.dt[targets.num > 1 & odds_ratio > 10]
  enh.any.res.go.dt[, type := "enhancer"]
  ### Enhancer Only ###
  #### Enhancer ####
  
  #### Combine all ####
  go.combined <- rbind(prom.res.go.dt, enh.any.res.go.dt)#, rna.res.go.dt)
  # go.combined[, desc := gsub("epidermal growth factor receptor", "EGFR", gsub(" \\(.*| via.*", "", name))]
  # go.combined[, desc := gsub("cyclin-dependent protein serine/threonine kinase", "CDK", desc)]
  # go.combined[, desc := gsub("transmembrane ", "", desc)]
  go.combined[direction == "down", neglogp := log(pval)]
  go.combined[direction == "up", neglogp := -log(pval)]
  go.combined[, go.id := gsub("\\)", "", gsub("^.*\\(", "", go.combined$name))]
  go.combined <- go.combined[order(pval),]
  go.combined[, go.id := factor(go.id, unique(go.id))]
  go.combined[, line.type := paste0(line, "_", type)]
  go.combined <- go.combined[targets.num > 2,]
  go.combined.subset <- go.combined[,.(line, direction, type, name, targets.paste, targets.num, pval, odds_ratio, combined_score, padj)]
  #fwrite(go.combined, "GO.combined.txt", sep = "\t")
  
  addWorksheet(wb, sheetName = database, zoom = 150)
  writeData(wb, database, go.combined.subset)
}
saveWorkbook(wb, paste0("enrichr_results.xlsx"), overwrite = TRUE)

# Common enhancer-reg in GO
intersect(unlist(go.combined[line.type == "lm_enhancer",get("targets")]),
          unlist(go.combined[line.type == "brm_enhancer",get("targets")]))
# Common enhancer-reg in general
intersect(enh.any.res[["lm"]][["up"]],enh.any.res[["brm"]][["up"]])

## Cluster the GOs
go.combined.cluster <- list()
element <- 1
total.targets <- unique(unlist(go.combined$targets))
for(element in 1:length(go.combined$targets)){
  tmp.bool <- as.integer(total.targets %in% go.combined$targets[[element]])
  names(tmp.bool) <- total.targets
  go.combined.cluster[[as.character(go.combined$go.id[[element]])]] <- tmp.bool
}
go.combined.cluster <- as.matrix(do.call(rbind, go.combined.cluster))
d <- dist(go.combined.cluster, method = "euclidean")
hc1 <- hclust(d, method = "complete")
plot(hc1)
## Cluster the GOs

go.combined$go.id <- factor(go.combined$go.id, hc1$labels[hc1$order])
go.combined <- go.combined[order(go.id)]
go.combined$desc <- factor(go.combined$desc, unique(go.combined$desc))
p <- ggplot(go.combined, aes(type, y = desc)) +
  geom_point(aes(size = targets.num, color = neglogp)) +
  #scale_color_brewer(palette = "RdBu") +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
                         high = "#B2182B", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p <- p + facet_grid(.~line)
p
ggsave("heatmap2.png", p, width = 6, height = 6)

# Some of the descriptions are pretty similar, just "positive" or "negative"
# So gotta remove those
desc.gsub <- gsub("positive |negative ", "", go.combined$desc)
desc.gsub <- desc.gsub[duplicated(desc.gsub)]
go.combined.filter <- go.combined[!(desc %in% desc.gsub)]
go.combined.filter[line == "lm", line.fixed := "LM2"]
go.combined.filter[line == "brm", line.fixed := "BrM2"]
go.combined.filter$line.fixed <- factor(go.combined.filter$line.fixed, c("LM2", "BrM2"))
go.combined.filter[type == "enhancer", type.fixed := "Enhancer"]
go.combined.filter[type == "promoter", type.fixed := "Promoter"]
go.combined.filter$type.fixed <- factor(go.combined.filter$type.fixed, c("Promoter", "Enhancer"))

p <- ggplot(go.combined.filter, aes(type.fixed, y = desc)) +
  geom_point(aes(size = targets.num, color = neglogp)) +
  #scale_color_brewer(palette = "RdBu") +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7",
                         high = "#B2182B", midpoint = 0, space = "Lab",
                         na.value = "grey50", guide = "colourbar", aesthetics = "colour",
                         name = "-log(P)") +
  scale_size_continuous(name = "Targets") +
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"))
p <- p + facet_grid(.~line.fixed)
p
ggsave("heatmap3_updated.png", p, width = 6, height = 6)


library(ComplexHeatmap)
pdf("yo.pdf", width = 10, height = 10)
Heatmap(go.combined.cluster)
dev.off()
library(RColorBrewer)
brewer.pal(11, "RdBu")
#### Combine all ####