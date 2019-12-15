library(Gviz)
library(data.table)
library(rtracklayer)
library(GenomicFeatures)

source("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/190720_load_files.R")

mainwd <- "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/"
inputfolder <- "input/"
outputfolder <- "output/IGV"
dir.create(file.path(mainwd, inputfolder), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainwd, outputfolder), recursive = TRUE, showWarnings = FALSE)
setwd(file.path(mainwd, outputfolder))

#library(biomaRt)
#hsapiensEnsembl <- makeTxDbFromBiomart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
#saveDb(hsapiensEnsembl, file = "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/rdata/hsapiensEnsembl.1toMTactive.txdb.sqlite")
hsapiensEnsembl <- loadDb(file = "/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/rdata/hsapiensEnsembl.1toMTactive.txdb.sqlite")
seqlevels(hsapiensEnsembl) <- seqlevels(hsapiensEnsembl)[1:25]

newscheme <- getScheme(name="default")
newscheme$GdObject$fontcolor.title <- "black"
newscheme$GdObject$background.title <- "white"
newscheme$GdObject$fontfamily.title <- "sans"
newscheme$GdObject$fontface.title <- 1
newscheme$GdObject$background.panel <- "transparent"
newscheme$GdObject$showAxis <- TRUE
newscheme$GdObject$col.axis <- "white"
newscheme$GdObject$cex.axis <- 0.2
newscheme$GdObject$rotation.title <- 90
newscheme$DataTrack$col.baseline <- "black"
newscheme$DataTrack$col <- "black"
newscheme$AnnotationTrack$showFeatureId <- TRUE
newscheme$AnnotationTrack$fontcolor.item <- "black"
newscheme$AnnotationTrack$fontsize <- 10
newscheme$AnnotationTrack$just.group <- "right"
newscheme$AnnotationTrack$just <- "right"
newscheme$AnnotationTrack$col <- "white"
newscheme$AnnotationTrack$fill <- "white"
newscheme$GeneRegionTrack$col <- "black"
newscheme$GeneRegionTrack$fill <- "black"

addScheme(newscheme, name = "custom")
options(Gviz.scheme = "custom")
options(ucscChromosomeNames=FALSE)

allgenes <- fread("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/resources/genomes/genes/output/ensembl_hgnc_map.txt")
load("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/output/inter.final.dt.all.RData")

symbol <- "IL1B"
chr <- "2"
start <- 113566813
end <- 113599447
inter.filter <- 0

#goi <- allgenes[hgnc_symbol == "SPARC", get("ensembl_gene_id")]
inter.final.all.dt.annote <- inter.final.all.dt
inter.final.all.dt.annote[is.na(annote1), annote1 := "unknown"]
inter.final.all.dt.annote[is.na(annote2), annote2 := "unknown"]
#inter.final.all.dt.annote <- inter.final.all.dt.annote[!(annote1 == "unknown" & annote2 == "unknown")]
inter.final.all.dt.annote[, gene1 := gsub("_.*", "", gene1)]
inter.final.all.dt.annote[, gene2 := gsub("_.*", "", gene2)]
inter.final.all.dt.annote <- inter.final.all.dt.annote[!(annote1 == "promoter" & annote2 == "promoter" & gene1 == gene2),]

transcript.of.interest <- "ENST00000407997"
IGVposition <- "chr22:39,457,286-39,485,270"
symbol <- "APOBEC3G"
igvplot <- function(IGVposition = NULL, symbol, type = "histogram",
                    inter.filter = 0, hichip.max = NULL,
                    connections = c("enhancer", "promoter"),
                    transcript.of.interest = NA,
                    atac.max = NA, h3k27ac.max = NA, h3k4me3.max = NA, rna.max = NA,
                    specific.links = NA, cell.lines = c("par", "lm", "brm")){
  
  # Parse IGVposition
  chr <- gsub(":.*|chr", "", IGVposition)
  start <- as.integer(gsub(",", "", gsub("^.*:|-.*", "", IGVposition)))
  end <- as.integer(gsub(",", "", gsub("^.*-", "", IGVposition)))
  
  #### Get gene by symbol ####
  ensembl.key <- allgenes[hgnc_symbol == symbol, get("ensembl_gene_id")][1]
  toi <- select(hsapiensEnsembl, keys = ensembl.key, columns = c("GENEID", "TXNAME"), keytype = "GENEID")
  if(is.na(transcript.of.interest)){
    trans.res <- transcripts(hsapiensEnsembl, filter=list(tx_name = toi$TXNAME))
    trans <- trans.res[which.max(width(trans.res))]$tx_name
  } else {
    trans <- transcript.of.interest
  }
  exon.test <- as.data.table(exons(hsapiensEnsembl, filter=list(tx_name = trans)))
  exon.test[,feature := "exon"]
  # cds.test <- as.data.table(cds(hsapiensEnsembl, filter=list(tx_name = longest.trans)))
  # cds.test[,feature := "protein_coding"]
  # final.text <- rbind(exon.test[,.(chromosome = seqnames, start, end, strand, feature)],
  #                     cds.test[,.(chromosome = seqnames, start, end, strand, feature)])
  final.text <- exon.test
  final.text[, gene := ensembl.key]
  final.text[, transcript := trans]
  final.text[, symbol := symbol]
  message(final.text)
  genetrack <- GeneRegionTrack(final.text, genome = "hg19", chromosome = chr, name = "Genes")#, transcriptAnnotation = "symbol")
  #plotTracks(genetrack)
  #### Get gene by symbol ####
  

  annote.width <- (end-start)/10
  sel = BigWigSelection(ranges = GRanges(seqnames = chr, ranges = IRanges(start = start, end = end)))
  atac.list <- list()
  h3k27ac.list <- list()
  h3k4me3.list <- list()
  rna.list <- list()
  list.maxes <- list()
  for(line in cell.lines){
    fname <- lines.dict[line.lower == line, get("line.filename")]
    atac.list[[line]] <- import.bw(paste0("/Users/wcai/Google_Drive/_Lab/Data/ATAC-seq/bigwig/MDAMB231_", fname, "_atac.merge.nomt.bs5.RPKM.bw"), selection = sel, as="GRanges")
    h3k27ac.list[[line]] <- import.bw(paste0("/Users/wcai/Google_Drive/_Lab/Data/ChIP-seq/bigwig/MDAMB231_", fname,"_H3K27ac_merge.mqsd.coverage.RPKM.bw"), selection = sel, as="GRanges")
    h3k4me3.list[[line]] <- import.bw(paste0("/Users/wcai/Google_Drive/_Lab/Data/ChIP-seq/bigwig/MDAMB231_", fname,"_H3K4me3_merge.mqsd.coverage.RPKM.bw"), selection = sel, as="GRanges")
    rna.list[[line]] <- import.bw(paste0("/Users/wcai/Google_Drive/_Lab/Data/RNA-seq/E190311_MDA-MB-231_RNA-seq/bw/", line,".rna.merge.bw"), selection = sel, as="GRanges")
    if(length(rna.list[[line]]) == 0){
      rna.list[[line]] <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
      
    }
    
    list.maxes$atac <- c(list.maxes$atac, max(atac.list[[line]]$score, na.rm = TRUE))
    list.maxes$h3k27ac <- c(list.maxes$h3k27ac, max(h3k27ac.list[[line]]$score, na.rm = TRUE))
    list.maxes$h3k4me3 <- c(list.maxes$h3k4me3, max(h3k4me3.list[[line]]$score, na.rm = TRUE))
    list.maxes$rna <- c(list.maxes$rna, max(rna.list[[line]]$score, na.rm = TRUE))
  }
  
  # Manual maxes
  for(signal.type in names(list.maxes)){
    type.max <- get(paste0(signal.type, ".max"))
    message(type.max)
    if(!is.na(type.max)) list.maxes[[signal.type]] <- rep(type.max, length(list.maxes[[signal.type]]))
  }
  
  atac.tracklist <- list()
  h3k27ac.tracklist <- list()
  h3k4me3.tracklist <- list()
  rna.tracklist <- list()
  line <- "par"
  title.color <- list()
  title.color[["par"]] <- "dimgray"
  title.color[["brm"]] <- "dodgerblue4"
  title.color[["lm"]] <- "darkorange1"
  for(line in cell.lines){
    atac.tracklist[[line]] <- DataTrack(atac.list[[line]], genome = "hg19", name = lines.dict[line.lower == line, get("line.fixed")], ylim = c(0,max(list.maxes$atac)), fill.mountain = c("white", "black"), fill.histogram = "black", col.histogram = NA, baseline = 0, background.title = title.color[[line]], col.title = title.color[[line]], col.axis = title.color[[line]])
    h3k27ac.tracklist[[line]] <- DataTrack(h3k27ac.list[[line]], genome = "hg19", name = lines.dict[line.lower == line, get("line.fixed")], ylim = c(0,max(list.maxes$h3k27ac)), col.baseline = "black", fill.mountain = c("white", "dodgerblue3"), fill.histogram = "dodgerblue3", col.histogram = NA, baseline = 0, background.title = title.color[[line]], col.title = title.color[[line]], col.axis = title.color[[line]])
    h3k4me3.tracklist[[line]] <- DataTrack(h3k4me3.list[[line]], genome = "hg19", name = lines.dict[line.lower == line, get("line.fixed")], ylim = c(0,max(list.maxes$h3k4me3)), col.baseline = "black", fill.mountain = c("white", "green4"), fill.histogram = "green4", col.histogram = NA, baseline = 0, background.title = title.color[[line]], col.title = title.color[[line]], col.axis = title.color[[line]])
    rna.tracklist[[line]] <- DataTrack(rna.list[[line]], genome = "hg19", name = lines.dict[line.lower == line, get("line.fixed")], ylim = c(0,max(list.maxes$rna)), col.baseline = "black", fill.mountain = c("white", "firebrick"), fill.histogram = "firebrick", col.histogram = NA, baseline = 0, background.title = title.color[[line]], col.title = title.color[[line]], col.axis = title.color[[line]])
  }
  
  if(!is.na(connections)){
    # HiChIP
    hichip <- import.bw("/Users/wcai/Google_Drive/_Lab/Data/Bioinformatics/190720_revamp/analysis/hichip/bigwig/merge.cut.nochr.combined.filterchrchr.slop.sort.bw", selection = sel, as="GRanges")
    if(is.null(hichip.max)) hichip.max <- max(hichip$score, na.rm = TRUE)
    hichip.track <- DataTrack(hichip, genome = "hg19", name = "HiChIP 1D", ylim = c(0, hichip.max), fill.mountain = c("white", "purple"), fill.histogram = "purple", col.histogram = NA, baseline = 0, background.title = "white", col.title = "black", col.axis = "white")
    
    # Sashimi
    goi.hichip <- rbind(inter.final.all.dt.annote[grep(ensembl.key, gene1),.(source = gene1, source.annote = annote1,
                                                                             target.annote = annote2, intercount.sum,
                                                                             source.peakid = atac1.peakid,
                                                                             target.peakid = atac2.peakid)],
                        inter.final.all.dt.annote[grep(ensembl.key, gene2),.(source = gene2, source.annote = annote2,
                                                                             target.annote = annote1, intercount.sum,
                                                                             source.peakid = atac2.peakid,
                                                                             target.peakid = atac1.peakid)])

    if(!is.na(specific.links)){
      message(specific.links)
      goi.hichip <- goi.hichip[target.peakid %in% specific.links & source.peakid %in% specific.links]
    }
    
    goi.hichip <- goi.hichip[!duplicated(goi.hichip[,.(source.peakid, target.peakid)]),]
    goi.hichip[, chr1 := as.character(sub("-.*", "", source.peakid))]
    goi.hichip[, chr2 := as.character(sub("-.*", "", source.peakid))]
    goi.hichip[, read1 := as.integer(sub("-.*", "", sub("^.+?-", "", source.peakid)))]
    goi.hichip[, read2 := as.integer(sub("-.*", "", sub("^.+?-", "", target.peakid)))]
    #goi.hichip <- goi.hichip[order(intercount.sum, decreasing = TRUE)]
    # Filter by parameters
    #connections <- "enhancer"
    goi.hichip <- goi.hichip[intercount.sum > inter.filter & target.annote %in% connections]#[atac1.peakid == "5-151066380-151066890" | atac2.peakid == "5-151066380-151066890" 
    goi.hichip[,newread1 := read1]
    goi.hichip[,newread2 := read2]
    goi.hichip[read1 > read2,newread1 := read2]
    goi.hichip[read1 > read2,newread2 := read1]
    sashimidata <- goi.hichip[,.(chr = chr1, start1 = newread1, end1 = newread1+1, start2 = newread2, end2 = newread2+1, weights = ceiling(intercount.sum/2))]
    # Filter for closer connections
    filterDist <- 1e5
    sashimidata <- sashimidata[start2 < filterDist + start1,]
    sashimidata[, cigar := paste0(end1-start1,"M",start2-end1,"N",end2-start2,"M")]
    sashimidata <- sashimidata[!duplicated(cigar),]
    sashimidata <- sashimidata[,.(chr, start1, end1, start2, end2, seqweights = seq(1:weights)), by = .(cigar)]
    track_ep_pairs <- AlignmentsTrack(range=with(sashimidata,GRanges(chr,IRanges(start1,end1),"*",cigar)),genome="hg19",isPaired=FALSE, name = "Linkage")
  }
  
  blank.annote <- AnnotationTrack(start = c(start+annote.width/2), width = c(annote.width), id = c(" "), chromosome = chr)
  atac.annote <- AnnotationTrack(start = c(start+annote.width/2), width = c(annote.width), id = c(paste0("0 - ", ceiling(max(list.maxes$atac, na.rm = TRUE)))), chromosome = chr)
  h3k27ac.annote <- AnnotationTrack(start = c(start+annote.width/2), width = c(annote.width), id = c(paste0("0 - ", ceiling(max(list.maxes$h3k27ac, na.rm = TRUE)))), chromosome = chr)
  h3k4me3.annote <- AnnotationTrack(start = c(start+annote.width/2), width = c(annote.width), id = c(paste0("0 - ", ceiling(max(list.maxes$h3k4me3, na.rm = TRUE)))), chromosome = chr)
  rna.annote <- AnnotationTrack(start = c(start+annote.width/2), width = c(annote.width), id = c(paste0("0 - ", ceiling(max(list.maxes$rna, na.rm = TRUE)))), chromosome = chr)
  
  axisTrack <- GenomeAxisTrack(labelPos = "below", scale = 0.5)

  num.cell.lines <- length(cell.lines)
  if(!is.na(connections)){
    all.tracks <- c(axisTrack, atac.annote, atac.tracklist, track_ep_pairs, #hichip.track,
                    blank.annote, h3k27ac.annote, h3k27ac.tracklist, blank.annote,
                    h3k4me3.annote, h3k4me3.tracklist, blank.annote, rna.annote, rna.tracklist,
                    blank.annote, genetrack)
    all.sizes <- c(4,4,rep(4,num.cell.lines),12,1,4,rep(4,num.cell.lines),1,4,rep(4,num.cell.lines),1,4,rep(4,num.cell.lines),1,2)
    all.types <- c(type, "sashimi")
  } else {
    all.tracks <- c(axisTrack, atac.annote, atac.tracklist, #hichip.track,
                    blank.annote, h3k27ac.annote, h3k27ac.tracklist, blank.annote,
                    h3k4me3.annote, h3k4me3.tracklist, blank.annote, rna.annote, rna.tracklist,
                    blank.annote, genetrack)
    all.sizes <- c(4,4,rep(4,num.cell.lines),1,4,rep(4,num.cell.lines),1,4,rep(4,num.cell.lines),1,4,rep(4,num.cell.lines),1,2)
    all.types <- c(type)
  }
  plotTracks(all.tracks, sizes = all.sizes, from = start, to = end, type = all.types)
  
}

#EMP2 chr16:10,372,237-11,139,646
symbol <- "EMP2"
IGVposition <- "chr16:10,661,351-10,681,691"
inter.filter <- 5
png(file = paste0("EMP2_3.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
igvplot(IGVposition = IGVposition, symbol = symbol, type = "histogram",
        inter.filter = inter.filter, hichip.max = 0.2, connections = "enhancer")
dev.off()


#### Manual plotting ####
#RIMS2
inter <- 0
png(file = paste0("RIMS2_2.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr8:105,204,661-105,270,562", "RIMS2", type = "histogram", inter.filter = inter,
        transcript.of.interest = "ENST00000507740", atac.max = 60, h3k27ac.max = 90)
dev.off()

#EMP2
inter <- 5
png(file = paste0("EMP2_2.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr16:10,660,438-10,680,778", "EMP2", type = "histogram", inter.filter = inter,
        transcript.of.interest = "ENST00000359543", atac.max = 60, h3k27ac.max = 90)
dev.off()

#SPARC chr5:151,063,234-151,068,324
png(file = paste0("SPARC.promoter.png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr5:151,063,234-151,068,324", "SPARC", type = "histogram", inter.filter = inter, connections = NA)
dev.off()

#PLCB1 chr20:8,111,277-8,114,823
png(file = paste0("PLCB1.promoter.png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr20:8,111,277-8,114,823", "PLCB1", type = "histogram", inter.filter = inter, connections = NA)
dev.off()

#Promoter for BrM2: CDH18 ENST00000382275, chr5:19,985,619-19,990,709
png(file = paste0("CDH18.promoter.png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr5:19,985,619-19,990,709", "CDH18", type = "histogram",
        transcript.of.interest = "ENST00000502796", inter.filter = inter, connections = NA)
dev.off()

#PTK7 chr6:43,039,781-43,116,827
# for(inter in 1:4){
#   png(file = paste0("PTK7.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
#   igvplot("chr6:43,039,781-43,116,827", "PTK7", type = "histogram", inter.filter = inter,
#           transcript.of.interest = "ENST00000230419")
#   dev.off()
# }

#SOCS2 chr12:93,940,948-94,001,215
# inter <- 2
# png(file = paste0("SOCS2.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
# igvplot("chr12:93,940,948-94,001,215", "SOCS2", type = "histogram", inter.filter = inter)
# dev.off()

#ANKFN1 chr17:54,212,083-54,285,200
# inter <- 0
# png(file = paste0("ANKFN1.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
# igvplot("chr17:54,212,083-54,285,200", "ANKFN1", type = "histogram", inter.filter = inter)
# dev.off()

#EVA1C chr21:33,768,952-33,855,635
# inter <- 0
# specific.links <- c("21-33784840-33785132",
#                     "21-33835482-33835724",
#                     "21-33835724-33836042",
#                     "21-33836117-33836708",
#                     "21-33836757-33836982")
# png(file = paste0("EVA1C.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
# igvplot("chr21:33,768,952-33,855,635", "EVA1C", type = "histogram", inter.filter = inter,
#         connections = "enhancer", specific.links = specific.links, h3k27ac.max = 120)
# dev.off()

#APOBEC3G chr22:39,457,286-39,485,270
inter <- 0
specific.links <- c("22-39472765-39473406", "22-39462574-39463086")
png(file = paste0("APOBEC3G.inter_gt", inter, ".png"), units = "in", width = 4, height = 7, res = 300)
igvplot("chr22:39,457,286-39,485,270", "APOBEC3G", type = "histogram", inter.filter = inter,
        connections = "enhancer", transcript.of.interest = "ENST00000407997",
        specific.links = specific.links, h3k27ac.max = 100)
dev.off()

#### Manual plotting ####
