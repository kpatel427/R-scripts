# to visualize copy number calls from CNV caller to see if there is a 17q amplification observed
#setwd("/Volumes/target_nbl_ngs/KP/KAT2A/")

#BiocManager::install("CopyNumberPlots")
library(CopyNumberPlots)

# s1.calls.file <- read.delim('cellLines/COGN426_NoIndex_L004.sorted.reheader.call.cns', header = T)
# s1.calls <- loadCopyNumberCalls(s1.calls.file)
# 
# s2.calls.file <- read.delim('cnvkit/cellLines/NB69_NoIndex_L008.sorted.reheader.call.cns', header = T)
# s2.calls <- loadCopyNumberCallsCNVkit('cnvkit/cellLines/NB69_NoIndex_L008.sorted.reheader.call.cns')

# kp <- plotKaryotype(chromosomes="chr17")
# plotCopyNumberCalls(kp, s2.calls)
# cn.cols <- getCopyNumberColors()
# legend("top", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols))
# By default we’ll see losses in green, 2n regions in gray and gains in yellow-orange-red.




files <- list.files(path = 'cnvkit/cellLines', pattern = '*.call.cns')
cn.calls <- list()

for(i in 1:length(files)){
  name <- files[i]
  name <- gsub('(.*?)_.*','\\1',name)
  assign(name,loadCopyNumberCallsCNVkit(paste0('cnvkit/cellLines/',files[i])))
  cn.calls[[as.character(name)]] <- get(name)
}


# for multiple Samples
cn.calls1 <- cn.calls[1:16]
pdf('CNCopyNumber17q_1.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosomes="chr17")
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")       
plotCopyNumberCalls(kp, cn.calls1, ymin = -1, ymax=10, r0=0, r=1, loh.height = 10, label.cex = 0.8)
dev.off()

cn.calls2 <- cn.calls[17:32]
pdf('CNCopyNumber17q_2.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosomes="chr17")
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")       
plotCopyNumberCalls(kp, cn.calls2, ymin = -1, ymax=10, r0=0, r=1, loh.height = 10, label.cex = 0.8)
dev.off()



#plotCopyNumberSummary(kp, cn.calls, r1=0.25)
#cn.cols <- getCopyNumberColors()
#legend("top", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols))



# script to generate CopyNumber plots for MDX project - for MYCN & c-MYC MDXs across phases
# setwd("/Volumes/target_nbl_ngs/KP/MDX_Minu/cnvkit/withoutControl")
library(CopyNumberPlots)
library(gridExtra)


# MYCN Samples ------------------
mycn.files <- list.files(path='.', pattern='MYCN.*.BQSR.call.cns')
cn.calls <- list()

for(i in 1:length(mycn.files)){
  name <- mycn.files[i]
  name <- gsub('.BQSR.call.cns','',name)
  assign(name,loadCopyNumberCallsCNVkit(mycn.files[i]))
  cn.calls[[as.character(name)]] <- get(name)
}


# for multiple Samples
# to order samples
cn.calls1 <- list("MYCN-0685-P0"=MYCN0685P0, "MYCN-3927-P1"=MYCN3927P1, "MYCN-6792-P1"=MYCN6792P1, "MYCN-289-P2"=MYCN289P2, "MYCN-2900-P2"=MYCN2900P2)

# running a loop to create plot for every chromsome
chrs <- c(seq(from = 1, to = 19, by = 1),'X','Y')
for(i in chrs){
  pdf(paste0('copyNumberPlots/MYCN_CNCopyNumber_',i,'.pdf'), width = 12, height = 9, onefile = F)
  plot.new()
  kp <- plotKaryotype(chromosome = paste0('chr',i), genome = 'mm10')
  kpAddCytobandLabels(kp)
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
  at <- autotrack(current.track = 1, total.tracks = 5)
  cn.cols <- getCopyNumberColors(colors = "red_blue")
  plotCopyNumberCalls(kp, cn.calls1, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
  dev.off()
}


# plotting across all chromosome - dividing into two plots
# plot1 = chr1-chr11
pdf('copyNumberPlots/MYCN_CNCopyNumber1to11_.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[1:11], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
dev.off()

# plot1 = chr12-chrY
pdf('copyNumberPlots/MYCN_CNCopyNumber12toY_.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
dev.off()


# plot.new()
# kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
# kpAddCytobandLabels(kp)
# cn.cols <- getCopyNumberColors(colors = "red_blue")
# plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
# legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
#        horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)



# CMYC Samples ------------------
cmyc.files <- list.files(path='.', pattern='c-Myc.*.BQSR.call.cns')
cn.calls.cmyc <- list()

for(i in 1:length(cmyc.files)){
  name <- cmyc.files[i]
  name <- gsub('.BQSR.call.cns','',name)
  assign(name,loadCopyNumberCallsCNVkit(cmyc.files[i]))
  cn.calls.cmyc[[as.character(name)]] <- get(name)
}


# for multiple Samples
# to order samples
cn.calls.cmyc <- list("c-Myc-W0172-P0"=`c-MycW0172P0`, "c-Myc-1340-P1"=`c-Myc1340P1`, "c-Myc-W0700-P1"=`c-MycW0700P1`, "c-Myc-7558-P2"=`c-Myc7558P2`)

# running a loop to create plot for every chromsome
chrs <- c(seq(from = 1, to = 19, by = 1),'X','Y')
for(i in chrs){
  pdf(paste0('copyNumberPlots/cMYC_CNCopyNumber_',i,'.pdf'), width = 15, height = 9, onefile = F)
  plot.new()
  kp <- plotKaryotype(chromosome = paste0('chr',i), genome = 'mm10')
  kpAddCytobandLabels(kp)
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
  at <- autotrack(current.track = 1, total.tracks = 5)
  cn.cols <- getCopyNumberColors(colors = "red_blue")
  plotCopyNumberCalls(kp, cn.calls.cmyc, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
  dev.off()
}


# plotting across all chromosome - dividing into two plots
# plot1 = chr1-chr11
pdf('copyNumberPlots/cMYC_CNCopyNumber1to11_.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[1:11], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls.cmyc, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
dev.off()

# plot1 = chr12-chrY
pdf('copyNumberPlots/cMYC_CNCopyNumber12toY_.pdf', width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls.cmyc, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
dev.off()

# adding gene annotations to plots ------------------
# 1. MYCN
gene.symbols <- c(goi$gene)

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", 
               "strand", "start_position", "end_position","gene_biotype"), 
             mart=ensembl)

# subset
genes.new <- genes[,c(2:6)]
# fix gene column name
genes.new$mgi_symbol <- toupper(genes.new$mgi_symbol)

# merge with genes of interest list
genes.new <- merge(goi, genes.new, by.x = 'gene', by.y = 'mgi_symbol', all.x = TRUE)
# fix Chromosome column
genes.new$chromosome_name <- gsub('(.*)','chr\\1', genes.new$chromosome_name)

# convert to GRanges 
genes.new <- na.omit(genes.new)
genes1 <- makeGRangesFromDataFrame(genes.new,
                                   keep.extra.columns=F,
                                   ignore.strand=TRUE,
                                   seqinfo=NULL,
                                   seqnames.field=c("seqnames", "seqname",
                                                    "chromosome", "chrom",
                                                    "chr", "chromosome_name",
                                                    "seqid"),
                                   start.field="start_position",
                                   end.field="end_position",
                                   strand.field="strand",
                                   starts.in.df.are.0based=FALSE)

values(genes1) <- data.frame(gene=genes.new$gene)
seqlevelsStyle(genes1) <- "UCSC"

head(genes)

#kp <- plotKaryotype(genome="hg38")
#kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol)

plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[1:11], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5) 
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)



plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5) 
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)



# 2. c-MYC
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[1:11], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls.cmyc, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)

plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "red_blue")
plotCopyNumberCalls(kp, cn.calls.cmyc, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)

