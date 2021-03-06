# script to generate CopyNumber plots for MDX project - for MYCN & c-MYC MDXs across phases
# setwd("/Volumes/target_nbl_ngs/KP/MDX_Minu/cnvkit/withoutControl")
library(CopyNumberPlots)
library(gridExtra)
library(biomaRt)


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
cn.calls1 <- list("MYCN-0685-P0"=MYCN0685P0, "MYCN-3927-P1"=MYCN3927P1,"MYCN-2900-P2"=MYCN2900P2, "MYCN-6792-P1"=MYCN6792P1, "MYCN-289-P2"=MYCN289P2)

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
pdf(paste0(Sys.Date(),'_copyNumberPlots/MYCN_CNCopyNumber1to11_.pdf'), width = 10, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[1:11], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "blue_red")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.5,1.5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
dev.off()

# plot1 = chr12-chrY
pdf(paste0(Sys.Date(),'_copyNumberPlots/MYCN_CNCopyNumber12toY_.pdf'), width = 10, height = 9, onefile = F)
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
# cn.cols <- getCopyNumberColors(colors = "green_orange_red")
# #cn.cols <- getCopyNumberColors(c("red", "#FFAAAA", "gray", "#AAAAFF", "#5555FF", "#0000FF"))
# plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5)
# legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,1.5), xpd = TRUE,
#        horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.5, x.intersp=0.5, xjust=0, yjust=0)
# 


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
# ...1. MYCN ----------
goi <- read.delim('/Volumes/target_nbl_ngs/KP/MDX_Minu/goi.txt', header = F)
names(goi)[1] <- 'gene'
goi[25,'gene'] <- 'SDHD'
goi[26,'gene'] <- 'SDHB'
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
cn.cols <- getCopyNumberColors(colors = "gree_orange_red")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5) 
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)



plot.new()
kp <- plotKaryotype(genome = 'mm10', plot.type = 1, chromosomes = paste('chr',chrs[12:21], sep = ''))
kpAddCytobandLabels(kp)
cn.cols <- getCopyNumberColors(colors = "green_orange_red")
plotCopyNumberCalls(kp, cn.calls1, cn.colors = cn.cols, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.5) 
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
legend("topright", legend=names(cn.cols), fill = cn.cols, ncol=length(cn.cols), inset = c(-0.2,5), xpd = TRUE, 
       horiz = TRUE, col = c(cn.cols), bty = "n", cex = 0.4, x.intersp=0.5, xjust=0, yjust=0)



# .....plotting for specific chrs ------------
pdf('MYCN_CNCopyNumber_chr4.pdf', width = 12, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosome = 'chr4', genome = 'mm10', zoom = toGRanges(data.frame(chr='chr4', start=140e6, end=142e6)))
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
at <- autotrack(current.track = 1, total.tracks = 5)
cn.cols <- getCopyNumberColors(colors = "green_orange_red")
plotCopyNumberCalls(kp, cn.calls1, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
kpRect(kp, chr="chr4", x0=140688514, x1=140706504, y0=0, y1=1, col=transparent("#ccb3ff"), data.panel="all", border=NA)

dev.off()

# coordinates for SDHB from biomart and IGV is different

pdf('MYCN_CNCopyNumber_chr9.pdf', width = 12, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosome = 'chr9', genome = 'mm10')
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
at <- autotrack(current.track = 1, total.tracks = 5)
cn.cols <- getCopyNumberColors(colors = "green_orange_red")
plotCopyNumberCalls(kp, cn.calls1, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
dev.off()



# ...2. c-MYC --------------
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




# .....plotting for specific chrs ------------
pdf('cMYC_CNCopyNumber_chr4.pdf', width = 15, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosome = 'chr4', genome = 'mm10')
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
at <- autotrack(current.track = 1, total.tracks = 5)
cn.cols <- getCopyNumberColors(colors = "green_orange_red")
plotCopyNumberCalls(kp, cn.calls.cmyc, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
dev.off()

pdf('cMYC_CNCopyNumber_chr9.pdf', width = 15, height = 9, onefile = F)
plot.new()
kp <- plotKaryotype(chromosome = 'chr9', genome = 'mm10')
kpAddCytobandLabels(kp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1, minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")   
at <- autotrack(current.track = 1, total.tracks = 5)
cn.cols <- getCopyNumberColors(colors = "green_orange_red")
plotCopyNumberCalls(kp, cn.calls.cmyc, ymin = -1, ymax=10, r0=0, r1=1, loh.height = 10, label.cex = 0.8, cn.colors = cn.cols)
kpPlotMarkers(kp, data=genes1, labels=genes1$gene,text.orientation = "horizontal",
              r1=0.5, cex=0.8, adjust.label.position = FALSE)
dev.off()



