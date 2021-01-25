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
