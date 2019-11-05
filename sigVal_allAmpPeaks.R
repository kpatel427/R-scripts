# 1. Keep only one peak (max signal value) for each gene.
# 2. Merge by genes for all amp lines
# 3. select genes which have peaks in all amp lines
# 4. compute median/mean signal value
# "~/KP/H3K27me3/masterFile_histoneMarks_MYCN/"

library(data.table)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# function definition
# takes cellLine name, histone/TF mark and working dir as arguments
peakMerge <- function(cellLine, histone, wd){
    # reading peakFile
    if(histone == 'NMYC'){
      df.name <- read.delim(paste0(wd,cellLine,'-',histone,'.macs2.SPMR.filtered.narrowPeak'), header = F)
    } else{
      df.name <- read.delim(paste0(wd,cellLine,'-',histone,'.macs2.SPMR.filtered.broadPeak'), header = F)
    }
    # subset peakflile = chr,start,end,peakID,signalValue
    df.name <- df.name[,c(1:4,7)]
    names(df.name) <-  c('Chr','Start','End','peakID','SignalValue')
   
    if(histone == 'NMYC'){
      df.name.anno <- read.delim(paste0(wd,cellLine,'-',histone,'.macs2.SPMR.filtered.narrowPeak.anno.txt'), header = T)
    } else {
      df.name.anno <- read.delim(paste0(wd,cellLine,'-',histone,'.macs2.SPMR.filtered.broadPeak.anno.txt'), header = T)    
    }  
      # subset anno file = chr,end,peakID,geneName
    df.name.anno <-  df.name.anno[,c(1,2,4,16)]
    names(df.name.anno) <-  c('peakID','Chr','End','Gene.Name')
    
    # merging peakInfo with metadata
    df.name.merge <- merge(df.name.anno,df.name, by = c('peakID','Chr','End'))
    
    # getting max signalValue for each gene
    final.df <- df.name.merge %>%
      select(Gene.Name, peakID, SignalValue) %>%
      group_by(Gene.Name) %>%
      #tally()
      mutate('Max_SignalValue' = max(SignalValue)) %>%
      select(-SignalValue, -peakID) %>%
      distinct()
  
  names(final.df)[2] <- paste0('Max_SignalValue',cellLine,'_',histone)
  return(final.df)
}


# function call for all histone/NMYC marks
cl <- c('COGN415','KELLY','LAN5','NB1643','NB69')
marks <- c('H3K27Ac','H3K27me3','H3K4me1','H3K4me3','NMYC')

for(y in marks){
  work_dir <- paste0('~/KP/H3K27me3/masterFile_histoneMarks_MYCN/',y,'_overlap/')
  for(x in cl){
    
    name <- paste0(x,'_',y)
    #print(name)
    assign(name, peakMerge(cellLine = x, histone = y, wd = work_dir))

  }
}


COG <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(COGN415_H3K27Ac,
                                                                                COGN415_H3K27me3,
                                                                                COGN415_H3K4me1,
                                                                                COGN415_H3K4me3,
                                                                                COGN415_NMYC))


KELLY <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(KELLY_H3K27Ac,
                                                                                  KELLY_H3K27me3,
                                                                                  KELLY_H3K4me1,
                                                                                  KELLY_H3K4me3,
                                                                                  KELLY_NMYC))

LAN5 <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(LAN5_H3K27Ac,
                                                                                 LAN5_H3K27me3,
                                                                                 LAN5_H3K4me1,
                                                                                 LAN5_H3K4me3,
                                                                                 LAN5_NMYC))

NB1643 <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(NB1643_H3K27Ac,
                                                                                   NB1643_H3K27me3,
                                                                                   NB1643_H3K4me1,
                                                                                   NB1643_H3K4me3,
                                                                                   NB1643_NMYC))


NB69 <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(NB69_H3K27Ac,
                                                                                 NB69_H3K27me3,
                                                                                 NB69_H3K4me1,
                                                                                 NB69_H3K4me3,
                                                                                 NB69_NMYC))


# read 46 genes
genes46 <- read.delim('~/KP/H3K27me3/motif-analysis/46_candidate\ genes/both-immune-inter-HR.txt',
                      header = F)

COG_46 <- merge(COG, genes46, by.x = 'Gene.Name', by.y = 'V1')
KELLY_46 <- merge(KELLY, genes46, by.x = 'Gene.Name', by.y = 'V1')
LAN5_46 <- merge(LAN5, genes46, by.x = 'Gene.Name', by.y = 'V1')
NB1643_46 <- merge(NB1643, genes46, by.x = 'Gene.Name', by.y = 'V1')
NB69_46 <- merge(NB69, genes46, by.x = 'Gene.Name', by.y = 'V1', all.y = TRUE)


all_merge <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(COG_46, KELLY_46, NB1643_46, LAN5_46, NB69_46))
all_merge[is.na(all_merge)] <- 0

#write.table(all_merge, file = paste0(Sys.Date(),'_signalVal_allMerge_46Genes.txt'), sep = '\t', col.names = T, quote = F)

# for plotting heatmap
all_merge <- all_merge %>%
  column_to_rownames(var = 'Gene.Name')
mat <- as.matrix(all_merge)

col_names <- as.data.frame(colnames(all_merge))
names(col_names)[1] <- 'col.names'

col_anno <- read.delim('col_annotation.txt', header = T)
col_anno <-  col_anno %>%
  column_to_rownames(var = 'col.names')
col_anno$cellLine <-  as.factor(col_anno$cellLine)
col_anno$mark <-  as.factor(col_anno$mark)

# row annotation
oy_grp <- read.delim('~/KP/H3K27me3/motif-analysis/46_candidate genes/PCA/orangeYellowGroup_PCA.txt', header = F)
oy_grp$group <- 'Orange-Yellow'

bg_grp <- read.delim('~/KP/H3K27me3/motif-analysis/46_candidate genes/PCA/greebBlueGroup_PCA.txt', header = F)
bg_grp$group <- 'Blue-Green'

row_anno <- rbind(oy_grp,bg_grp)

row_anno <-  merge(genes46, row_anno, by = 'V1', all.x = TRUE)
row_anno[is.na(row_anno)] <-  'outliers'

row_anno <-  row_anno %>%
  column_to_rownames(var = 'V1')
row_anno$group <- as.factor(row_anno$group)

anno_colors <- list(
  cellLine = c(COGN415 = '#8c662e',KELLY = '#ce9643', LAN5 = '#f7efb3',NB1643 = '#a1d379',NB69 = '#3b712d'),
  mark = c(H3K27Ac = '#cc4125',H3K27me3 = '#5a5b9f',H3K4me1 = '#f6b265',H3K4me3 = '#ffd966', MYCN = '#93c47d'),
  group = c(`Blue-Green` = 'royalblue',`Orange-Yellow` = 'orange', `outliers` = 'grey')
)

hm <- pheatmap(mat, 
         annotation_col = col_anno, 
         annotation_row = row_anno, 
         annotation_colors = anno_colors,
         color = colorRampPalette(c("white", "red"))(30),
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = F)
ggsave(hm, file = paste0(Sys.Date(),'_46Genes_signalVal.pdf'), height = 10, width = 12)




# seperating into smaller dfs as per the marks
H3K27Ac <- all_merge[,c(1,2,7,12,17)]
H3K27me3 <- all_merge[,c(1,3,8,13,18)]
H3K4me1 <- all_merge[,c(1,4,9,14,19)]
H3K4me3 <- all_merge[,c(1,5,10,15,20)]
NMYC <- all_merge[,c(1,6,11,16,21)]

H3K27Ac <-  H3K27Ac %>%
  column_to_rownames(var = 'Gene.Name')
H3K27Ac$Max_H3K27Ac <- apply(H3K27Ac, 1, 'max')
H3K27Ac <-  H3K27Ac %>%
  rownames_to_column(var = 'Gene.Name')
H3K27Ac <- H3K27Ac[,c(1,6)]

H3K27me3 <-  H3K27me3 %>%
  column_to_rownames(var = 'Gene.Name')
H3K27me3$MAx_H3K27me3 <- apply(H3K27me3, 1, 'max')
H3K27me3 <-  H3K27me3 %>%
  rownames_to_column(var = 'Gene.Name')
H3K27me3 <- H3K27me3[,c(1,6)]


H3K4me1 <- H3K4me1 %>%
  column_to_rownames(var = 'Gene.Name')
H3K4me1$Max_H3K4me1 <- apply(H3K4me1, 1, 'max')
H3K4me1 <- H3K4me1 %>%
  rownames_to_column(var = 'Gene.Name')
H3K4me1 <-  H3K4me1[,c(1,6)]

H3K4me3 <-  H3K4me3 %>%
  column_to_rownames(var = 'Gene.Name')
H3K4me3$Max_H3K4me3 <-  apply(H3K4me3, 1, 'max')
H3K4me3 <-  H3K4me3 %>%
  rownames_to_column(var = 'Gene.Name')
H3K4me3 <-  H3K4me3[,c(1,6)]


NMYC <-  NMYC %>%
  column_to_rownames(var = 'Gene.Name')
NMYC$Max_NMYC <-  apply(NMYC, 1, 'max')
NMYC <-  NMYC %>%
  rownames_to_column(var = 'Gene.Name')
NMYC <-  NMYC[,c(1,6)]


final <- Reduce(function(x, y) merge(x, y, by = c("Gene.Name"), all = TRUE), list(NMYC, H3K27Ac, H3K27me3, H3K4me1, H3K4me3))

#write.table(final, file = paste0(Sys.Date(),'_signalVal_peaks_46Genes.txt'), sep = '\t', col.names = T, quote = F)




