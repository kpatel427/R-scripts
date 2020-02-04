# visualize MAST results
library(reshape)
library(reshape2)
library(dplyr)
library(tidyverse)
library(stringr)
library(ComplexHeatmap)

mast_res <- motif_refseq
mast_res <- mast_res[!duplicated(mast_res),]


# calculating motif_distances to TSS
mast_res <-  mast_res %>%
  mutate('Distance.to.TSS.motifs' = ifelse(strand == '+', TSS_start - motif_start, motif_end - TSS_start))


# subset
mast_res <-  mast_res[,c(2,12,13)]
mast_res <- mast_res[!duplicated(mast_res),]

# grouping the distances
# 0 = NA
# 1 = within_promoter_region
# 2 = close_to_promoter_region
mast_res$group <- ifelse(mast_res$Distance.to.TSS.motifs < -5000,0,
                         ifelse(mast_res$Distance.to.TSS.motifs < -2000, 1,
                                ifelse(mast_res$Distance.to.TSS.motifs < 200, 2,
                                       ifelse(mast_res$Distance.to.TSS.motifs < 500, 1,0))))

# mast_res$group <- ifelse(mast_res$Distance.to.TSS.motifs < -2000, '',
#                          ifelse(mast_res$Distance.to.TSS.motifs < 200, 'within_promoter_region;',''))

# reshaping dataframe
mast_res <- mast_res %>%
  distinct(gene, MOTIF_NAME, group)

# filtering for duplicated rows after grouping
mast_res <- mast_res %>% 
  group_by(gene, MOTIF_NAME) %>%
  top_n(1, abs(group)) 

# replacing group numbers with string
mast_res$group <- ifelse(mast_res$group == 2, 'close_to_promoter_region;',
                         ifelse(mast_res$group == 1, 'within_promoter_region;', ''))


# spreading data frame
mast_res <- melt(mast_res)
reshaped_MAST <- dcast(mast_res, MOTIF_NAME ~ gene, value.var = 'group')  

# shortening motif names
reshaped_MAST$MOTIF_NAME <-  gsub('-ChIP-Seq\\(.*\\)/Homer','', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Hoxd12(Homeobox)/ChickenMSG-Hoxd12.Flag','Hoxd12(Homeobox)/Hoxd12.Flag', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Prop1(Homeobox)/GHFT1-PROP1.biotin','Prop1(Homeobox)/GHFT1-PROP1', reshaped_MAST$MOTIF_NAME, perl = TRUE)

reshaped_MAST <-  reshaped_MAST %>%
  column_to_rownames(var = 'MOTIF_NAME')

reshaped_MAST[is.na(reshaped_MAST)] <-  ''


# converting to matrix
mat <- as.matrix(reshaped_MAST)

# creating oncoprint
pdf('Oncoprint_21Genes.pdf', width = 15, height = 8)
print(oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
                alter_fun = list(
                  #background = function(x, y, w, h) NULL, # for white background
                  'close_to_promoter_region' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "cadetblue", col = NA)),
                  'within_promoter_region' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "goldenrod2", col = NA))
                ), col = c('close_to_promoter_region' = "cadetblue", 'within_promoter_region' = "goldenrod2"),
                show_column_names = TRUE))
dev.off()





# --------------------------------------------------------------------------------------------------#
mast_res <- motif_refseq
mast_res <- mast_res[!duplicated(mast_res),]
mast_res$motif_present = '1;'


# distance between motif-end and peak-start
mast_res$motif_end <-  as.integer(mast_res$motif_end)
mast_res$peak_start <- as.integer(mast_res$peak_start)
mast_res$motif_peak_dis <- mast_res$motif_end - mast_res$peak_start


# subset
mast_res <-  mast_res[,c(2,6,12,14)]

mast_res <- mast_res %>%
  group_by(gene,seq_name)


mast_res$concat <-  paste0(mast_res$gene, '_', mast_res$seq_name)
mast_res <-  mast_res[,c(5,3,4)]
mast_res <- mast_res[!duplicated(mast_res),]

mast_res$concat <-  gsub('_hg19_ct_UserTrack_3545','', mast_res$concat)


# grouping distances
mast_res$group <- ifelse(mast_res$motif_peak_dis <= 25,'within_25bp;',
                         ifelse(mast_res$motif_peak_dis <=50, 'within_50bp;',
                                ifelse(mast_res$motif_peak_dis <= 75, 'within_75bp;',
                                       ifelse(mast_res$motif_peak_dis <=100, 'within_100bp;',
                                              ifelse(mast_res$motif_peak_dis <=150, 'within_150bp;',
                                                     ifelse(mast_res$motif_peak_dis <= 200, 'within_200bp;','more_than_200bp_away;'))))))


# subset
mast_res <- mast_res[,c(1,2,4)]
# removing duplicate rows
mast_res <- mast_res[!duplicated(mast_res),]

# some of the motifs occur twice in the same track for the same gene
# Error: Each row of output must be identified by a unique combination of keys.
# Keys are shared for 18 rows:
#   * 142, 143
# * 515, 516
# * 49, 50
# * 53, 54
# * 71, 72
# * 75, 76
# * 99, 100
# * 113, 114
# * 274, 275
mast_res <-  mast_res[c(-143,-50,-54,-72,-76,-516,-100,-114,-275),]

reshaped_MAST <- mast_res %>%
  distinct(concat, MOTIF_NAME, group) %>%
  spread(key = 'concat', value = 'group')


# shortening motif names
reshaped_MAST$MOTIF_NAME <-  gsub('-ChIP-Seq\\(.*\\)/Homer','', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Hoxd12(Homeobox)/ChickenMSG-Hoxd12.Flag','Hoxd12(Homeobox)/Hoxd12.Flag', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Prop1(Homeobox)/GHFT1-PROP1.biotin','Prop1(Homeobox)/GHFT1-PROP1', reshaped_MAST$MOTIF_NAME, perl = TRUE)

# converting MOTIF_NAME col to rownames
reshaped_MAST <-  reshaped_MAST %>%
  column_to_rownames(var = 'MOTIF_NAME')

reshaped_MAST[is.na(reshaped_MAST)] <-  ''

# converting to matrix
mat <- as.matrix(reshaped_MAST)

#pdf(paste0(Sys.Date(),'_Oncoprint_39Genes_motifs_peakRegion.pdf'), width = 15, height = 10)
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = list(
            #background = function(x, y, w, h) NULL, # for white background
            'within_25bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "cadetblue", col = NA)),
            'within_50bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "goldenrod2", col = NA)),
            'within_75bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "brown", col = NA)),
            'within_100bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "darkolivegreen", col = NA)),
            'within_150bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "dodgerblue", col = NA)),
            'within_200bp' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "tomato", col = NA)),
            'more_than_200bp_away' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "violetred", col = NA))
          ), col = c('within_25bp' = "cadetblue", 'within_50bp' = "goldenrod2", 'within_75bp' = "brown", 'within_100bp' = "darkolivegreen", 'within_150bp' = "dodgerblue",
                     'within_200bp' = "tomato", 'more_than_200bp_away' = "violetred"),
          show_column_names = TRUE)
#dev.off()

# ----------------------------------------------------------------------------------------------------------------------#

mast_res <- motif_refseq 
# remove duplicate rows
mast_res <- mast_res[!duplicated(mast_res),]

# distance between motif-end and peak-start
mast_res$motif_end <-  as.integer(mast_res$motif_end)
mast_res$peak_start <- as.integer(mast_res$peak_start)
mast_res$motif_start <-  as.integer(mast_res$motif_start)
mast_res$peak_end <- as.integer(mast_res$peak_end)

mast_res <-  mast_res %>%
  mutate('peak_mid' = round((peak_start + peak_end)/2)) %>%
  mutate('motif_mid' = round((motif_start + motif_end)/2))
  
# calculating motif offset from center of the peak
mast_res$motif_offset <- round(mast_res$peak_mid - mast_res$motif_mid)


# subset
mast_res <-  mast_res[,c(2,6,12,15)]

mast_res <- mast_res %>%
  group_by(gene,seq_name)


mast_res$concat <-  paste0(mast_res$gene, '_', mast_res$seq_name)
mast_res <-  mast_res[,c(5,3,4)]
mast_res <- mast_res[!duplicated(mast_res),]

mast_res$concat <-  gsub('_hg19_ct_UserTrack_3545','', mast_res$concat)


# grouping distances
mast_res$group <- ifelse(mast_res$motif_offset < -100,'< -100;',
                         ifelse(mast_res$motif_offset < -50, '-100 < motif < -50;',
                                ifelse(mast_res$motif_offset < -20, '-50 < motif < -20;',
                                       ifelse(mast_res$motif_offset < -10, '-20 < motif < -10;',
                                              ifelse(mast_res$motif_offset < 10, '-10 < motif < 10;',
                                                     ifelse(mast_res$motif_offset < 20, '10 < motif < 20;',
                                                            ifelse(mast_res$motif_offset < 50, '20 < motif < 50;',
                                                                   ifelse(mast_res$motif_offset < 100, '50 < motif < 100;','> 100;'))))))))
              

# subset
mast_res <- mast_res[,c(1,2,4)]
# removing duplicate rows
mast_res <- mast_res[!duplicated(mast_res),]
  
# some of the motifs occur twice in the same track for the same gene
# Error: Each row of output must be identified by a unique combination of keys.
# Keys are shared for 4 rows:
# * 21, 22
# * 8, 9
mast_res <-  mast_res[c(-22,-9),]

  
reshaped_MAST <- mast_res %>%
  distinct(concat, MOTIF_NAME, group) %>%
  spread(key = 'concat', value = 'group')



# shortening motif names
reshaped_MAST$MOTIF_NAME <-  gsub('-ChIP-Seq\\(.*\\)/Homer','', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Hoxd12(Homeobox)/ChickenMSG-Hoxd12.Flag','Hoxd12(Homeobox)/Hoxd12.Flag', reshaped_MAST$MOTIF_NAME, perl = TRUE)
reshaped_MAST$MOTIF_NAME <-  gsub('Prop1(Homeobox)/GHFT1-PROP1.biotin','Prop1(Homeobox)/GHFT1-PROP1', reshaped_MAST$MOTIF_NAME, perl = TRUE)

# converting MOTIF_NAME col to rownames
reshaped_MAST <-  reshaped_MAST %>%
  column_to_rownames(var = 'MOTIF_NAME')

reshaped_MAST[is.na(reshaped_MAST)] <-  ''

# converting to matrix
mat <- as.matrix(reshaped_MAST)

pdf(paste0(Sys.Date(),'_Oncoprint_24Genes_motifs_peakRegion.pdf'), width = 15, height = 10)
print(oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = list(
            #background = function(x, y, w, h) NULL, # for white background
            '< -100' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "cadetblue", col = NA)),
            '-100 < motif < -50' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "darkmagenta", col = NA)),
            '-50 < motif < -20' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "brown", col = NA)),
            '-20 < motif < -10' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "darkolivegreen", col = NA)),
            '-10 < motif < 10' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "dodgerblue", col = NA)),
            '10 < motif < 20' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "tomato", col = NA)),
            '20 < motif < 50' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "goldenrod2", col = NA)),
            '50 < motif < 100' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "darkorange2", col = NA)),
            '> 100' = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "violetred", col = NA))
          ), col = c('< -100' = "cadetblue", 
                     '-100 < motif < -50' = "darkmagenta", 
                     '-50 < motif < -20' = "brown", 
                     '-20 < motif < -10' = "darkolivegreen", 
                     '-10 < motif < 10' = "dodgerblue",
                     '10 < motif < 20' = "tomato", 
                     '20 < motif < 50' = "goldenrod2",
                     '50 < motif < 100' = "darkorange2",
                     '> 100' = "violetred"),
          show_column_names = TRUE))
dev.off()
