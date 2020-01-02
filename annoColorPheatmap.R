
library(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggpubr)
library(pheatmap)
library(dendextend)
library(RColorBrewer)

# Plotting Heatmap --------
hm <- broad_expr_ccle_genelist_NBL_merge[,c(1:3,5)]
# create a unique column
hm$concat <- paste0(hm$gene,'_',hm$type)
hm <- hm[,c(1,3,5)]
# remove duplicates
hm <-  hm[!duplicated(hm),]
# reshaping df
hm <- hm %>%
  spread(key = concat, value = ceres_score)

# fixing column names
colnames(hm) <-  gsub('_.*','', colnames(hm), perl = TRUE)
hm <- hm %>%
  column_to_rownames(var = 'cellLine')

# replacing NAs with 0
hm[is.na(hm)] <-  0

# clustering
my_hclust_gene <- hclust(dist(t(hm)), method = "complete")
#plot(my_hclust_gene)
# plotting dendogram
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = FALSE)

# getting cluster info
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4)


cluster_info <- data.frame(gene = names(my_gene_col),
                   cluster = as.data.frame(my_gene_col))

names(cluster_info)[2] <- 'cluster'

# writing cluster info to a table
write.table(cluster_info, file = 'cluster_info_genes.txt', col.names = T, row.names = F, sep = '\t')


# customizing scale
breaksList = seq(-2.5, 0.9, by = 0.5)


# add column annotations
col_anno <- broad_expr_ccle_genelist_NBL_merge[,c(2,5)]
# remove duplicates
col_anno <-  col_anno[!duplicated(col_anno),]

# If present = 1
col_anno$n <- 1

# spread data
col_anno <- col_anno %>%
  spread(key = 'type', value = 'n')

# replacing NAs with 0
col_anno[is.na(col_anno)] <- 0

# convert column to rownames
col_anno <-  col_anno %>%
  column_to_rownames(var = 'gene')


# change annotation colors
col_anno$mb0 <-  as.factor(col_anno$mb0)
col_anno$mb1 <-  as.factor(col_anno$mb1)
col_anno$mb3b <-  as.factor(col_anno$mb3b)
col_anno$mb4 <-  as.factor(col_anno$mb4)

anno_colors <- list(
  mb0 = c('0' = 'azure', '1' = 'royalblue'),
  mb1 = c('0' = 'azure', '1' = 'brown1'),
  mb3b = c('0' = 'azure', '1' = 'darkolivegreen4'),
  mb4 = c('0' = 'azure', '1' = 'orange'))

# heatmap
pheatmap(hm,
         cutree_cols = 4,
         annotation_col = col_anno,
         annotation_colors = anno_colors,
         cluster_cols = my_hclust_gene,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         fontsize_col = 8)

ggsave(filename = paste0(Sys.Date(),'_heatmap_ceres_scores_NBL.pdf'), p, width = 15, height = 10)
