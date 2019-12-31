library(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)
library(readxl)
library(ggpubr)
library(pheatmap)
library(dendextend)
library(RColorBrewer)

hm <- broad_expr_ccle_genelist_NBL_merge[,c(1:3,5)]

# create a unique column
hm$concat <- paste0(hm$gene,'_',hm$type)
hm <- hm[,c(1,3,5)]
# remove duplicates
hm <-  hm[!duplicated(hm),]
# reshaping df
hm <- hm %>%
  spread(key = concat, value = ceres_score)

# only 17 NBL cell lines have dependency data
# check!
#unique(intersect(broad_ccle$X, nbl_cl$cellLine))

colnames(hm) <-  gsub('_.*','', colnames(hm), perl = TRUE)
hm <- hm %>%
  column_to_rownames(var = 'cellLine')

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


# column annotation
col_anno <- broad_expr_ccle_genelist_NBL_merge[,c(2,5)]
col_anno <-  col_anno %>%
  column_to_rownames(var = 'gene')
  
# customizing scale
breaksList = seq(-2.5, 0.9, by = 0.5)


# heatmap
p <- pheatmap(hm,
         cutree_cols = 4,
         cluster_cols = my_hclust_gene,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)
