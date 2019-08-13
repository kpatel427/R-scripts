# this function uses pheatmap package to make heatmaps
# Takes input matrix and row annotation
# normalizes matrix values
# outputs heatmap with clustered rows (k=3) with row annotation
plotHeatMap <-  function(mat, gene_anno){
  
  my_gene_col <-  data.frame()
  gene_anno <-  data.frame()
  
  
  # data wrangling
  mat <- mat %>%
    mutate('concat_col' = paste0(Gene.Name,'_', row_number())) %>%
    select(-c(Gene.Name,Annotation)) %>%
    column_to_rownames(var = 'concat_col')
  
  
  # convert to matrix
  mat <- as.matrix(mat)
  
  # normalization
  mat <- t(apply(mat, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  # clustering rows
  my_hclust_gene <- hclust(dist(mat), method = "complete")
  my_gene_col <- as.data.frame(cutree(tree = as.dendrogram(my_hclust_gene), k = 3))
  colnames(my_gene_col)[1] <- 'cluster'
  
  gene_anno <- merge(gene_anno, my_gene_col, by = 'row.names')
  gene_anno <-  gene_anno %>% column_to_rownames(var = 'Row.names')
  
  pheatmap(mat, 
           annotation_row = gene_anno)
  
  #ggsave(filename = paste0(Sys.Date(),'_heatmap_MaxPvalOnly_histone_mycn.pdf'), p, width = 10, height = 10)
  
}
