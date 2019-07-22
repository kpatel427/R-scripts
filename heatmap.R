library(tidyverse)
library(dplyr)
library(data.table)
library(pheatmap)
library(dendextend)

# ---- Reading expression data ----
expr <- read.delim("CCLE_expression_full.csv", sep = ",", header = T)

# ---- gather expr data ----
expr<- expr %>%
  gather(key = "gene", value = "log2TPM", -(X))

# ---- Replace ENGS string with gene name to just get gene names ----
expr$gene <-  gsub("..ENSG.*", "", expr$gene, perl = T)
test.expr <-  expr[1:100,]


# ---- Reading metadata ----
metadata <-  read.delim("sample_meta.csv", sep =",", header =T)
NBL.cellines <- subset(metadata, metadata$disease_sutype == "neuroblastoma", select = c("DepMap_ID","stripped_cell_line_name","disease"))

# ---- Merging NBL.cellines with expr data; merge.expr only contains expr data for NBL cell lines ----
merge.expr <-  merge(NBL.cellines, expr, by.x = 'DepMap_ID', by.y = 'X')


# ---- CRISPR data ----
achilles.data <-  read.delim("Achilles_gene_effect.csv", sep=",", header = T)
crispr <- achilles.data %>%
  gather(key = "gene", value = "CERES.score", -(X))
crispr$gene <-  gsub("\\..[0-9].*", "", crispr$gene, perl = T)

# ---- get MYCN CERES scores for all cell lines ----
mycn.achilles <- subset(achilles.data, select = c(X,grep("MYCN", colnames(achilles.data))))

# ---- QC: since number of NBL cell lines in expr data matrix does not equal to cell lines in crispr data matrix ----
# there are total 32 NBL cell lines
NBL.cellines <-  merge(NBL.cellines, mycn.achilles, by.x = "DepMap_ID", by.y = "X")
setdiff(unique(NBL.cellines$DepMap_ID), mycn.achilles$X)
# there are 15 cell lines that do not have crispr MYCN knockout data
# [1] "ACH-000136" "ACH-000149" "ACH-000203" "ACH-000345" "ACH-000446" "ACH-001188" "ACH-001338" "ACH-001344" "ACH-001355" "ACH-001366" "ACH-001548" "ACH-001603"
# [13] "ACH-001674" "ACH-001716" "ACH-002083"


# ---- splitting cell lines into MYCN depleted cell lines and MYCN enriched cell lines ----
mycn.depleted <-  subset(mycn.achilles,  mycn.achilles$MYCN..4613. < 0, select = c(colnames(mycn.achilles)))
mycn.enriched <-  subset(mycn.achilles,  mycn.achilles$MYCN..4613. > 0, select = c(colnames(mycn.achilles)))


# ---- merging merge.expr data with crispr data ----
mycn.depleted.expr <- merge(mycn.depleted, merge.expr, by.x = "X", by.y = "DepMap_ID")
mycn.depleted.expr$status <- "depleted"
mycn.enriched.expr <- merge(mycn.enriched, merge.expr, by.x = "X", by.y = "DepMap_ID")
mycn.enriched.expr$status <- "enriched"

# ---- cbind depleted and enriched data ----
final.merge <- rbind(mycn.depleted.expr, mycn.enriched.expr)
setnames(final.merge,"MYCN..4613.","CERES.score")

# remove genes withh 0 TPM expression
final.merge <- final.merge %>%
  filter(final.merge$log2TPM !=0)

# filter duplicated rows
final.merge <- final.merge[!duplicated(final.merge),]

# get cell lines and their status into seperate data frame
status <-  final.merge[,c(3,7)]
status <- status[!duplicated(status),]

ggplot(final.merge, aes(x = X, y = log2TPM, fill = status)) +
  geom_boxplot(aes(fill = factor(status))) +
  theme_bw()

# ---- preparing data for heatmap ----
test <-  final.merge[,c(3,5,7,6)]

# QC: finding duplicates in test
test.dups <-  test[,1:2]
dups <- test.dups[duplicated(test.dups),]
dups <-  dups[!duplicated(dups$gene),]

# Keeping only rows with max expression for duplicate genes
test <- test %>%
  group_by(stripped_cell_line_name,gene) %>%
  summarize(max.log2TPM = max(log2TPM))

# filtering out genes with less than 2.321928 expression i.e. less than 5 TPM expression (log2(5))
test <-  test %>%
  filter(max.log2TPM > 2.321928)

# widening the data
temp <- dcast(test, gene ~ stripped_cell_line_name)
temp[is.na(temp)] <- 0

# ---- Plotting heatmap ----
# 1. input data should be a data matrix
# 2. data must be numeric
# 3. row.names should be present

# working sample
hm <- temp[100:200,]

# creating column annotation
cellLine_col <- data.frame(status)
row.names(cellLine_col) <- cellLine_col[,1]
cellLine_col$stripped_cell_line_name <-  NULL

row.names(hm) <- hm$gene
hm$gene <-  NULL

# clustering
cluster = kmeans(hm, 4)
cluster.genes <-  as.data.frame(cluster$cluster)
setnames(cluster.genes, "cluster$cluster","cluster")

hm <-  data.matrix(hm)

pheatmap(hm,
         show_rownames=T,
         cluster_cols= T,
         cluster_rows= T,
         annotation_col = cellLine_col,
         annotation_row = cluster.genes,
         fontsize = 8)


# working sample ends #

# Plotting heatmap for the actual dataframe
row.names(temp) <-  temp$gene
temp$gene <-  NULL

# clustering
cluster = kmeans(temp, 3)
cluster.genes <-  as.data.frame(cluster$cluster)
setnames(cluster.genes, "cluster$cluster","cluster")


temp <-  data.matrix(temp)
pheatmap(temp,
         show_rownames=T,
         cluster_cols= T,
         cluster_rows=T,
         annotation_col = cellLine_col,
         annotation_row = cluster.genes,
         fontsize = 8)

# ---- get genes from clusters 2 & 3 ----
filter_cluster <- cluster.genes %>%
  rownames_to_column(var = "gene")

gene.list <- subset(filter_cluster, filter_cluster$cluster != 1, select = c(colnames(filter_cluster)))

#write.table(gene.list$gene, file = paste0(Sys.Date(), "-gene-list.txt"), quote = F, row.names = F, col.names = F)

# ---- overlap gene list withh repressed and final list ----
repressed.genes <- read.delim("~/KP/differentially-exp-gene-lists-MYCN/repressed_genes.txt", header = F, sep = "\t")

overlap.genes <- as.data.frame(unique(intersect(repressed.genes$V1, gene.list$gene)))
#write.table(overlap.genes, file = paste0(Sys.Date(), "-overlap-with-repressed-genes.txt"), quote = F, row.names = F, col.names = F)

