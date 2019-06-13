library(tidyverse)
library(dplyr)
library(data.table)
library(pheatmap)
library(reshape)
library(dendsort)
library(heatmaply)

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

# widening the data
temp <- dcast(test, gene ~ stripped_cell_line_name)
temp[is.na(temp)] <- 0

# ---- Plotting heatmap ----
# 1. input data should be a data matrix
# 2. data must be numeric
# 3. row.names should be present

# working sample
hm <- temp[1:100,]

# creating column annotation
cellLine_col <- data.frame(status)
row.names(cellLine_col) <- cellLine_col[,1]
cellLine_col$stripped_cell_line_name <-  NULL

# removing genes with less than 1 log2TPM expression in all cell lines
#hm <- hm[rowSums(hm < 1) ==17 , , drop = FALSE]

row.names(hm) <- hm$gene
hm$gene <-  NULL
hm <-  data.matrix(hm)
pheatmap(hm,
         show_rownames=T,
         cluster_cols= T,
         cluster_rows=T,
         annotation_col = cellLine_col,
         fontsize = 8)
