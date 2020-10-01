# script to perform GSEA
# 1. all genes from linear regression (invivo) group - ranked according to p-value
# 2. DEG list ranked according to their fold-changes
# setwd("~/KP/GSEA")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("fgsea")
library(fgsea)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggridges)


# DEG invivo ------------------------------------
# Using logFC as ranks
DEG_invivo <- read.delim('../2020-07-15_DEG_invivo_allGenes.txt', header = T)
DEG_invivo <- na.omit(DEG_invivo)
DEG_invivo <- DEG_invivo[,c(1,4)]


# ...convert gene symbol to EntrezID -----------
hs <- org.Hs.eg.db
my.symbols <- c(as.character(DEG_invivo$gene))
df <- AnnotationDbi::select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")


DEG_invivo <- merge(DEG_invivo, df, by.x = 'gene', by.y = 'SYMBOL', all.x = TRUE)
DEG_invivo <- na.omit(DEG_invivo)
DEG_invivo <- DEG_invivo[is.finite(DEG_invivo$log2FC_sensitive_resistant),]



# ...ranked genelist ------------
rankGenes <- DEG_invivo$log2FC_sensitive_resistant
names(rankGenes) <- DEG_invivo$ENTREZID

# plot ranks 
barplot(sort(rankGenes, decreasing = T))

# ...Perform GSEA -----------------------
# we would want to use reactome pathways ----------
source("https://bioconductor.org/biocLite.R")
biocLite("reactome.db")

#my_pathways <- reactomePathways(names(rankGenes))
# using hallmark pathways
h_pathwyas_invivo <- gmtPathways("/Volumes/target_nbl_ngs/KP/miscellaneous/h.all.v7.2.entrez.gmt")


fgseaRes <- fgsea(pathways = h_pathwyas_invivo, 
                  stats = rankGenes,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)

# ...Looking at top 10 results -----------------
head(fgseaRes[order(padj, -abs(NES)), ], n=10)

# ...Tidy the results ---------
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

write.table(fgseaRes[,c(1:5,7)], file = paste0(Sys.Date(),'_GSEA_DEG_invivo_hallmarkPathways.txt'), col.names = T, row.names = F, quote = T, sep = '\t')


# ...enrichment plot ------------
plotEnrichment(my_pathways[["DNA Damage/Telomere Stress Induced Senescence"]], rankGenes)

# ...Plot GSEA table ----------------
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
deg_reactome_table <- plotGseaTable(my_pathways[topPathways], rankGenes, fgseaRes, 
              gseaParam = 0.5)
deg_hallmark_table <- plotGseaTable(h_pathwyas_invivo[topPathways], rankGenes, fgseaRes, 
              gseaParam = 0.5)

ggsave(deg_reactome_table, filename = 'DEG_reactome_table.pdf', width = 20, height = 10)
ggsave(deg_hallmark_table, filename = 'DEG_hallmark_table.pdf', width = 20, height = 10)


# ...-----
# invivo LR - ALL genes ------------------
# using slope as ranks
invivo.lr <- read.delim('2020-01-21_linearRegression_PFS_expression_invivo_PDX_ALLgenes.txt', header = T)
invivo.lr <- invivo.lr[invivo.lr$term == 'pdx_FPKM',]
invivo.lr <- invivo.lr[,c(1,3)]
names(invivo.lr)[2] <- 'slope'
invivo.lr.ordered <- invivo.lr[order(invivo.lr$slope, decreasing = FALSE),]

# ...convert gene symbol to EntrezID -----------
hs <- org.Hs.eg.db
my.symbols.invivo <- c(as.character(invivo.lr.ordered$gene))
df.invivo <- AnnotationDbi::select(hs, 
                            keys = my.symbols.invivo,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")

# merge entrezIDs with invivo.lr.ordered
invivo.lr.ordered <- merge(invivo.lr.ordered, df.invivo, by.x = 'gene', by.y = 'SYMBOL', all.x = TRUE)
invivo.lr.ordered <- na.omit(invivo.lr.ordered)
invivo.lr.ordered <- invivo.lr.ordered[order(invivo.lr.ordered$slope),]

# ...ranked genelist ------------
rankGenes.invivo <- invivo.lr.ordered$slope
names(rankGenes.invivo) <- invivo.lr.ordered$ENTREZID

# plot ranks 
barplot(sort(rankGenes.invivo, decreasing = T))

# ...Perform GSEA -----------------------
#my_pathways_invivo <- reactomePathways(names(rankGenes.invivo))
# using hallmark pathways
h_pathwyas_invivo <- gmtPathways("/Volumes/target_nbl_ngs/KP/miscellaneous/h.all.v7.2.entrez.gmt")



fgseaResInvivo <- fgsea(pathways = h_pathwyas_invivo, 
                  stats = rankGenes.invivo,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)

head(fgseaResInvivo[order(pval, -abs(NES)), ], n=10)

# ...Tidy the results ---------
fgseaResTidyInvivo <- fgseaResInvivo %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidyInvivo %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

write.table(fgseaResInvivo[,c(1:5,7)], file = paste0(Sys.Date(),'_GSEA_linearRegression_invivo_hallmarkPathways.txt'), col.names = T, row.names = F, quote = T, sep = '\t')



# ...Plot GSEA table ----------------
topPathwaysUp <- fgseaResInvivo[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaResInvivo[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# lr_reactome_table <- plotGseaTable(my_pathways[topPathways], rankGenes.invivo, fgseaResInvivo, 
#               gseaParam = 0.5)
lr_hallmark_table <- plotGseaTable(h_pathwyas_invivo[topPathways], rankGenes.invivo, fgseaResInvivo, 
                                   gseaParam = 0.5)

ggsave(lr_hallmark_table, filename = 'LR_hallmark_table.pdf', width = 20, height = 10)



# ...-----
# PLOT DATA -----------------
deg_reactome <- read.delim('2020-10-01_GSEA_DEG_invivo_reactomePathways.txt')
deg_hallmark <- read.delim('2020-10-01_GSEA_DEG_invivo_hallmarkPathways.txt')
lr_reactome <- read.delim('2020-10-01_GSEA_linearRegression_invivo_reactomePathways.txt')
lr_hallmark <- read.delim('2020-10-01_GSEA_linearRegression_invivo_hallmarkPathways.txt')

# deg_hallmark
p1 <- ggplot(deg_hallmark, aes(reorder(pathway, NES), NES, fill = pval <0.05)) +
  #geom_point(aes(size = size)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = '') +
  theme_bw()

ggsave(p1, filename = 'barplot_DEG_hallmarkPathway.pdf', width = 10, height = 8)

# deg_reactome
ggplot(deg_reactome, aes(reorder(pathway, NES), NES, color = pval <0.05)) +
  geom_point(aes(size = size)) +
  coord_flip() +
  theme_bw()

deg_reactome <- deg_reactome[deg_reactome$pval < 0.05,]
p2 <- ggplot(deg_reactome, aes(reorder(pathway, NES), NES, fill = pval <0.05)) +
  geom_bar(stat = 'identity') +
  labs(x = '') +
  coord_flip() +
  scale_fill_manual(values = c("#00BFC4")) +
  theme_bw()
ggsave(p2, filename = 'barplot_DEG_reactomePathway.pdf', width = 10, height = 8)


# lr_hallmark
p3 <- ggplot(lr_hallmark, aes(reorder(pathway, NES), NES, fill = pval <0.05)) +
  geom_bar(stat = 'identity') +
  labs(x = '') +
  coord_flip() +
  theme_bw()
ggsave(p3, filename = 'barplot_LR_hallmarkPathway.pdf', width = 10, height = 8)


# lr_reactome
ggplot(lr_reactome, aes(reorder(pathway, NES), NES, color = pval <0.05)) +
  geom_point(aes(size = size)) +
  coord_flip() +
  theme_bw()

lr_reactome <- lr_reactome[lr_reactome$pval<0.05,]
p4 <- ggplot(lr_reactome, aes(reorder(pathway, NES), NES, fill = pval <0.05)) +
  geom_bar(stat = 'identity') +
  labs(x = '') +
  coord_flip() +
  scale_fill_manual(values = c("#00BFC4")) +
  theme_bw()
ggsave(p4, filename = 'barplot_LR_reactomePathway.pdf', width = 10, height = 8)
