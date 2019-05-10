library(dplyr)
library(tidyverse)
library(ggplot2)
library(DT)
library(DESeq2)

# reading patient data
patient <- read.delim("2018-08-22-DiffExp_MYCN_adjp.10.txt", sep = "\t", header = T)
patient$Gene <-  rownames(patient)

# reading cell line data
cellLine <-  read.delim("2018-09-04-CellLineDiffExp_MYCN_adjp.10.txt", sep = "\t", header = T)
cellLine$Gene <-  rownames(cellLine)

# merging for genes found in both
overlapping_genes <-  merge(cellLine,patient, by = "Gene")
overlapping_genes <-  overlapping_genes[,1:7]

# -log10 transforming p values
overlapping_genes$pvalue.log <-  -log10(overlapping_genes$pvalue.x)


overlapping_genes <- overlapping_genes %>%
  mutate(threshold = ifelse(log2FoldChange.x > 0,"Up-regulated", ifelse(log2FoldChange.x < 0 , "Down-regulated", "C")))


ggplot(overlapping_genes, aes(x=log2FoldChange.x, y=-log10(pvalue.x))) +
  geom_point(aes(colour = threshold), size=2.5, alpha = 1/4) +
  scale_colour_manual(values = c("Up-regulated"= "orange", "Down-regulated"="steelblue4",  "C"= "black")) +
  labs(title = "Differentially expressed genes in Patient and CellLine data", x = "log2FoldChange", y="-log10(Pvalue)") +
  theme_bw() +
  theme(legend.title = element_blank(),
      legend.text = element_text(size=12),
      axis.title.x = element_text(size=15),
      axis.title.y = element_text(size=15),
      axis.text = element_text(size=12))


ggsave(filename = paste0(Sys.Date(),"-volcanoPlot-TARGET_CellLine.pdf"),plot = last_plot(),width = 8, height=6, device = "pdf")


#----------------------------------------------------------------------------------------#
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhancedVolcano(overlapping_genes,
                lab = overlapping_genes$Gene,
                x = "log2FoldChange.x",
                y = "pvalue.x",
                transcriptPointSize = 1,
                transcriptLabSize = 3,
                FCcutoff = 0,
                legend = c('NS','NS','P value','p value & Log (base 2) fold-change'),
                legendPosition = 'top')
