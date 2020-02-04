# script to create MYCN coverage across cellLines
# setwd("~/KP/ExomeDepth")

library(tidyverse)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(grid)
library(ggpubr)
library(corrplot)
library(ggrepel)


#mycnDepth <- read.delim('mycn_depths.txt', sep = ' ', header = F)
mycnDepth <- read.delim('mycnCov_bedtools.txt', sep = ' ', header = F)
# subset
mycnDepth <- mycnDepth[,c(1,5)]
mycnDepth$V5 <- as.numeric(mycnDepth$V5)
mycnDepth$V2 <- gsub('_.+','',mycnDepth$V1, perl = TRUE)
names(mycnDepth) <- c('name','MYCN_depth', 'CellLine') 


# reading cellLine data
#load('~/KP/RShiny/data/cellLine_mData.RData')

# get mycn status for our cell lines
mycnDepth <- merge(mycnDepth, cellline_mData, by = 'CellLine', all.x = TRUE)
mycnDepth <- mycnDepth[,c(1:4)]
#mycnDepth[is.na(mycnDepth)] <- 'missing_info'
mycnDepth <- na.omit(mycnDepth)

#write.table(mycnDepth, file = 'mycnDepth.txt', row.names = F, col.names = T, sep = '\t', quote = F)


# read average whole exome data
exomCov <- read.delim('exomeCov_bedtools.txt', header = F, sep = ' ')
# subset
exomCov = exomCov[,c(1,5)]
exomCov$V5 <- as.numeric(exomCov$V5)
exomCov$V1 <- gsub('_.+','',exomCov$V1, perl = TRUE)
names(exomCov) <- c('CellLine','exomCov')

exomCov <- merge(exomCov,mycnDepth, by = 'CellLine', all.x = TRUE)
exomCov <- na.omit(exomCov)

exomCov <- exomCov %>%
  group_by(CellLine) %>%
  mutate(mean_MYCN_dept = mean(MYCN_depth)) %>%
  select(-name, -MYCN_depth) %>%
  mutate(mean_exomCov = mean(exomCov)) %>%
  select(-exomCov) %>%
  distinct()


exomCov$MYCN_Status <- ifelse(exomCov$CellLine == 'Kelly', 'Amplified', exomCov$MYCN_Status)

exomCov$coverageX <- exomCov$mean_MYCN_dept/exomCov$mean_exomCov
exomCov$coverageX <- round(exomCov$coverageX)

p <- ggplot(exomCov, aes(x = reorder(CellLine, coverageX), y = mean_exomCov, fill = MYCN_Status)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_point(shape = 1, size = 4, alpha = 0.3) +
  geom_line(aes(group = 1), alpha = 0.3) +
  theme_bw() +
  #guides(colour = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=paste0(as.character(coverageX),'x')), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(title = 'MYCN depth in Exome data across various cell lines', x = 'Cell Lines', y = 'Exome Depth', fill = 'MYCN Status')


ggsave(p, filename = paste0(Sys.Date(),'_mycn_depth_exome.pdf'), width = 8, height = 7)



# reshape
exomCov <- exomCov %>%
  gather(key = 'name', value = 'depth', -c(CellLine,MYCN_Status, coverageX))

exomCov$name <- ifelse(exomCov$name == 'mean_MYCN_dept', 'MYCN depth', 'exome depth')
# convert CellLine names to upper
exomCov$CellLine <- toupper(exomCov$CellLine)


# p <- ggplot(mycnDepth, aes(reorder(mycnDepth$CellLine, avg_depth), mycnDepth$avg_depth, fill = mycnDepth$MYCN_Status)) +
#   geom_bar(stat = 'identity') +
#   labs(title = 'MYCN depth across various cell lines', x = 'Cell Lines', y = 'Depth', fill = 'MYCN Status') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         plot.title = element_text(hjust = 0.5))
# 
# ggsave(p, filename = 'mycn_depth_exome.pdf', width = 8, height = 7)


ggplot(exomCov, aes(reorder(exomCov$CellLine, depth), depth, fill = forcats::fct_rev(name), color = MYCN_Status)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme_bw() +
  scale_fill_manual(values = c('darkcyan','goldenrod2')) +
  geom_text(aes(label=ifelse(name == 'MYCN depth',paste0(as.character(coverageX),'x'),'')), position=position_dodge(width=0.9), vjust = -1, hjust = 0.2, show.legend = FALSE) +
  labs(title = 'MYCN depth in Exome data across various cell lines', x = 'Cell Lines', y = 'Depth', fill = 'Depth', color = 'MYCN Status') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(a, filename = paste0(Sys.Date(),'_stackedBar_mycn_depth_exome.pdf'), width = 10, height = 8)


# ggplot(exomCov, aes(x = reorder(CellLine, depth), y = depth)) +
#   geom_bar(data = exomCov, aes(x = reorder(CellLine, depth), y = depth, fill = MYCN_Status), stat = "identity", position = position_dodge()) +
#   geom_point(aes(group = factor(name), colour = factor(name)), shape = 1, size = 4) +
#   geom_line(aes(group = factor(name), colour = factor(name))) +
#   theme_bw() +
#   guides(colour = FALSE) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         plot.title = element_text(hjust = 0.5)) +
#   geom_text(aes(label=ifelse(factor(name)=='MYCN depth',as.character(coverageX),'')),hjust=1,vjust=1) +
#   labs(title = 'MYCN depth in Exome data across various cell lines', x = 'Cell Lines', y = 'Depth', fill = 'Depth', colour = 'MYCN Coverage')


# ------------------------------------- Correlation of CNA expr vs MYCN coverage(x) -------------------------------------------- 
data <- read.delim('/Volumes/target_nbl_ngs/KP/RShiny/data/datasets_desc.txt', header = T)
load('/Volumes/target_nbl_ngs/KP/RShiny/data/STAR_FPKM_40cells_genes.RData')

# reshape
star.reshape <- STAR_FPKM_40cells_genes %>%
  rownames_to_column(var = 'gene') %>%
  gather(key = 'CellLine', value = 'FPKM', -gene) %>%
  filter(gene == 'MYCN')

names(star.reshape)[3] <- 'MYCN_FPKM'

# merge
exomCov <- merge(exomCov, star.reshape, by = 'CellLine', all.x = T)
# SKNKAN has no expression data in STAR FPKM 40 cell lines dataset
exomCov <- na.omit(exomCov)


# plot correlation
corrPlot <- ggplot(exomCov, aes(coverageX, MYCN_FPKM, color = MYCN_Status)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  theme_bw() +
  xlab('Coverage (x)') + ylab('FPKM') +
  labs(title = 'MYCN', color = 'MYCN Status') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = CellLine)) +
  stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~")))

ggsave(corrPlot, filename = paste0(Sys.Date(),'_correlation_MYCN_expr_depth.pdf'), width = 8, height = 7)



# plotting MYCN coverage and expression together
top = ggplot(exomCov, aes(reorder(CellLine, depth), depth, fill = forcats::fct_rev(name), color = MYCN_Status)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme_bw() +
  ylim(0,50000) +
  guides(color = FALSE) +
  scale_fill_manual(values = c('darkcyan','goldenrod2')) +
  geom_text(aes(label=ifelse(name == 'MYCN depth',paste0(as.character(coverageX),'x'),'')), position=position_dodge(width=0.9), vjust = -1, hjust = 0.2, show.legend = FALSE) +
  labs(title = 'MYCN', x = '', y = 'Depth', fill = 'Depth', color = 'MYCN Status') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5,size = 18),
        plot.margin=unit(c(1,1,-0.5,1),"cm"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))


bottom = ggplot(exomCov, aes(reorder(CellLine, depth), MYCN_FPKM, fill = MYCN_Status)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() +
  labs(title = '', x = 'Cell Lines', y = 'FPKM', fill = 'MYCN Status') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(hjust = 0.5,size = 18),
          plot.margin=unit(c(-0.5,0.8,1.5,1.3),"cm"),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15))

plots <- grid.arrange(top,bottom)
ggsave(plots, filename = paste0(Sys.Date(),'_MYCN_cov_expr.pdf'), width = 15, height = 12)


# ------------------ reading CMYC depths --------------- 
cmyc_depth <- read.delim('cmycCov_bedtools.txt', sep = ' ', header = F)
cmyc_depth <- cmyc_depth[,c(1,5)]
cmyc_depth$V1 <- gsub('_.+','', cmyc_depth$V1)
cmyc_depth$V1 <- toupper(cmyc_depth$V1)
names(cmyc_depth) <- c('CellLine','cmyc_depth')

cmyc_depth <- cmyc_depth %>%
  group_by(CellLine) %>%
  summarise(mean_cmyc_dept = mean(cmyc_depth)) %>%
  distinct()

exomCov <- read.delim('exomeCov_bedtools.txt', header = F, sep = ' ')
# subset
exomCov = exomCov[,c(1,5)]
exomCov$V5 <- as.numeric(exomCov$V5)
exomCov$V1 <- gsub('_.+','',exomCov$V1, perl = TRUE)
names(exomCov) <- c('CellLine','exomCov')
exomCov$CellLine <- toupper(exomCov$CellLine)

exomCov <-  exomCov %>%
  group_by(CellLine) %>%
  summarise(mean_exome_depth = mean(exomCov)) %>%
  distinct()

# merge cmyc + exome
cmyc_depth <- merge(cmyc_depth, exomCov, by = 'CellLine', all.x = T)


cmyc_depth$cmyc_CovX <- cmyc_depth$mean_cmyc_dept/cmyc_depth$mean_exome_depth
names(cmyc_depth) <- c('CellLine','cmyc_depth','exome_depth','coverageX')

write.table(cmyc_depth, file = 'cmyc_coverages.txt', col.names = T, row.names = F, sep = '\t', quote = F)



