# script to create a plot showing MYCN depth in Exome data across various cellLines
# creates a barplot overlayed with a lineplot showing different MYCN coverages
# setwd("~/KP/ExomeDepth")

library(tidyverse)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(grid)


#mycnDepth <- read.delim('mycn_depths.txt', sep = ' ', header = F)
mycnDepth <- read.delim('mycnCov_bedtools.txt', sep = ' ', header = F)
# subset
mycnDepth <- mycnDepth[,c(1,5)]
mycnDepth$V5 <- as.numeric(mycnDepth$V5)
mycnDepth$V2 <- gsub('_.+','',mycnDepth$V1, perl = TRUE)
names(mycnDepth) <- c('name','MYCN_depth', 'CellLine') 


# reading cellLine data
#load('~/KP/cellLine_mData.RData')

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
  labs(title = 'MYCN depth in Exome data across various cell lines', x = 'Cell Lines', y = 'Depth', fill = 'MYCN Status')


ggsave(p, filename = paste0(Sys.Date(),'_mycn_depth_exome.pdf'), width = 8, height = 7)

