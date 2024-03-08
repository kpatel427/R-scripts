# script to wrangle data and plot correlation
library(data.table)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(ggforce)
library(ggplot2)


# function to wrangle data
dataWrangle <- function(data, mData){
  df.name <- data %>%
    rownames_to_column(var = 'gene') %>%
    gather(key = 'sample', value = 'FPKM', -c(gene)) %>%
    group_by(gene)
  #mutate('mean_FPKM' = mean(FPKM)) %>%
  #mutate('median_FPKM' = median(FPKM))

  # merge with metadata
  df.name <-  merge(df.name, mData, by.x = 'sample', by.y = 'Title')
  names(df.name)[4] <- 'MYCN_Status'
  
  # get expression values for genes in genelist
  res <- merge(genes8, df.name, by.x = 'V1', by.y = 'gene')
  
  # reshaping data
  res <- res %>%
    spread(key = 'V1', value = 'FPKM') %>%
    gather(key = 'gene', value = 'FPKM',-c(sample, MYCN_Status, MYCN, RISK))
  
  return(res)
}

# function to plot data
plotCorr <- function(df,colorby) {
  # plot correlation
  ggplot(df, aes(x = MYCN, y = FPKM, color = eval(parse(text = colorby)))) +
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    xlab("MYCN expression") +
    ylab("Gene expression") +
    theme_bw() +
    #stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
     stat_cor(method = "pearson", aes(label = paste(after_stat(r.label),after_stat(rr.label),after_stat(p.label), sep = "~` `~` `~"))) +
    facet_wrap(. ~ gene, scales = "free") +
    theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
  
}
