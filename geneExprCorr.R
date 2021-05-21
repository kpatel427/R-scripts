# script to correlate MYCN & c-MYC with BPTF across cell lines, patients and PDXs
# also check correlation with CRC (Core Regulatory Circuit) genes 
# setwd("~/KP/BPTF")

library(tidyverse)
library(ggpubr)

# READ DATA ---------------
# ...read data: Cell Line data -----------

# ...read data: Patient data -----------
# TARGET

# GMKF

# ...read data: PDX data -----------

# to keep only NBL PDXs
PPTC_FPKM_hg19_Wheeler_subtracted_mData <- PPTC_FPKM_hg19_Wheeler_subtracted_mData[PPTC_FPKM_hg19_Wheeler_subtracted_mData$CANCER_TYPE_DETAILED == "Neuroblastoma",]


# READ Gene list ---------------
CRC = c("MYCN","PHOX2B","TBX2","ASCL1","HAND2","LMO1","GATA3","ISL1", "MYC", "BPTF","GATA2")


# WRANGLE & SUBSET DATA ---------------
# subset metadata
TARGET_NBL_FPKM_PST_mData <- TARGET_NBL_FPKM_PST_mData[,c(2,14)]
names(TARGET_NBL_FPKM_PST_mData)[2] <- 'MYCN_Status'
PPTC_FPKM_hg19_Wheeler_subtracted_mData <- PPTC_FPKM_hg19_Wheeler_subtracted_mData[,c(1,6)]
gm <- gm %>%
  column_to_rownames(var = 'gene')
gabby.clinical <- gabby.clinical %>%
  rownames_to_column(var = 'Sample') %>%
  select(Sample, mycn_status)

# function definition
wrangleData <- function(df, df.meta, keyname){
  new <- df %>%
    rownames_to_column(var = 'Gene') %>%
    filter(Gene %in% CRC) %>%
    gather(key = !!keyname, value = 'FPKM', -Gene)
  
  # merge with metadata
  new <- merge(new, df.meta, by = eval(keyname))
  
  # reshape data
  new <- new %>%
    spread(key = 'Gene', value = 'FPKM')
  
  return(new)
}

# function call
cl <- wrangleData(df = STAR_FPKM_40cells_genes,
                  df.meta = cellline_mData,
                  keyname = 'CellLine')

target <- wrangleData(df = TARGET_NBL_FPKM_PST_data,
                      df.meta = TARGET_NBL_FPKM_PST_mData,
                      keyname = 'TARGET.ID')

gmkf <- wrangleData(df = gm,
                    df.meta = gabby.clinical,
                    keyname = 'Sample')

pdx <- wrangleData(df = PPTC_FPKM_hg19_Wheeler_subtracted_data,
                   df.meta = PPTC_FPKM_hg19_Wheeler_subtracted_mData,
                   keyname = 'MODEL')


# reshape data frames
cl <- cl %>%
  gather(key = 'Genes', value = 'FPKM', -c(1:4,6))

target <- target %>%
  gather(key = 'Genes', value = 'FPKM', -c(1,2,4))


gmkf <- gmkf %>%
  gather(key = 'Genes', value = 'FPKM', -c(1,2,4))

pdx <- pdx %>%
  gather(key = 'Genes', value = 'FPKM', -c(1,2,4))


# Correcting MYCN status column ------------
target$MYCN_Status <- gsub('single_copy', 'Non-amplified', target$MYCN_Status)
target$MYCN_Status <- gsub('amplified', 'Amplified', target$MYCN_Status)

gmkf$mycn_status <- gsub('not amplified', 'Non-amplified', gmkf$mycn_status)
gmkf$mycn_status <- gsub('amplified', 'Amplified', gmkf$mycn_status)
gmkf$mycn_status <- gsub('unknown', 'Unknown', gmkf$mycn_status)

pdx$MYCN_Status <- gsub('amp', 'Amplified', pdx$MYCN_Status)
pdx$MYCN_Status <- gsub('non-amp', 'Non-amplified', pdx$MYCN_Status)

# plot correlations ---------
# ...cell lines ----------
ggplot(cl, aes(reorder(CellLine,BPTF), BPTF, fill = MYCN_Status)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# 1. MYCN ~ BPTF correlation
# fixing order of facets
cl$Genes = factor(cl$Genes, levels=c("MYCN","MYC","ASCL1","GATA3","GATA2","HAND2","ISL1" ,"LMO1", "PHOX2B", "TBX2"))
p1 <- ggplot(cl, aes(FPKM, BPTF, color = MYCN_Status)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
  labs(x = '', y = 'BPTF expression') +
  facet_wrap(~ Genes, ncol = 3, scales = 'free_x') +
  theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
ggsave(p1, filename = paste0(Sys.Date(),'_correlation_cl_CRC_bptf_mycnStatus.pdf'), width = 10, height = 10)


#' #' 2. cMYC ~ BPTF correlation
#' p2ggplot(cl, aes(MYC, BPTF, color = MYCN_Status)) + 
#'   geom_point() +
#'   theme_bw() +
#'   geom_smooth(method=lm, se=FALSE) +
#'   stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
#'   labs(x = '', y = 'BPTF expression') +
#'   facet_wrap(~ Genes, ncol = 2, scales = 'free_x') +
#'   theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
#' ggsave(p2, filename = paste0(Sys.Date(),'_correlation_cl_cmyc_bptf_mycnStatus.pdf'), width = 10, height = 7)


# ...Patients = TARGET --------------
# 1. MYCN ~ BPTF correlation
target$Genes = factor(target$Genes, levels=c("MYCN","MYC","ASCL1","GATA3","GATA2","HAND2","ISL1" ,"LMO1", "PHOX2B", "TBX2"))
p3 <- ggplot(target, aes(FPKM, BPTF, color = MYCN_Status)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
  labs(x = '', y = 'BPTF expression') +
  facet_wrap(~ Genes, ncol = 3, scales = 'free_x') +
  theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
ggsave(p3, filename = paste0(Sys.Date(),'_correlation_target_mycn_bptf_mycnStatus.pdf'), width = 10, height = 10)


#' 2. cMYC ~ BPTF correlation
# p4 <- ggplot(target, aes(MYC, BPTF, color = MYCN_Status)) + 
#   geom_point() +
#   theme_bw() +
#   xlim(0,300) +
#   geom_smooth(method=lm, se=FALSE) +
#   stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
#   labs(x = 'MYC expression', y = 'BPTF expression') +
#   theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
# ggsave(p4, filename = paste0(Sys.Date(),'_correlation_target_cmyc_bptf_mycnStatus.pdf'), width = 10, height = 7)


# ...Patients = GMKF ----------------
# 1. MYCN ~ BPTF correlation
gmkf$Genes = factor(gmkf$Genes, levels=c("MYCN","MYC","ASCL1","GATA3","GATA2","HAND2","ISL1" ,"LMO1", "PHOX2B", "TBX2"))
p5 <- ggplot(gmkf, aes(FPKM, BPTF, color = mycn_status)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
  labs(x = '', y = 'BPTF expression') +
  facet_wrap(~ Genes, ncol = 3, scales = 'free_x') +
  theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
ggsave(p5, filename = paste0(Sys.Date(),'_correlation_GMKF_mycn_bptf_mycnStatus.pdf'), width = 10, height = 10)



#' 2. cMYC ~ BPTF correlation
# p6 <- ggplot(gmkf, aes(MYC, BPTF, color = mycn_status)) + 
#   geom_point() +
#   theme_bw() +
#   xlim(0,150) +
#   geom_smooth(method=lm, se=FALSE) +
#   stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
#   labs(x = 'MYC expression', y = 'BPTF expression') +
#   theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
# ggsave(p6, filename = paste0(Sys.Date(),'_correlation_GMKF_cmyc_bptf_mycnStatus.pdf'), width = 10, height = 7)


# ...PDX ----------------
# COG-N-589x is non-amp
# 1. MYCN ~ BPTF correlation
pdx$Genes = factor(pdx$Genes, levels=c("MYCN","MYC","ASCL1","GATA3","GATA2","HAND2","ISL1" ,"LMO1", "PHOX2B", "TBX2"))
p7 <- ggplot(pdx, aes(FPKM, BPTF, color = MYCN_Status)) + 
  geom_point() +
  theme_bw() +
  geom_smooth(method=lm, se=FALSE) +
  stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
  labs(x = '', y = 'BPTF expression') +
  facet_wrap(~ Genes, ncol = 3, scales = 'free_x') +
  theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
ggsave(p7, filename = paste0(Sys.Date(),'_correlation_PDX_mycn_bptf_mycnStatus.pdf'), width = 10, height = 7)


#' 2. cMYC ~ BPTF correlation
# p8 <- ggplot(pdx, aes(MYC, BPTF, color = MYCN_Status)) + 
#   geom_point() +
#   theme_bw() +
#   geom_smooth(method=lm, se=FALSE) +
#   stat_cor(method = "pearson", aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~` `~` `~"))) +
#   labs(x = 'MYC expression', y = 'BPTF expression') +
#   theme(strip.text = element_text(size=8), legend.position = "bottom", legend.title = element_text(size=5))
# ggsave(p8, filename = paste0(Sys.Date(),'_correlation_PDX_cmyc_bptf_mycnStatus.pdf'), width = 10, height = 7)
# 











