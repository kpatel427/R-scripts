library(tidyr)
library(tidyverse)
library(dplyr)
library(grid)

# color mapping
# to color cell lines with specific colors
fillPalette <-  c("#104A7F","#C9C9C9","#97D1A9","#FDFDAC","#FDD6A2","#DBCDE9")
cl <- c("SK-N-FI","SK-N-AS","NB-69","NGP","Kelly","SK-N-BE(2)")
colorPalette <-  c("#47A1A7","#E97D72")
mycn <- c("Non-amplified","Amplified")

names(fillPalette) <- cl
names(colorPalette) <- mycn



gapdh <- read.delim("RT_qPCR_GAPDH.txt", sep = "\t", header = T)
gapdh <-  gapdh %>%
  gather(.,variable,value,-c("X")) %>%
  mutate(mycn_status = ifelse(X %in% c("Kelly","NGP","SK-N-BE(2)"),"Amplified","Non-amplified"))

gapdh <-  subset(gapdh, gapdh$X != 'NBL-S', select = c(colnames(gapdh)))



ggplot(gapdh, aes(x=X,y=value, fill =X, colour = mycn_status)) +
  geom_bar(stat="identity") +
  facet_wrap( ~ variable, scales = "free", ncol = 4) +
  scale_fill_manual(values=fillPalette) +
  scale_color_manual(values = colorPalette) +
  labs(title = "Normalized to GAPDH", x="Cell lines",y="expression values") +
  guides(fill = guide_legend("Cell Line"), colour = guide_legend("MYCN Status")) +
  theme_bw() +
  theme(axis.text = element_text(size=10, angle = 90),
        axis.title = element_text(size = 15))


ggsave(filename = paste0(Sys.Date(),"-RTqPCR-GAPDH-barplot.pdf"),plot = last_plot(),width = 10, height=15, device = "pdf")

  
hprt1 <- read.delim("RT_qPCR_HPRT1.txt", sep = "\t", header = T)
hprt1 <-  hprt1 %>%
  gather(.,variable,value,-c("X")) %>%
  mutate(mycn_status = ifelse(X %in% c("Kelly","NGP","SK-N-BE(2)"),"Amplified","Non-amplified"))

hprt1 <-  subset(hprt1, hprt1$X != 'NBL-S', select = c(colnames(hprt1)))


ggplot(hprt1, aes(x=X,y=value, fill =X, colour = mycn_status)) +
  geom_bar(stat="identity") +
  facet_wrap( ~ variable, scale = "free", ncol = 4) +
  scale_fill_manual(values=fillPalette) +
  scale_color_manual(values = colorPalette) +
  #scale_fill_brewer(palette="Pastel2") +
  #scale_color_brewer(palette = "Dark2") +
  labs(title = "Normalized to HPRT1", x="Cell lines",y="expression values") +
  guides(fill = guide_legend("Cell Line"), colour = guide_legend("MYCN Status")) +
  theme_bw() +
  theme(axis.text = element_text(size=10, angle = 90),
        axis.title = element_text(size = 15))

ggsave(filename = paste0(Sys.Date(),"-RTqPCR-HPRT1-barplot.pdf"),plot = last_plot(),width = 10, height=15, device = "pdf")



final.gene <- unique(gapdh$variable)
write.table(final.gene, file = "2019-05-10-final-gene-list.txt", sep = "\t", row.names = F, col.names = F, quote = F)


#-------------------- Patient Distribution for GSE49711 --------------------------#

pd <- read.delim("../patient_dist_GSE49711.txt", sep = "\t", header = TRUE)
colnames(pd) <-  c("Stage","MYCN_NonAmp","MYCN_Amp","ND","Total","Total_%")

pd <-  pd %>%
  select(colnames(pd)[1:4]) %>%
  gather(key, value, -c(Stage))


pd$key <-  gsub("MYCN_NonAmp","Non-Amplified", pd$key)
pd$key <-  gsub("MYCN_Amp","Amplified", pd$key)
pd$key <-  gsub("ND","N/A", pd$key)
  

ggplot(pd, aes(Stage,value, fill = key)) +
  geom_bar(stat="identity", position = position_dodge()) +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() + 
  scale_fill_manual(values = c("#00BA38","#619CFF","#F8766D")) +
  labs(title = "Patient Distribution GSE49711", x = "Stage", y = "Number of Patients") +
  guides(fill = guide_legend("MYCN Status")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size=10))

ggsave(filename = paste0(Sys.Date(),"-patient-distribution-GSE49711-barplot.pdf"),plot = last_plot(),width = 7, height=5, device = "pdf")
