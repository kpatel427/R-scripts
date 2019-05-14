library(tidyr)
library(tidyverse)
library(dplyr)
library(grid)

# color mapping
# to color cell lines with specific colors
fillPalette <-  c("#104A7F","#C9C9C9","#97D1A9","#FDFDAC","#FDD6A2","#D49DC7")
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


ggplot(hprt1, aes(x=X,y=value, fill =X, colour = mycn_status)) +
  geom_bar(stat="identity") +
  facet_wrap( ~ variable, scale = "free", ncol = 4) +
  scale_fill_brewer(palette="Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Normalized to HPRT1", x="Cell lines",y="expression values") +
  guides(fill = guide_legend("Cell Line"), colour = guide_legend("MYCN Status")) +
  theme_bw() +
  theme(axis.text = element_text(size=10, angle = 90),
        axis.title = element_text(size = 15))

ggsave(filename = paste0(Sys.Date(),"-RTqPCR-HPRT1-barplot.pdf"),plot = last_plot(),width = 10, height=15, device = "pdf")



final.gene <- unique(gapdh$variable)
write.table(final.gene, file = "2019-05-10-final-gene-list.txt", sep = "\t", row.names = F, col.names = F, quote = F)
