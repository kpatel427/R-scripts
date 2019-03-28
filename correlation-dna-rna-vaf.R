# correlation between rna.vaf and dna.vaf
library(ggpubr)
library(dplyr)
library(plyr)

df1 %>%
  group_by(df1$TSB) %>%
  ggscatter(., x = "VAF", y = "rna.vaf",
          shape = 1, alpha = 0.6,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "DNA.VAF", ylab = "RNA.VAF") +
          theme_bw()



# Use R2 instead of R
df1 %>%
  group_by(df1$TSB) %>%
  ggscatter(., x = "VAF", y = "rna.vaf",shape = 1, alpha = 0.6,
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "DNA.VAF", ylab = "RNA.VAF") +
            theme_bw()+
            stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                label.x = 0.15, label.y = 1) +
            # Save as PDF
            ggsave(paste0(Sys.Date(),"-correlation-cosmic-DNA-RNA-VAF.pdf"), width=25, height=12, device = "pdf", units = "cm")


