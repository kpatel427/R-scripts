library(ggplot2)
library(ggforce)
data <- read.delim("FPKM_3_transcripts_per_exon.txt")

sub_dat <- subset(data, data$Sample == "COG-N-100x-D", c(colnames(data)))
tx_dat <- subset(sub_dat, sub_dat$transcript == "ENST000000000.1", c(colnames(sub_dat)))


a <- ggplot(sub_dat, aes(exon_n, FPKM)) +
  geom_bar(aes(fill = transcript),stat = "identity") +
  scale_x_continuous(breaks = pretty(sub_dat$exon_n, n = 35)) +
  labs(title = "FPKM for exons in model COG-N-100x") +
  facet_wrap(. ~ sub_dat$transcript, nrow = 3) +
  theme_bw()

ggsave(paste0(Sys.Date(),"-exon-number-all-TXs-COG-N-100x-.pdf"), a, width=45, height=20, units = "cm")


# Grepping pattern to get all rows matching the string
nbl <- data[grep("COG-", data$Model), ]

unique(nbl$Model) # 26 Models

ggplot(nbl, aes(exon_n, FPKM)) +
  geom_bar(aes(fill = transcript),stat = "identity") +
  scale_x_continuous(breaks = pretty(nbl$exon_n, n = 35)) +
  labs(title = "FPKM for exons in all NBL models") +
  facet_wrap(. ~ transcript, nrow = 3) +
  theme_bw()



# total number of Models = 26
n_pages_needed = 13
pdf("NBL-Tx-ATRX.pdf", width=15, height=9)

for (i in seq_len(n_pages_needed)) {
  print(ggplot(nbl,aes(nbl$exon_n, nbl$FPKM, group = transcript, color = transcript)) +
          geom_bar(aes(fill = transcript), stat = "identity", width = 0.5) +
          scale_x_continuous(breaks = pretty(nbl$exon_n, n = 36)) +
          theme_bw()+
          labs(title = "FPKM for exons in all NBL models", x = "exon number", y = "FPKM") + 
          theme(legend.text=element_text(size=10)) +
          facet_grid_paginate(nbl$transcript ~ nbl$Model, scales = "fixed", ncol = 2 , nrow = 3, page = i) + 
          theme(strip.text.y = element_blank()) )

}
dev.off()
