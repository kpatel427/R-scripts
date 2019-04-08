library(dplyr)

# read matrix to get FPKM expression values
matrix <- read.delim("GSE49711_SEQC_NB_TUC_G_log2.txt", header = T)
matrix <- matrix %>%
  gather("Sample_title",value,-c("X00gene_id"))
  
# convert from log2(FPKM + 1) to FPKM
matrix$value <- ( ( (matrix$value) ^2) - 1)

# splitting sample_titles
setDT(matrix)[, c("V1","V2","V3") := tstrsplit(matrix$Sample_title, "_", fixed=TRUE)]
setDT(matrix)[,"new_sample_title" := paste0(matrix$V1,"_",matrix$V2)]

matrix <-  subset(matrix, select = c("X00gene_id","new_sample_title","value"))
colnames(matrix)[colnames(matrix)=="new_sample_title"] <- "Sample_title"


# reading in soft file to get risk and mycn status for sample_titles
family.soft <-  read.delim("GSE49711_family.soft", header = F)
# Sample_titles
sample.title <-  as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_title" ), ])
colnames(sample.title) <- "V1"
setDT(sample.title)[, c("V1", "sample_title") := tstrsplit(sample.title$V1, "=", fixed=TRUE)]

# high_risk
risk <-  as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_characteristics_ch1 = high risk"), ])
colnames(risk) <- "V1"
setDT(risk)[, c("V1", "high_risk") := tstrsplit(risk$V1, ":", fixed=TRUE)]

# mycn_status
mycn_status <- as.data.frame(family.soft[ (family.soft$V1 %like% "!Sample_characteristics_ch1 = mycn status"), ])
colnames(mycn_status) <- "V1"
setDT(mycn_status)[, c("V1", "mycn_status") := tstrsplit(mycn_status$V1, ":", fixed=TRUE)]



# combining sample with mycn status and risk
sample.mycn.risk <- as.data.frame(cbind(sample.title$sample_title, mycn_status$mycn_status, risk$high_risk))
colnames(sample.mycn.risk) <-  c("Sample_title","mycn_status","high_risk")
# Sample titles were of type integer, changing to character
sample.mycn.risk$Sample_title <- as.character(sample.mycn.risk$Sample_title)
# removing blank space before the sample_title
sample.mycn.risk$Sample_title <- gsub(" ","",sample.mycn.risk$Sample_title, fixed = T)



# GOI gene list
gene <- scan("identified_genelist.txt", what = "", sep = ",")
df.gene <- as.data.frame(gene)

# to get only those genes from the matrix, present in the list of identified genes
matrix.genelist <- merge(matrix,df.gene, by.x = "X00gene_id", by.y = "gene")
#merging to get risk and mycn status information
matrix.genelist <-  merge(matrix.genelist,sample.mycn.risk, by = "Sample_title")

# filter out rows with FPKM < 1.0
matrix.genelist <- matrix.genelist %>%
  filter(matrix.genelist$value > 1.00)

matrix.genelist <-  plot.matrix.genelist

# Renaming factor levels
levels(matrix.genelist$high_risk) <- c("Low Risk","High Risk")
levels(matrix.genelist$mycn_status) <-  c("Non Amplified", "Amplified", "N/A")


# boxplots showing FPKM expression at different risk levels for all genes in gene list
ggplot(matrix.genelist, aes(x = X00gene_id, y = value, fill = high_risk)) +
  geom_boxplot(aes(fill = factor(high_risk))) +
  labs(title = "FPKM expression at different risk levels", x = "Genes", y = "FPKM") +
  guides(fill=guide_legend(title= "Risk")) +
  theme_bw()

ggsave("boxplot-FPKM-risk-levels.pdf", width = 15, height = 9)
  
# boxplots showing FPKM expression for mycn_status for all genes in gene list
ggplot(matrix.genelist, aes(x = X00gene_id, y = value, fill = mycn_status)) +
  geom_boxplot(aes(fill = factor(mycn_status))) +
  labs(title = "FPKM expression for mycn_status", x = "Genes", y = "FPKM") +
  guides(fill=guide_legend(title= "MYCN status")) +
  theme_bw()
ggsave("boxplot-FPKM-mycn-status.pdf", width = 15, height = 9)


# correlation of MYCN expression with the rest of the genes
matrix.genelist.subset <- subset(matrix.genelist, matrix.genelist$X00gene_id != "MYCN", select = c(colnames(matrix.genelist)) )
mycn.expr <- subset(matrix.genelist, matrix.genelist$X00gene_id == "MYCN", select = c("Sample_title","value") )
colnames(mycn.expr)[colnames(mycn.expr)=="value"] <- "mycn_expression_value"

matrix.genelist.subset <- merge(matrix.genelist.subset, mycn.expr, all.x = T ,by = "Sample_title")

# Renaming factor levels
levels(matrix.genelist.subset$high_risk) <- c("Low Risk","High Risk")
levels(matrix.genelist.subset$mycn_status) <-  c("Non Amplified", "Amplified", "N/A")


# plotting correlation ~ mycn_status
ggscatter(matrix.genelist.subset, x = "value", y = "mycn_expression_value",
          shape = 1, alpha = 0.6,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Genes", ylab = "Mycn expression") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 0.15, label.y = 1) +
  facet_wrap(~mycn_status, scales="free_x") +
  labs(x = "Genes of interest", y = "MYCN expression") + 
  theme_bw()
ggsave("correlation-mycn-status.pdf", width = 10, height = 6)


# plotting correlation ~ Risk
ggscatter(matrix.genelist.subset, x = "value", y = "mycn_expression_value",
          shape = 1, alpha = 0.6,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Genes", ylab = "Mycn expression") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 0.15, label.y = 1) +
  facet_wrap(~high_risk, scales="free_x") +
  labs(x = "Genes of interest", y = "MYCN expression") + 
  theme_bw()
ggsave("correlation-risk-levels.pdf", width = 10, height = 6)


# ggplot(matrix.genelist.subset) +
#   geom_jitter(aes(value, mycn_expression_value, colour=mycn_status),) + 
#   geom_smooth(aes(value, mycn_expression_value, colour=mycn_status), method=lm, se=FALSE) +
#   facet_wrap(~mycn_status, scales="free_x") +
#   labs(x = "Genes of interest", y = "MYCN expression")


#------------------------------------------ Wilcox's test -----------------------------------------------------------#
# Comparing the means of two independent groups (high risk and low risk)
# by Risk := comparing high risk (amp/non-amp) vs low risk (amp/non-amp)

#----------------------- Risk ------------------------------#
# 1. High risk amplified : Low risk amplified
# there is no low risk amplified expression for the set of genes provided

# 2. High risk non-amplified : Low risk non-amplified

risk.nonamp <- subset(matrix.genelist, matrix.genelist$mycn_status == " 0", select = c(colnames(matrix.genelist)) )

res1 <-  wilcox.test(value ~ high_risk, data = risk.nonamp,
            exact = FALSE)
res1$p.value

# p-value = 0.4607 which is not less than 0.05, 
# so there is no significant difference in median expression levels of genes between High risk non-amplified & Low risk non-amplified

# Renaming factor levels
levels(risk.nonamp$high_risk) <- c("Low Risk","High Risk")

# plot
ggboxplot(risk.nonamp, x = "high_risk", y = "value",
          color = "high_risk", palette = "jco") + 
          stat_compare_means()+ 
          stat_compare_means(method = "wilcox.test") +
          labs(title = "For non-amplified samples", x = "Risk", y = "FPKM") +
          theme(legend.title = element_blank(), legend.position = "none")

ggsave("wilcox-highrisk.pdf", width = 7, height = 5)


#----------------------- Amplification status ------------------------------#
# Comparing the means of two independent groups (Amp and non-amp)
# by Amplification status := comparing Amp (high risk/low risk) vs non-amp (high risk/low risk)

# 3. High risk amplified : High risk non-amplified
amp.highrisk <- subset(matrix.genelist, matrix.genelist$high_risk == " 1" & matrix.genelist$mycn_status != " N/A", select = c(colnames(matrix.genelist)) )

res2 <-  wilcox.test(value ~ mycn_status, data = amp.highrisk,
                     exact = FALSE)
res2$p.value
# p-value = 0.000000000004442 which is significantly less than 0.05, 
# so there is significant difference in median expression levels of genes between high rosk amplified and non-amplified patients

# Renaming factor levels
levels(amp.highrisk$mycn_status) <- c("Non-amplified","Amplified","N/A")

# plot
ggboxplot(amp.highrisk, x = "mycn_status", y = "value",
          color = "mycn_status", palette = "jco") + 
          stat_compare_means()+ 
          stat_compare_means(method = "wilcox.test") +
          labs(title = "For High-risk samples", x = "MYCN status", y = "FPKM")+
          theme(legend.title = element_blank(), legend.position = "none")

ggsave("wilcox-mycn-status.pdf", width = 7, height = 5)


# 4. low risk amplified : low risk non-amplified
# there is no low risk amplified expression for the set of genes provided


#--------------------------------------------------------------------------#
# wilcox's test pairwise for genes (high risk-low risk)

risk.wilcox <-  matrix.genelist %>%
  group_by(X00gene_id) %>% 
  summarise(p=wilcox.test(value~high_risk)$p.value)


# wilcox's test pairwise for genes (mycn status)
mycn.wilcox <-  matrix.genelist %>%
  filter(mycn_status != "N/A") %>%
  group_by(X00gene_id) %>% 
  summarise(p=wilcox.test(value~mycn_status)$p.value)


