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



# Kristen's gene list
gene <- scan("identified_genelist.txt", what = "", sep = ",")
df.gene <- as.data.frame(gene)

# to get only those genes from the matrix, present in the list of identified genes
matrix.genelist <- merge(matrix,df.gene, by.x = "X00gene_id", by.y = "gene")
#merging to get risk and mycn status information
matrix.genelist <-  merge(matrix.genelist,sample.mycn.risk, by = "Sample_title")

# filter out rows with FPKM < 1.0
matrix.genelist <- matrix.genelist %>%
  filter(matrix.genelist$value > 1.00)

plot.matrix.genelist <- matrix.genelist

# Renaming factor levels
levels(matrix.genelist$high_risk) <- c("Low Risk","High Risk")
levels(matrix.genelist$mycn_status) <-  c("Non Amplified", "Amplified", "N/A")


# boxplots showing FPKM expression at different risk levels for all genes in gene list
ggplot(matrix.genelist, aes(x = X00gene_id, y = value, fill = high_risk)) +
  geom_boxplot(aes(fill = factor(high_risk))) +
  labs(title = "FPKM expression at different risk levels", x = "Genes", y = "FPKM") +
  guides(fill=guide_legend(title= "Risk")) +
  theme_bw()
  
# boxplots showing FPKM expression for mycn_status for all genes in gene list
ggplot(matrix.genelist, aes(x = X00gene_id, y = value, fill = mycn_status)) +
  geom_boxplot(aes(fill = factor(mycn_status))) +
  labs(title = "FPKM expression for mycn_status", x = "Genes", y = "FPKM") +
  guides(fill=guide_legend(title= "MYCN status")) +
  theme_bw()




