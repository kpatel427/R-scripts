# obtain expression values for the given list of genes from R2 genomics
# using soft file to get all annotation data esp - EFS/OS survival data

library(data.table)

data <-  read.delim("GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt", header = T)
data <- data %>%
  gather(variable,value, -c("RefSeqID"))

soft.file <- soft.file <- read.delim("GSE62564_family.soft", header = F)
# Sample_ID
sample.ID <-  as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_title" ), ])
colnames(sample.ID) <- "V1"
setDT(sample.ID)[, c("V1", "sample_title") := tstrsplit(sample.ID$V1, "=", fixed=TRUE)]
setDT(sample.ID)[, c("sample_ID") := tstrsplit(sample.ID$sample_title, "[2]", fixed=TRUE)]
sample.ID <- subset(sample.ID,select = c("sample_ID") )

# Sample_titles
sample.title <-  as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_geo_accession" ), ])
colnames(sample.title) <- "V1"
setDT(sample.title)[, c("V1", "sample_title") := tstrsplit(sample.title$V1, "=", fixed=TRUE)]


#efs_day 
efs.day <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = efs day" ), ])
colnames(efs.day) <- "V1"
setDT(efs.day)[, c("V1", "sample_title") := tstrsplit(efs.day$V1, "=", fixed=TRUE)]
setDT(efs.day)[, c("sample_title","EFS day") := tstrsplit(efs.day$sample_title, ":", fixed=TRUE)]
efs.day <- subset(efs.day, select = c("EFS day"))



#efs_bin
efs.bin <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = efs bin" ), ])
colnames(efs.bin) <- "V1"
setDT(efs.bin)[, c("V1", "sample_title") := tstrsplit(efs.bin$V1, "=", fixed=TRUE)]
setDT(efs.bin)[, c("sample_title","EFS bin") := tstrsplit(efs.bin$sample_title, ":", fixed=TRUE)]
efs.bin <- subset(efs.bin, select = c("EFS bin"))

#os_day
os.day <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = os day" ), ])
colnames(os.day) <- "V1"
setDT(os.day)[, c("V1", "sample_title") := tstrsplit(os.day$V1, "=", fixed=TRUE)]
setDT(os.day)[, c("sample_title","OS day") := tstrsplit(os.day$sample_title, ":", fixed=TRUE)]
os.day <- subset(os.day, select = c("OS day"))


#os_bin
os.bin <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = os bin" ), ])
colnames(os.bin) <- "V1"
setDT(os.bin)[, c("V1", "sample_title") := tstrsplit(os.bin$V1, "=", fixed=TRUE)]
setDT(os.bin)[, c("sample_title","OS bin") := tstrsplit(os.bin$sample_title, ":", fixed=TRUE)]
os.bin <- subset(os.bin, select = c("OS bin"))


# high_risk
high.risk <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = high risk:" ), ])
colnames(high.risk) <- "V1"
setDT(high.risk)[, c("V1", "sample_title") := tstrsplit(high.risk$V1, "=", fixed=TRUE)]
setDT(high.risk)[, c("sample_title","High risk") := tstrsplit(high.risk$sample_title, ":", fixed=TRUE)]
high.risk <- subset(high.risk, select = c("High risk"))

# Sex
gender <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = Sex:" ), ])
colnames(gender) <- "V1"
setDT(gender)[, c("V1", "sample_title") := tstrsplit(gender$V1, "=", fixed=TRUE)]
setDT(gender)[, c("sample_title","Gender") := tstrsplit(gender$sample_title, ":", fixed=TRUE)]
gender <- subset(gender, select = c("Gender"))

# training/validation:
training_validation <- as.data.frame(soft.file[ (soft.file$V1 %like% "!Sample_characteristics_ch1 = training/validation:" ), ])
colnames(training_validation) <- "V1"
setDT(training_validation)[, c("V1", "sample_title") := tstrsplit(training_validation$V1, "=", fixed=TRUE)]
setDT(training_validation)[, c("sample_title","training/validation") := tstrsplit(training_validation$sample_title, ":", fixed=TRUE)]
training_validation <- subset(training_validation, select = c("training/validation"))



survival.data <-  as.data.frame(cbind(sample.ID$sample_ID,sample.title$sample_title, efs.day$`EFS day`, efs.bin$`EFS bin`,
                                      os.day$`OS day`, os.bin$`OS bin`, high.risk$`High risk`, gender$Gender, training_validation$`training/validation`))
colnames(survival.data) <-  c("Sample ID","Sample_title","EFS day", "EFS bin","OS day","OS bin", "High Risk","Gender","training/validation")
survival.data$Sample_title <-  tolower(survival.data$Sample_title)
survival.data$Sample_title <- gsub(" ","",survival.data$Sample_title, fixed = T)
survival.data$`Sample ID` <-  as.character(survival.data$`Sample ID`)
survival.data$`Sample ID` <- gsub(" ","",survival.data$`Sample ID`, fixed = T)

survival.data <- merge(survival.data,data, by.x = "Sample ID", by.y = "variable")


# merging gse62564 (downloaded from R2 genomics) with survival data (soft file) by RefSeq ID
gse62564 <-  read.delim("GSE62564.txt", header = T)
gse62564 <-  gse62564 %>%
  gather(variable,value, -c("X.H.hugo","probeset"))

gse62564 <-  subset(gse62564, select = c("X.H.hugo","probeset"))
gse62564 <- gse62564[unique(gse62564$X.H.hugo),]

gse62564 <-  merge(gse62564, survival.data, by.x = "probeset", by.y = "RefSeqID")

colnames(gse62564) <-  c("probeset","Gene Symbol","Sample ID", "Sample_title","EFS day","EFS bin","OS day","OS bin",
                         "High Risk", "Gender","training/validation","log2RPM")


# saving results
save(gse62564, file = paste0(Sys.Date(), "-R2_SEQC_dataset.rda") )
write.table(gse62564, file = paste0(Sys.Date(),"-R2_SEQC_dataset.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
