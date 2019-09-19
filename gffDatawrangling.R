library(data.table)
library(stringr)

# reading the file
gff <- read.delim('test.txt', header = F)

# adding another column with exon_number NA where exon_number is missing
gff$V15 <- ifelse(grepl('exon_number',gff$V9),'','exon_number NA;')

# splitting before coverage because we want exon_number to be present before coverage
gff[,c('V10','V11')] <- tstrsplit(gff$V9,'; cov')

# concatenating 'exon_number NA before coverage'
gff$V9 <- ifelse(!is.na(gff$V15), paste0(gff$V10,' ; ',gff$V15,'; cov', gff$V11), '')

# subsetting
gff <- gff[,c(1:9)]

# splitting again
gff[,c('V10','V11','V12','V13','V14','V15', 'V16')] <- tstrsplit(gff$V9, ';', fixed = TRUE)
