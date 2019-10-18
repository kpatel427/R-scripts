# script to analyze MAST results: mast.txt
# get gene and motif information into seperate files
setwd('~/KP/H3K27me3/motif-analysis/46_candidate genes/MAST/knownMotifs_MAST')

library(data.table)
library(dplyr)
library(tidyverse)

# 1. Data ---------- 
# _1.1 reading hg19 tss data ---------
tss_hg19 <-  read.delim('~/KP/hg19/refGene_hg19_TSS.bed.txt', header = F)
# change column type
tss_hg19$V1 <-  as.character(tss_hg19$V1)

# _1.2 reading refseq IDS ~ 46 Genes -----------
refseqID <-  read.delim('~/KP/H3K27me3/motif-analysis/46_candidate genes/refseqIDS_46Genes.txt', header = F, sep = '.')
refseqID <-  as.data.frame(refseqID[,c(1)])
names(refseqID)[1] <-  'ID'

# _1.3 reading 46 candidate genes --------------
genes46 <- read.delim('~/KP/H3K27me3/motif-analysis/46_candidate genes/both-immune-inter-HR.txt', header = F)

# _1.4 reading hg19 refGene data --------------
refGene_hg19 <- read.delim('~/KP/hg19/hg19.refGene.txt', header = F) 
refGene_hg19 <-  refGene_hg19[,c(2,3,4,7,8,13)]
colnames(refGene_hg19) <-  c('ID','Chr','strand','start','stop','gene')


# 2.1 TSS + 46 genes -------------
getTSS <- merge(genes46, tss_hg19, by.x = 'V1', by.y = 'V5', all.x = T)
getTSS$V4 <-  NULL
getTSS <-  getTSS[!duplicated(getTSS),]
colnames(getTSS) <-  c('gene','Chr','TSS_start','TSS_stop','strand')


# 2.2 refGene + refseq IDS -----------------
getrefGene <-  merge(refseqID, refGene_hg19, by = 'ID', all.x = T)
getrefGene <-  na.omit(getrefGene)
names(getrefGene)[1] <-  'refGene_ID'

# 3. TSS + refGene -------------
refGene_TSS <- merge(getTSS, getrefGene, by = c('gene','Chr','strand'))

# _3.1 calculating distance of TSS from gene ----------------
refGene_TSS <- refGene_TSS %>%
  mutate('Distance.to.TSS' = ifelse(strand == '+', TSS_start - start, stop - TSS_stop))




# 4. reading motif diagram files ------------------
motif_diagram <- read.delim('motif_diagrams.txt', header = F, sep = ';')
colnames(motif_diagram) <-  c('seq_name','Evalue','motif_diagram')
motif_diagram$motif_diagram <-  as.character(motif_diagram$motif_diagram)

# _4.1 reading motif_info ---------------------
motif_info <-  read.delim('motifInfo.txt', header = T, sep = ':')
motif_info$MOTIF <-  as.character(motif_info$MOTIF)

# handing outlier cases
l1 <- gsub("\\]\\-\\[","\\]\\-0\\-\\[",motif_diagram$motif_diagram, perl = TRUE)
l1 <- gsub("\\-\\<\\-|\\-\\<\\+","\\-[\\-",l1, perl = TRUE)
l1 <- gsub("\\>\\-","]\\-",l1, perl = TRUE)

# trying to get string positions in one column and motif locations in the next column
string_pos <- strsplit(l1,'\\-\\[[+-][0-9]+\\]\\-', perl = T) # splitting at motif locations
motif_location <- strsplit(l1,'[0-9]+\\-|\\-[0-9]+\\-', perl = T) # splitting at string locations

for(i in 1:length(motif_location)){
  #print(l1[i])
  motif_location[[i]] <- motif_location[[i]][-1]
  motif_location[[i]] <- append(motif_location[[i]],"")

}

# QC
for(i in 1:length(motif_location)){
  if( length(string_pos[[i]]) != length(motif_location[[i]]) ){
    print(paste0("not matching at: ",i))
  }
}

# binding lists into dataframe
df <- do.call(rbind,lapply(1:length(motif_location), function(i)
  data.frame(motif_num = i,
            string_pos=unlist(string_pos[[i]]),
            motif_loc=unlist(motif_location[[i]]))))

# subsetting
motif_diagram <-  motif_diagram[,c(1,2)]
# adding row number
motif_diagram$motif_num <- seq.int(nrow(motif_diagram))

# merging motif_diagram with df (string position and motif location)
motif_diagram <-  merge(motif_diagram, df, by = 'motif_num')
motif_diagram$motif_loc <-  gsub('\\]\\-[0-9]+','\\]', motif_diagram$motif_loc)

# manipulating and splitting columns to get motif and strand
motif_diagram$motif_loc <-  gsub("\\[\\-","\\[\\-\\:", motif_diagram$motif_loc, perl = T)
motif_diagram$motif_loc <-  gsub("\\[\\+","\\[\\+\\:", motif_diagram$motif_loc, perl = T)
setDT(motif_diagram)[,c('strand','motif'):=tstrsplit(motif_diagram$motif_loc,':')]
motif_diagram$strand <- gsub('\\[','',motif_diagram$strand)
motif_diagram$motif <- gsub('\\]','',motif_diagram$motif)

# subsetting
motif_diagram <- motif_diagram[,c(2:4,6,7)]
motif_diagram[is.na(motif_diagram)] <- 'EOS' 

setDT(motif_diagram)[,c('V1','V2','V3'):= tstrsplit(motif_diagram$seq_name,'_', fixed = TRUE)]
motif_diagram <-  motif_diagram[,c(6,1:5)]

# adding motif_info to motif_diagram
# subset motif_info
motif_info_subset <-  motif_info[,c(1,3,4)]

# create dictionary: key = motif# value = motif width
dict_motifInfo <- motif_info_subset$WIDTH
names(dict_motifInfo) <- motif_info_subset$MOTIF
# adding motif_width to motif_diagram df
for(i in 1:length(motif_diagram$motif)){
  if(motif_diagram$motif[i] != 'EOS'){
    motif_diagram$motif_width[i] <- dict_motifInfo[[motif_diagram$motif[i]]]
  } else{
    motif_diagram$motif_width[i] <- ''
  }

}

# splitting column V1 to get gene coordinates
setDT(motif_diagram)[,c('misc','coordinates'):=tstrsplit(motif_diagram$V1,'=')]
setDT(motif_diagram)[,c('Chr','start','end'):=tstrsplit(motif_diagram$coordinates,':|-')]

# subset
motif_diagram <-  motif_diagram[,c(1:7,10:12)]

backup <- motif_diagram
motif_diagram <-  backup

# 5. calculating actual motif locations -----------------
motif_diagram$string_pos <-  as.character(motif_diagram$string_pos)
motif_diagram$string_pos <-  as.numeric(motif_diagram$string_pos)
motif_diagram$motif_width <-  as.numeric(motif_diagram$motif_width)

motif_diagram <- transform(motif_diagram, start = as.numeric(start))
motif_diagram <- transform(motif_diagram, end = as.numeric(end))

motif_diagram$motif_start[1] <- motif_diagram$start[1] + motif_diagram$string_pos[1]
motif_diagram$motif_end[1] <- motif_diagram$motif_start[1] + motif_diagram$motif_width[1]



for(x in 2:nrow(motif_diagram)){
  if(motif_diagram$motif[x] == 'EOS') {
    print(paste0('row with EOS:', x))
    motif_diagram$motif_start[x] <- 0
    motif_diagram$motif_end[x] <- 0
    
  } else if(motif_diagram$motif[x-1] == 'EOS'){
    print(paste0('row after EOS:', x))
    motif_diagram$motif_start[x] <- motif_diagram$start[x] + motif_diagram$string_pos[x]
    motif_diagram$motif_end[x] <- motif_diagram$motif_start[x] + motif_diagram$motif_width[x]
  }
  else{
    motif_diagram$motif_start[x] <- motif_diagram$motif_end[x-1] + motif_diagram$string_pos[x]
    motif_diagram$motif_end[x] <- motif_diagram$motif_start[x] + motif_diagram$motif_width[x]

    }
}

#write.table(motif_diagram, file = paste0(Sys.Date(),'_motif_start_stop_locations.txt'), col.names = T, row.names = F, sep = '\t', quote = F)

# 6. converting seq.name to refseq ID -------------------------------
seqname_refseqID_mapping <-  read.delim('~/KP/H3K27me3/motif-analysis/46_candidate genes/seqName_refseqID.txt', header = F, sep = ';')
seqname_refseqID_mapping <-  seqname_refseqID_mapping[!duplicated(seqname_refseqID_mapping),]

# merging seq name from motif_diagram with refseq IDS
motif_refseq <- merge(motif_diagram, seqname_refseqID_mapping, by.x = 'V1', by.y = 'V2', allow.cartesian=TRUE)
motif_refseq <-  motif_refseq[!duplicated(motif_refseq),]
# renaming columns
names(motif_refseq)[13] <- 'refseq_ID'
# dropping a column
motif_refseq$seq_name <-  NULL
# renaming column
names(motif_refseq)[1] <- 'seq_name'

# 7. Mapping gene symbols, start, stop, distance.to.TSS to motifs ----------------------------
motif_refseq <-  merge(motif_refseq, refGene_TSS, by.x = 'refseq_ID', by.y = 'refGene_ID', all = TRUE, allow.cartesian=TRUE)
motif_refseq <-  na.omit(motif_refseq) # 21/46 genes have been found to have motif binding
# subsetting
motif_refseq <-  motif_refseq[,c(-14,-15,-18,-19)]
# renaming columns
names(motif_refseq)[5] <- 'motif_strand'
names(motif_refseq)[15] <- 'strand'
names(motif_refseq)[8] <-  'Chr'
names(motif_refseq)[9] <-  'start'

#write.table(motif_refseq, file = paste0(Sys.Date(),'_motif_refseq_details.txt'), col.names = T, row.names = F, sep = '\t', quote = F)

