# super enhancer lists have been divided into separate cell line files
# create pairwise overlaps such that a gene is present in atleast two cell lines
# create a union list of all such individual pairwise comparisons

amp_cl <- c('COGN415','KELLY', 'NB1643', 'LAN5', 'NGP', 'SKNBE2C')
nonAmp_cl <- c('NB69', 'NBLS', 'SKNAS', 'SKNFI')

# ------------ FOR AMPLIFIED LINES -----------------
# reading in individual cell line SE files
for(x in amp_cl){
  df <- x
  assign(df,read.delim(paste0('/Volumes/target_nbl_ngs/KP/enhancer-rank-list-H3K27Ac/',x,'-H3K27Ac-SEs-anno-rank-sorted.bed'), header = T))
}

pairwise_ampLst <- c()

# loop for pairwise comparisons
for(i in 1:length(amp_cl)){
  for(j in 1:length(amp_cl)){

    intersect_cl <- paste0(amp_cl[i],'_',amp_cl[j])
    
    # eliminating same cell line comparisons
    if(i != j){
      print(paste0('Processing...', amp_cl[i],'_', amp_cl[j]))
      
      # appending to list
      pairwise_ampLst <- c(pairwise_ampLst, paste0(amp_cl[i],'_', amp_cl[j]))
      
      # fetching the dfs
      one <- get(amp_cl[i])
      two <- get(amp_cl[j])
      
      # intersecting
      assign(intersect_cl, intersect(one[,6], two[,6]))

    }
  } 
}

# creating unique union list
l1 <- lapply(pairwise_ampLst, function(x) get(x))
union_l1 <- unique(Reduce(union, l1))

write.table(union_l1, file = paste0(Sys.Date(),'_ampLines_SE_pairwise_list.txt'), col.names = F, quote = F, row.names = F)

# ------------ FOR NON-AMPLIFIED LINES -----------------
# reading in individual cell line SE files
for(x in nonAmp_cl){
  df <- x
  assign(df,read.delim(paste0('/Volumes/target_nbl_ngs/KP/enhancer-rank-list-H3K27Ac/',x,'-H3K27Ac-SEs-anno-rank-sorted.bed'), header = T))
}


pairwise_nonAmpLst <- c()

for(i in 1:length(nonAmp_cl)){
  for(j in 1:length(nonAmp_cl)){
    
    intersect_cl <- paste0(nonAmp_cl[i],'_',nonAmp_cl[j])
    
    if(i != j){
      print(paste0('Processing...', nonAmp_cl[i],'_', nonAmp_cl[j]))
      
      # appending to list
      pairwise_nonAmpLst <- c(pairwise_nonAmpLst, paste0(nonAmp_cl[i],'_', nonAmp_cl[j]))
      
      # fetching the dfs
      one <- get(nonAmp_cl[i])
      two <- get(nonAmp_cl[j])
      
      # intersecting
      assign(intersect_cl, intersect(one[,6], two[,6]))
      
    }
  } 
}


# creating unique union list
l2 <- lapply(pairwise_nonAmpLst, function(x) get(x))
union_l2 <- unique(Reduce(union, l2))

write.table(union_l2, file = paste0(Sys.Date(),'_nonAmpLines_SE_pairwise_list.txt'), col.names = F, quote = F, row.names = F)
