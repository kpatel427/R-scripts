library(boxr)

box_auth() 

box_ls(dir_id = 87809132620) # NBSD GPC2 CAR escaped tumors_RNA seq_9_2019 folder
box_ls(dir_id = 87943295083) # 2019-09-23-RNAseq-Fastq-files folder 

# set working dir
box_setwd(dir_id = 87943295083)

# download specific files
files <- box_search(query = "_001.fastq.gz", ancestor_folder_ids = 87943295083)

files[[1]]$id

files_to_download <- c()

for(x in 1:length(files)){
  print('Getting file ids....')
  
  files_to_download <- c(files_to_download, ifelse(grep("_001.fastq.gz",files[[x]]$name, perl = T), paste0(files[[x]]$id,':',files[[x]]$name), ''))
 
}


# download files
for(x in files_to_download){
  
  id = strsplit(x,':')[[1]][1]
  file_name = strsplit(x,':')[[1]][2]
  
  print(paste0('Processing ',file_name))
  
  # downloading files
  file_name <- as.data.frame(box_dl(file_id = id))
  write.table(get(file_name), file = paste0('/mnt/isilon/maris_lab/target_nbl_ngs/KP/CAR_T_RNASeq/',file_name))
  
  print(paste0(file_name,' file written!'))
  
}




