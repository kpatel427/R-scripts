# script to download files from a folder on Box using Box api
library(boxr)

box_auth() 

box_ls(dir_id = 8xxxxxx)
box_ls(dir_id = 8xxxxxx) 

# set working dir
box_setwd(dir_id = 10xxxxxx)

# download specific files
files <- box_search(query = "_001.fastq.gz", ancestor_folder_ids = 10xxxxxx)

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
  file_name <- box_dl(file_id = id, local_dir = getwd())
  
  print(paste0(file_name,' file written!'))
  
}

