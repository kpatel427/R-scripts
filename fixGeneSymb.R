# function to replace gene symbols once they have been opened in excel


replaceSymb <- function(df){
  df$gene <- gsub('1-Mar','MAR1',df$gene)
  df$gene <- gsub('10-Sep','SEPT10',df$gene)
  df$gene <- gsub('11-Mar','MARCH11',df$gene)
  df$gene <- gsub('11-Sep','SEPT11',df$gene)
  df$gene <- gsub('15-Sep','SEP15',df$gene)
  df$gene <- gsub('2-Mar','MAR2',df$gene)
  df$gene <- gsub('2-Sep','SEP2',df$gene)
  df$gene <- gsub('3-Sep','SEP3',df$gene)
  df$gene <- gsub('5-Mar','MAR5',df$gene)
  df$gene <- gsub('6-Mar','MAR6',df$gene)
  df$gene <- gsub('6-Sep','SEPT6',df$gene)
  df$gene <- gsub('7-Sep','SEPT7',df$gene)
  df$gene <- gsub('7-Mar','MAR7',df$gene)
  df$gene <- gsub('8-Mar','MAR8',df$gene)
  df$gene <- gsub('8-Sep','SEPT8',df$gene)
  df$gene <- gsub('9-Mar','MAR9',df$gene)
  df$gene <- gsub('9-Sep','SEPT9',df$gene)
  
  return(df)
}
