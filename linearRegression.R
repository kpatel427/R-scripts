library(data.table)
library(tidyverse)
library(dplyr)
library(reshape)
library(reshape2)
library(matrixStats)
library(ggpubr)
library(lme4)
library(broom)


# ..LINEAR REGRESSION ----------------------------------------------------------------------------------
# PFS vs expression
# PFS = response variable
# expression = dependent variable

# To predict response with expression data


# to perform ULR for every group individually
fits <- lmList(PFS ~ pdx_FPKM | gene, data=pdx_sensitive)
# summary(fits$AACSP1)
# summary(fits$AARS)
# summary(fits$AARS2)



# # function definition
# lmStats <- function(x) {
#   # getting stats into a dataframe
#   df.res <- data.frame(gene = names(fits)[x])
#   
#   df1 <- tidy(summary(fits[[x]])) # coefficients and p vals
#   df2 <- glance(summary(fits[[x]]))[1:2] # r.squared and adjusted r.squared
#   
#   df.res <-  cbind(df.res, df1, df2)
#   
#   return(df.res)
# }


# to extract regression stats from model fits
final <- data.frame()
for(x in 1:length(unique(pdx_sensitive$gene))) {
  
  df.res <- data.frame(gene = names(fits)[x])
  
  df1 <- tidy(summary(fits[[x]])) # coefficients and p vals
  df2 <- glance(summary(fits[[x]]))[1:2] # r.squared and adjusted r.squared
  
  df.res <-  cbind(df.res, df1, df2)
  
  final <- rbind(final, df.res)
}


#write.table(final, file = paste0(Sys.Date(),'_linearRegression_PFS_expression_invivo_sensitive_DEgenes.txt'), col.names = T, row.names = F, sep = '\t', quote = F)


# filtering for genes with p.val =< 0.025 and adjusted r2 >= 0.8 (less stringent)
final_filter_ls <- final[(final$p.value <= 0.025 & final$adj.r.squared >= 0.8),]
write.table(final_filter_ls, file = paste0(Sys.Date(),'_linearRegression_PFS_expression_invivo_sensitive_DEgenes_filter_less_stringent.txt'), col.names = T, row.names = F, sep = '\t', quote = F)


# filtering for genes with p.val =< 0.025 and adjusted r2 >= 0.9 (more stringent)
final_filter_ms <- final[(final$p.value <= 0.025 & final$adj.r.squared >= 0.9),]
write.table(final_filter_ms, file = paste0(Sys.Date(),'_linearRegression_PFS_expression_invivo_sensitive_DEgenes_filter_more_stringent.txt'), col.names = T, row.names = F, sep = '\t', quote = F)

