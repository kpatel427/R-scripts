# function to normalize/scale values between any range

# Parameters
# df = data frame
# x = name of the column to be scaled
# a = start of a range
# b = end of a range
# returns a df one column containing scaled values

normalizeValues <- function(df, x, a, b){
  min_logfc <- min(df[,x])
  max_logfc <- max(df[,x])
  
  scaled_log2FC <- as.data.frame(sapply(df[,x], function(x) ( (b-a) * (x - min_logfc)/(max_logfc - min_logfc)) + a ))
  return(scaled_log2FC)
}
