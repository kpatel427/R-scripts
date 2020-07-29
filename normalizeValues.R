# function to normalize/scale values between any range

# Parameters
# x = name of the column to be scaled
# a = start of a range
# b = end of a range

normalizeValues <- function(x, a, b){
  min_logfc <- min(model.hm[,x])
  max_logfc <- max(model.hm[,x])
  
  scaled <- as.data.frame(sapply(model.hm[,x], function(x) ( (b-a) * (x - min_logfc)/(max_logfc - min_logfc)) + a ))
  return(scaled)
}
