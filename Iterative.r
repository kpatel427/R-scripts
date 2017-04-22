#** INSTRUCTIONS **#
#----Create a 21717 folder in your home directory (~/)
#----Change the path in the first line (read.csv) to your working directory.
#----Make sure that you have set your current directory as your working directory beforeyou start working.
#----Change the number of columns in the range x <- (1:total number of cols). Check master file to find out the total number of columns
#----Your files will be created in the same directory as you are currently working


data1 <- read.csv("~/R/21717/master_21717.csv")
x <- (1:18)


for (i in x)
{
  if (i %% 2 == 0)next  #skipping alternate iterations
  a <- i+1        
 
  print(i)
  print(a)
  
  x = colnames(data1)[i]
  y = colnames(data1)[a] 

  # write these cols in a new csv
  seq1 <- data1[[x]]
  seq2 <- data1[[y]]
  
  
  
  myfile <- file.path("~/R/21717", paste0(y, ".csv"))
  
  # creating new temp CSVs
  X <- data.frame(seq1, seq2)
  write.table(X, file = myfile, row.names=FALSE,col.names=c("x","y"), append = TRUE, na = "", sep = ",")
  
  
  
  #----------------------------SEPERATE DISTRIBUTIONS------------------------------#
  data <- read.csv(myfile)
  data<-na.omit(data)
  
  
  
  library(mixtools)
  x11()
  
  wait = data$y
  mixmdl = normalmixEM(wait, k=3, epsilon = 1e-08,  maxit = 1000)
  plot(mixmdl,which=2)
  lines(density(wait), lty=2, lwd=2)

  #mypath <- file.path("~/R/21717",paste("myplot_", y, ".jpg", sep = ""))
  
  #jpeg(file=mypath)
  
  
}









