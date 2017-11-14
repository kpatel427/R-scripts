# reads data from CSV files and plots and overlays distribution using ggplot2 package
# Each CSV file has different number of rows


fill <- "#4271AE"
line <- "#1F3552"

data1 <- read.csv("BS.csv", header = TRUE)
data2 <- read.csv("PF.csv", header = TRUE)
data3 <- read.csv("EC.csv", header = TRUE)
data4 <- read.csv("PP.csv", header = TRUE)
data5 <- read.csv("VH.csv", header = TRUE)
data6 <- read.csv("VC.csv", header = TRUE)
data7 <- read.csv("PA.csv", header = TRUE)

library(ggplot2);library(reshape2)

mean_d1=mean(data5$diff) 

# merge the two data sets
d <- rbind(data1, data2, data3, data4)
d$diff <- as.factor(d$diff)

# update the aesthetic
p1 <- ggplot(data = data1, aes(x=diff))+
  geom_density() +
  # ----------------------------
  geom_density(data=data2,position = "stack", fill="purple",alpha=0.25) +
  labs(title = "Distribution of CpGs Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  # ----------------------------
  geom_density(data=data3, position = "stack", fill="green",alpha=0.25) +
  labs(title = "Distribution of CpGs Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  # ----------------------------
  geom_density(data=data4,position = "stack", fill="blue",alpha=0.25) +
  labs(title = "Distribution of CpGs Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  # ----------------------------

  geom_density(data=data6,position = "stack", fill="yellow",alpha=0.25) +
  labs(title = "Distribution of CpGs Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  # ----------------------------
  geom_density(data=data5,position = "stack", fill="cyan",alpha=0.25) +
  labs(title = "Distribution of CpGs Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  # ----------------------------
  geom_density(data=data7,position = "stack", fill="black",alpha=0.25) +
  labs(title = "Distribution of CpG Lengths")+
  labs(y="Density")+
  labs(x="CpG Lengths")+
  xlim(0,6000)


plot(p1)
