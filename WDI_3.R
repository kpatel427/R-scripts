set.seed(12345)
library(ggplot2)
library(ggrepel)
library(reshape)
library(dplyr)
library(ggpubr)
library(tidyr)

# reading data
dat1 <- read.csv("~/Desktop/DISH/WDI_csv/WDIData.csv")
#colnames(dat1)
dat1[is.na(dat1)] = 0

patterns <- c("^IT.CEL.SETS.P2*","^IT.MLT.MAIN.P2","^SP.POP.*.FE.ZS$","^SP.POP.*.MA.ZS$")
# for total population ages - extracting using indicator code
q3 <- dat1[grep(paste(patterns,collapse = "|"), dat1$Indicator.Code),]

# Re-structring data
#melting the dataframe and subsetting
tmp <- melt(q3)
tmp.2 <- subset(tmp, tmp$variable == "X2016", select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value"))

tmp.2 <- subset(tmp.2,tmp.2$Indicator.Name != "Population, male (% of total)")
tmp.2 <- subset(tmp.2,tmp.2$Indicator.Name != "Population, female (% of total)")

df1 <- subset(tmp.2, tmp.2$Indicator.Name == "Fixed telephone subscriptions (per 100 people)" | tmp.2$Indicator.Name == "Mobile cellular subscriptions (per 100 people)"
              , select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value"))

df1 <- cast(df1,Country.Name + Country.Code + value ~ Indicator.Name)
df1[is.na(df1)] <- 0


df1.TL <- subset(df1,df1$`Fixed telephone subscriptions (per 100 people)` != 0,select = c("Country.Name","Country.Code","Fixed telephone subscriptions (per 100 people)") )
df1.Cell <- subset(df1,df1$`Mobile cellular subscriptions (per 100 people)` != 0,select = c("Country.Name","Country.Code","Mobile cellular subscriptions (per 100 people)") )

df2 <- subset(tmp.2, grepl("^SP.POP.*.ZS$",tmp.2$Indicator.Code), select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value"))

# merging dataframes
merge.df1.cell <- merge(df1.Cell,df2, by = "Country.Name")
colnames(merge.df1.cell)[colnames(merge.df1.cell) == 'Mobile cellular subscriptions (per 100 people)'] <- 'cellular.subscriptions'

merge.df1.TL <- merge(df1.TL,df2, by = "Country.Name")
colnames(merge.df1.TL)[colnames(merge.df1.TL) == 'Fixed telephone subscriptions (per 100 people)'] <- 'telephone.subscriptions'


# to see if there is an association between age distribution and mobile phone adoption

#install.packages("corrplot")
library(corrplot)

merge.df1.cell <- cast(merge.df1.cell, Country.Name + cellular.subscriptions + value ~ Indicator.Name)
merge.df1.cell <- subset(merge.df1.cell, select = -c(value))
merge.df1.cell[is.na(merge.df1.cell)] <- 0

# generating a correlation matrix

M<-cor(merge.df1.cell[,2:8])

# for headers ion corr Matrix M
names<- colnames(merge.df1.cell)
colnames(M) <- names[2:ncol(merge.df1.cell)]
rownames(M) <- names[2:ncol(merge.df1.cell)]

# plotting corr matrix
corrplot(M, method = "circle", type = "upper")


# to see if there is an association between age distribution and telephone adoption
merge.df1.TL <- cast(merge.df1.TL, Country.Name + telephone.subscriptions + value ~ Indicator.Name)
merge.df1.TL <- subset(merge.df1.TL, select = -c(value))
merge.df1.TL[is.na(merge.df1.TL)] <- 0

# generating a correlation matrix

M1<-cor(merge.df1.TL[,2:8])

# for headers ion corr Matrix M
names<- colnames(merge.df1.TL)
colnames(M1) <- names[2:ncol(merge.df1.TL)]
rownames(M1) <- names[2:ncol(merge.df1.TL)]

# plotting corr matrix
corrplot(M1, method = "circle", type = "upper")
