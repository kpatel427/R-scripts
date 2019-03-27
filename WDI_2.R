set.seed(12345)
library(ggplot2)
library(ggrepel)
library(reshape)
library(dplyr)

# reading data
dat1 <- read.csv("~/Desktop/DISH/WDI_csv/WDIData.csv")
#colnames(dat1)
dat1[is.na(dat1)] = 0

patterns <- c("^SP.POP.*.FE.ZS$","^SP.POP.*.MA.ZS$")
# for total population ages - extracting using indicator code
q2 <- dat1[grep(paste(patterns,collapse = "|"), dat1$Indicator.Code),]
q2[is.na(q2)] = 0
q2$gender <- ifelse(grepl( "^SP.POP.*.FE.ZS$",q2$Indicator.Code),"F","M")


# Re-structring data
#melting the dataframe and subsetting
tmp <- melt(q2)
tmp.2 <- subset(tmp, tmp$variable == "X2016", select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value","gender"))


ggplot(tmp.2, aes(Country.Code,value, fill = factor(tmp.2$Indicator.Name) )) + 
  geom_bar(aes(fill = gender), stat = "identity", position = "dodge") +
  coord_flip()+ 
  #geom_text_repel(aes(label = ifelse( tmp.2$label > 0 ,paste0(tmp.2$Country.Name),"")),size = 2,point.padding = 1,segment.color = "black", show.legend = FALSE) + 
  theme_bw() +
  labs(title = "Age distribution across all age groups in 2016", x = "female vs male", y = "Country") +
  theme(axis.text = element_text(size=5))

# since it is impossible to visualize which countries have the largest difference 
# it is better to cut down & visualize the top 50 countries with largest difference

df <- subset(tmp.2,tmp.2$Indicator.Name != "Population, female (% of total)")
df <- subset(df,df$Indicator.Name != "Population, male (% of total)")

df <- df %>%
  group_by(Country.Name) %>%
  #arrange(Indicator.Name) %>%
  mutate(diff = value - lag(value, default = first(value)))

for(i in seq(1, nrow(df), by = 2)){ 
  df$diff[i] = 0
  df$diff = abs(df$diff)
}

test <- head(df[with(df,order(-diff)),],50)

subset.tmp.2 <- subset(tmp.2,tmp.2$Country.Code %in% test$Country.Code, select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value","gender"))

labels <- c(F = "Female", M = "Male")

# only plotting top 50 countries with largest difference in male and female populations
p <- ggplot(subset.tmp.2, aes(Country.Name,value, fill = factor(Indicator.Name) )) + 
  geom_bar(aes(fill = Indicator.Name), stat = "identity", position = "dodge") +
  coord_flip()+ 
  #geom_text_repel(aes(label = ifelse( tmp.2$label > 0 ,paste0(tmp.2$Country.Name),"")),size = 2,point.padding = 1,segment.color = "black", show.legend = FALSE) + 
  theme_bw() +
  labs(title = "Age distribution across all age groups in 2016", x = "Country", y = "Percent Population") +
  theme(axis.text = element_text(size=10)) +
  facet_wrap(gender ~., scales = "free", labeller = labeller(gender = labels))


# Save as PDF
ggsave(paste0(Sys.Date(),"-Q_2",".pdf"), p, width=40, height=20, device = "pdf", units = "cm")
