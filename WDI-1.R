set.seed(12345)
library(ggplot2)
library(ggrepel)
library(reshape)
library(dplyr)

# reading data
dat1 <- read.csv("~/Desktop/DISH/WDI_csv/WDIData.csv")
#colnames(dat1)
dat1[is.na(dat1)] = 0

# for total population ages - extracting using indicator code
q1 <- dat1[grep('^SP.POP.*.TO.ZS$', dat1$Indicator.Code),]
q1[is.na(q1)] = 0

# Re-structring data
#melting the dataframe and subsetting
tmp <- melt(q1)
tmp.2 <- subset(tmp, tmp$variable == "X2016", select = c("Country.Name","Country.Code","Indicator.Name","Indicator.Code","variable","value"))
tmp.2$label <- 0

#colnames(tmp.2)

# for labelling outliers for each group
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

p <- tmp.2 %>%
  group_by(Indicator.Name) %>%
  mutate(outlier = ifelse(is_outlier(value), paste0(tmp.2$Country.Name), as.numeric(NA))) %>%
  ggplot(., aes(x = factor(Indicator.Name), y = value, fill = Indicator.Name)) +
  geom_boxplot(outlier.size=5, alpha=0.1) +
  geom_text_repel(aes(label = outlier),nudge_y = 0.5, na.rm = TRUE, hjust = 0.5, size = 2, segment.size = 0.3) +
  theme_bw() +
  labs(title = "Age distribution across all age groups in 2016", x = "Population divided into age groups", y = "count")


# Save as PDF
ggsave(paste0(Sys.Date(),"-Q_1",".pdf"), p, width=40, height=20, device = "pdf", units = "cm")




