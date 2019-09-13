# script to extend RColorbrewer palette

library(RColorBrewer)


colourCount = length(unique(mtcars$hp))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

ggplot(mtcars) + 
  geom_histogram(aes(factor(hp)), fill=getPalette(colourCount)) + 
  theme(legend.position="right")


ggplot(mtcars) + 
  geom_histogram(aes(factor(hp), fill=factor(hp)), stat = 'count') + 
  scale_fill_manual(values = getPalette(colourCount))
