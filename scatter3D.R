set.seed(1234)

library(ggplot2)
library(tidyverse)
library(data.table)
library(scatterplot3d)


# ....plotting  ---------------------------------------
final.merged$group <- as.factor(final.merged$group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(final.merged$group)]


scatterplot3d(final.merged$rank, y=final.merged$p.value, z=final.merged$slope, 
              angle = 55,
              color = colors,
              grid=TRUE, 
              box=T,
              main="3D Plot for PDX models",
              xlab = "Rank",
              ylab = "P-value",
              zlab = "Slope")

legend("bottom", legend = levels(final.merged$group),
       col =  c("#E69F00", "#56B4E9"),
       inset = -0.55, 
       xpd = T,
       par(mar=c(5, 4, 4, 2) + 0.1),
       horiz = TRUE,
       bty = "n",
       pt.cex = 1,
       lty = 1,
       pch=1, cex = 0.75,
       x.intersp = 0.1,
       text.width = 0.75)


dev.copy(pdf,'scatter3DPlot_top10K_binarySlope_invivo_PDX.pdf')
dev.off()
