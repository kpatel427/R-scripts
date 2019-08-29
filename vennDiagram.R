library(VennDiagram)
library(ggplot2)

# for mycn-amp 
# reading in the gene lists
mycn.amp.list <- read.delim('2019-08-26-MYCN.amp.genelist.txt', header = F)
mycn.Nonamp.list <- read.delim('2019-08-26-MYCN.Nonamp.genelist.txt', header = F)

versteeg.list <- read.delim('Versteeg_paper_gene_list.txt', header = F)
boeva.list <- read.delim('Boeva_paper_gene_list.txt', header = F)

# to get length of lists for each area
mycnAmp_counts <- length(mycn.amp.list$MYCN_amplified)
versteeg_counts <- length(versteeg.list$Versteeg)
boeva_counts <- length(boeva.list$Boeva)

mycnAmp_versteeg_counts <- length(unique(intersect(mycn.amp.list$MYCN_amplified, versteeg.list$Versteeg)))
versteeg_boeva_counts <- length(unique(intersect(versteeg.list$Versteeg, boeva.list$Boeva)))
mycnAmp_boeva_counts <- length(unique(intersect(mycn.amp.list$MYCN_amplified, boeva.list$Boeva)))
mycnAmp_versteeg_boeva_counts <- length(Reduce(function(x, y) unique(intersect(x, y)), list(mycn.amp.list$MYCN_amplified, versteeg.list$Versteeg, boeva.list$Boeva)))

# plotting venn for mycn amp
grid.newpage()
venn.plot <- draw.triple.venn(area1 = mycnAmp_counts, 
                 area2 = versteeg_counts, 
                 area3 = boeva_counts, 
                 n12 = mycnAmp_versteeg_counts, 
                 n23 = versteeg_boeva_counts, 
                 n13 = mycnAmp_boeva_counts, 
                 n123 = mycnAmp_versteeg_boeva_counts, 
                 category = c("MYCN Amplified gene list", "Versteeg Gene list", "Boeva Gene list"), 
                 lty = "blank",
                 scaled = TRUE,
                 cat.fontfamily = "sans",
                 fontfamily = 'sans',
                 alpha = c(0.4,0.4,0.4),
                 fill = c("skyblue", "mediumorchid", "orange"),
                 cat.cex = 1.2, # category name size
                 cex = 0.9, # areas' labels font size
                 resolution = 300,
                 output = TRUE,
                 height = 12,
                 width = 10,
                 units = 'in'
);

#grid.draw(venn.plot);
# saving plot
ggsave(venn.plot, file="MYCNAmp_versteeg_boeva.tiff", device = "tiff", width = 12, height = 10)



# for mycn non-amp

# calculating lengths of list
mycnNonamp_counts <- length(mycn.Nonamp.list$V1)

mycnNonamp_versteeg_counts <- length(unique(intersect(mycn.Nonamp.list$V1, versteeg.list$Versteeg)))
versteeg_boeva_counts <- length(unique(intersect(versteeg.list$Versteeg, boeva.list$Boeva)))
mycnNonamp_boeva_counts <- length(unique(intersect(mycn.Nonamp.list$V1, boeva.list$Boeva)))
mycnNonamp_versteeg_boeva_counts <- length(Reduce(function(x, y) unique(intersect(x, y)), list(mycn.Nonamp.list$V1, versteeg.list$Versteeg, boeva.list$Boeva)))

# plotting venn for mycn non-amp
grid.newpage()
venn.plot <- draw.triple.venn(area1 = mycnNonamp_counts, 
                              area2 = versteeg_counts, 
                              area3 = boeva_counts, 
                              n12 = mycnNonamp_versteeg_counts, 
                              n23 = versteeg_boeva_counts, 
                              n13 = mycnNonamp_boeva_counts, 
                              n123 = mycnNonamp_versteeg_boeva_counts, 
                              category = c("MYCN Non-amp gene list", "Versteeg Gene list", "Boeva Gene list"), 
                              lty = "blank",
                              #lty = 1,
                              scaled = TRUE,
                              cat.fontfamily = "sans",
                              fontfamily = 'sans',
                              alpha = c(0.4,0.4,0.4),
                              fill = c("skyblue", "mediumorchid", "orange"),
                              cat.cex = 1.2, # category name size
                              margin = 0.05,
                              cex = 0.9, # areas' labels font size
                              resolution = 300,
                              output = TRUE,
                              height = 12,
                              width = 10,
                              units = 'in'
);

#grid.draw(venn.plot);
#dev.off()

# saving plot
ggsave(venn.plot, file="MYCN_Nonamp_versteeg_boeva.tiff", device = "tiff", width = 12, height = 10)
