# script to get Ensemble gene ID's for gene symbols
library(biomaRt)

df <-  read.delim('genelist.txt', header = F)
genes <- df$V1

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#filters: is a vector of filters that one wil use as input to the query.
filters = listFilters(ensembl)
filters[1:5,]

# attributes: is a vector of attributes that one wants to retrieve (= the output of the query).
attributes = listAttributes(ensembl)
attributes[1:5,]

geneIDs <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', 'chromosome_name',
                              'start_position', 'end_position'), 
      filters = 'hgnc_symbol', 
      values = genes, 
      mart = ensembl)


write.table(geneIDs, file = "SE_filtered_gene_biomaRT.txt", quote = F, col.names = F, row.names = F)
