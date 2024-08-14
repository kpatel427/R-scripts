# script to get all genes orthologous to human in mouse

library(tidyverse)
library(homologene)

# Retrieve all orthologs between human and mouse
all_human_mouse_orthologs <- homologeneData2 |> 
  filter(Taxonomy %in% c(9606, 10090))

# Split human and mouse data
human_genes <- all_human_mouse_orthologs |>  filter(Taxonomy == 9606)
mouse_genes <- all_human_mouse_orthologs |>  filter(Taxonomy == 10090)


# Merge human and mouse data by HID to get ortholog pairs
human_mouse_pairs <- merge(human_genes, mouse_genes, by = "HID", suffixes = c("_human", "_mouse"))
