#To provide a list of mouse genes

library(biomaRt)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listFilters(mouse)

getBM(attributes=c("ensembl_gene_id", "mgi_symbol"), filters= "mgi_symbol", mart=mouse)