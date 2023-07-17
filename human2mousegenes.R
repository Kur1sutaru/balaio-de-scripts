#convert human to mouse genes 
setwd()
library(biomaRt)
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  no_mouse_genes <- length(x)
  no_human_genes <- length(humanx)
  
  if(no_human_genes != no_mouse_genes){
    print("Some genes could not be translated!")
    genes_not_trans <- setdiff(x,genesV2$HGNC.symbol)
    print("These genes could not be translated:")
    print(genes_not_trans)
    print(paste("A total number of ",length(genes_not_trans),"genes could not be translated!"),sep=" ")
  }else{
    print("All genes were translated successfully!")
  }
  
  # Print all gene names that could not be translated and the number of genes that were not translated
  
  return(humanx)
}