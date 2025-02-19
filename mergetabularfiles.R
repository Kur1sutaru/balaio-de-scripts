

library(dplyr)
library(tidyr)

# Listar todos os arquivos .tabular na pasta
file_list <- list.files(
  path = "C:/Users/Cristal/Downloads/1. GSE162694 - PRONTO PRO R-20250219T221835Z-001/1. GSE162694 - PRONTO PRO R",
  pattern = "*.tabular", full.names = TRUE
)

# Criar uma lista para armazenar os dataframes
df_list <- lapply(file_list, function(f) {
  df <- read.delim(f, header = TRUE)  # Ler o arquivo tabular
  
  # Garantir que há pelo menos 2 colunas
  if (ncol(df) >= 2) {
    df <- df[, 1:2]  # Seleciona apenas as duas primeiras colunas
    colnames(df) <- c("entrez_id", basename(f))  # Define os nomes de colunas
    
    # Renomear a coluna do arquivo sem a extensão ".tabular"
    colnames(df)[2] <- gsub("\\.tabular$", "", colnames(df)[2])
    return(df)
  } else {
    warning(paste("O arquivo", f, "não possui colunas suficientes!"))
    return(NULL)
  }
})

# Unir os arquivos por "entrez_id" (Join Progressivo)
merged_df <- Reduce(function(x, y) full_join(x, y, by = "entrez_id"), df_list)

# Visualizar os primeiros registros
head(merged_df)
write.csv(merged_df, "GSE162694alltabular.csv")
