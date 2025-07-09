setwd("C:/Users/Cristal/Downloads/pati_mirna")
# Load packages
library(tidyverse)
library(limma)
library(ggplot2)
library(pheatmap)

# Read the long-format table
df <- read.csv("C:/Users/Cristal/Downloads/pati_mirna/mirna_expression_long_annotated.csv")  # This file must have columns: Feature.ID, Expression.values, sample, group

# Optional: inspect structure
str(df)

# Convert long format to wide expression matrix
expr_matrix <- df %>%
  select(Feature.ID, sample, Expression.values) %>%
  pivot_wider(names_from = sample, values_from = Expression.values) %>%
  column_to_rownames("Feature.ID")

write.csv(expr_matrix, "expr_matrix_mirna_countswithannotations.csv")

# Create metadata
sample_metadata <- df %>%
  select(sample, group) %>%
  distinct() %>%
  column_to_rownames("sample")

# Ensure column order matches
expr_matrix <- expr_matrix[, rownames(sample_metadata)]

# Create design matrix
design <- model.matrix(~0 + sample_metadata$group)
colnames(design) <- gsub("sample_metadata\\$group", "", colnames(design))

# Apply voom normalization and linear modeling
v <- voom(expr_matrix, design)
fit <- lmFit(v, design)

# Define comparisons
contr.matrix <- makeContrasts(
  DCD_vs_DBD = pgd_dcd - pgd_dbd,
  PGD_DCD_vs_nonPGD_DCD = pgd_dcd - non_pgd_dcd,
  PGD_DBD_vs_nonPGD_DBD = pgd_dbd - non_pgd_dbd,
  levels = design
)

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

# Get top differentially expressed miRNAs
top_dcd_vs_dbd <- topTable(fit2, coef = "DCD_vs_DBD", number = Inf)
top_dcd_vs_nonpgd <- topTable(fit2, coef = "PGD_DCD_vs_nonPGD_DCD", number = Inf)
top_dbd_vs_nonpgd <- topTable(fit2, coef = "PGD_DBD_vs_nonPGD_DBD", number = Inf)


# Heatmap of top miRNAs
top_genes <- rownames(top_dcd_vs_dbd)[1:30]
pheatmap(expr_matrix[top_genes, ], annotation_col = sample_metadata)


write.csv(top_dcd_vs_dbd, "DCD_vs_DBD_DEGs.csv")
write.csv(top_dcd_vs_nonpgd, "PGD_DCD_vs_nonPGD_DCD_DEGs.csv")
write.csv(top_dbd_vs_nonpgd, "PGD_DBD_vs_nonPGD_DBD_DEGs.csv")

# To macth the columns - top_dcd_vs_dbd
top_dcd_vs_dbdwithannot <- merge(DCD_vs_DBD_DEGs, mirna_expression_long_annotated,
                           by.x = "X", by.y = "Feature.ID" )

write.csv(top_dcd_vs_dbdwithannot, "top_dcd_vs_dbdwithannot.csv")



# To macth the columns - top_dcd_vs_nonpgd
PGD_DCD_vs_nonPGD_DCD_DEGswithannot <- merge(PGD_DCD_vs_nonPGD_DCD_DEGs, mirna_expression_long_annotated,
                                 by.x = "X", by.y = "Feature.ID" )

write.csv(PGD_DCD_vs_nonPGD_DCD_DEGswithannot, "PGD_DCD_vs_nonPGD_DCD_DEGswithannot.csv")


# To macth the columns - top_dbd_vs_nonpgd
PGD_DBD_vs_nonPGD_DBD_DEGswithannot <- merge(PGD_DBD_vs_nonPGD_DBD_DEGs, mirna_expression_long_annotated,
                                             by.x = "X", by.y = "Feature.ID" )

write.csv(PGD_DBD_vs_nonPGD_DBD_DEGswithannot, "PGD_DBD_vs_nonPGD_DBD_DEGswithannot.csv")



