### pathway analysis ###
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)  # Use "org.Hs.eg.db" for human genes
library(pathview)
library(biomaRt)
library(pathfindR)

# Load DEG results
deg_results <- read.csv("DEG_results_fixed.csv", header = TRUE)

# Connect to Ensembl for mouse genes
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# Convert SYMBOL to ENTREZID
gene_conversion <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                         filters = "external_gene_name",
                         values = deg_results$Gene_name, 
                         mart = mart)

# Rename columns
colnames(gene_conversion) <- c("Gene_name", "ENTREZID")

# Merge with DEG results
deg_kegg <- merge(deg_results, gene_conversion, by="Gene_name")

# Remove NA values
deg_kegg <- na.omit(deg_kegg)

# Check results
head(deg_kegg)

## if biomart is offline, try annotation dbi
# Convert SYMBOL to ENTREZID
gene_conversion <- select(org.Mm.eg.db, 
                          keys = deg_results$Gene_name, 
                          keytype = "SYMBOL", 
                          columns = "ENTREZID")

# Merge with DEG results
deg_kegg <- merge(deg_results, gene_conversion, by.x="Gene_name", by.y="SYMBOL")

# Remove NA values
deg_kegg <- na.omit(deg_kegg)

# Check results
head(deg_kegg)

# Create a named vector for KEGG analysis
gene_fc <- deg_kegg$logFC
names(gene_fc) <- deg_kegg$ENTREZID

# Sort gene list by log2 Fold Change
gene_fc <- sort(gene_fc, decreasing = TRUE)

# Run KEGG enrichment analysis
kegg_enrichment <- enrichKEGG(
  gene = names(gene_fc),   # Use Entrez IDs
  organism = "mmu",        # "mmu" for mouse, "hsa" for human
  keyType = "kegg",
  pvalueCutoff = 0.05
)

# View top enriched pathways
head(kegg_enrichment)

# Select a KEGG pathway to visualize (use pathway ID from kegg_enrichment)
pathway_id <- "mmu04917"  # Example: Cell cycle pathway

# Generate KEGG pathway map
pathview(
  gene.data = gene_fc,   # Log2 Fold Change values
  pathway.id = pathway_id,
  species = "mmu",  # "hsa" for human, "mmu" for mouse
  out.suffix = "KEGG_Map",
  kegg.native = TRUE
)

write.csv(kegg_enrichment, "KEGG_enrichment_results.csv", row.names = FALSE)


