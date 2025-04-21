# R SCRIPT: Supervised & Unsupervised Clustering, Volcano Plot, GO/KEGG/Reactome Enrichment Analysis (Mouse)
setwd("/Users/cristalvillalba/Desktop/Lee lab bioinfo/keren/ruvseq")
# === 1. Load Required Libraries ===
library(tidyverse)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(enrichplot)

# === 2. Load Data ===
expr <- read_csv("genecountsedsxudn.csv")

# Ensure gene names are kept
gene_symbols <- expr$Gene.symbol
expr <- expr %>% select(-Gene.symbol)

# === 3. Define Sample Groups ===
eds_samples <- grep("^sample_", colnames(expr), value = TRUE)
udn_samples <- grep("^UDN", colnames(expr), value = TRUE)
group <- c(rep("EDS", length(eds_samples)), rep("UDN", length(udn_samples)))

# === 4. Supervised Clustering Heatmap (Top 100 Variable Genes) ===
expr_matrix <- as.matrix(expr[, c(eds_samples, udn_samples)])
rownames(expr_matrix) <- gene_symbols

# Z-score normalization
gene_var <- apply(expr_matrix, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:100]
z_expr <- t(scale(t(expr_matrix[top_genes, ])))

# Sample group annotation
annotation_col <- data.frame(Group = factor(group))
rownames(annotation_col) <- colnames(z_expr)

# Heatmap
pheatmap(z_expr, annotation_col = annotation_col, show_rownames = TRUE,
         show_colnames = TRUE, main = "Supervised Clustering - Top 100 Genes")

# === 5. Unsupervised Clustering Heatmap ===
pheatmap(z_expr, cluster_cols = TRUE, cluster_rows = TRUE, show_rownames = TRUE,
         main = "Unsupervised Clustering - Top 100 Genes")

# === 6. Differential Expression (limma) ===
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- levels(factor(group))
fit <- lmFit(expr_matrix[, c(eds_samples, udn_samples)], design)
contrast <- makeContrasts(UDN - EDS, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number = Inf, adjust = "fdr")
results$Gene <- gene_symbols

# === 7. Volcano Plot ===
EnhancedVolcano(results,
                lab = results$Gene,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot: UDN vs EDS',
                subtitle = 'Differential Gene Expression',
                pointSize = 3.0,
                labSize = 3.5,
                drawConnectors = TRUE,
                widthConnectors = 0.5
)

# === 8. Enrichment Analysis (GO, KEGG, Reactome) ===
# Get significant genes
sig_genes <- results %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>% pull(Gene) %>% unique()

# Mouse gene symbols use first-letter capitalization
sig_genes <- gsub("^(.)", "\\U\\1", tolower(sig_genes), perl = TRUE)

# Convert to Entrez IDs
entrez_ids <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO BP
ego <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(ego, showCategory = 15, title = "GO Biological Processes")

# KEGG
ekegg <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)
dotplot(ekegg, showCategory = 15, title = "KEGG Pathways")

# Reactome
reactome <- enrichPathway(
  gene = entrez_ids$ENTREZID,
  organism = "mouse",
  pvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(reactome, showCategory = 15, title = "Reactome Pathways")

# === 9. Save Enrichment Results ===
write.csv(as.data.frame(ego), "GO_mouse_results.csv")
write.csv(as.data.frame(ekegg), "KEGG_mouse_results.csv")
write.csv(as.data.frame(reactome), "Reactome_mouse_results.csv")