# Load DESeq2 package (install if necessary)
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

# Inspect the first few rows of the dataframe
head(genecountsd14)
#Convert gene names to character (if they aren’t already)
genecountsd14$Gene <- as.character(genecountsd14$Gene)

# Option 1: Remove rows with missing or invalid gene names
genecountsd14 <- genecountsd14[!is.na(genecountsd14$Gene) & genecountsd14$Gene != "#N/A", ]

# Option 2: Make gene names unique (even if duplicates exist)
rownames(genecountsd14) <- make.unique(genecountsd14$Gene)
# Convert gene names to character if they aren't already
genecountsd14$Gene <- as.character(genecountsd14$Gene)

# Remove rows where the gene name ends with ".1"
genecountsd14 <- genecountsd14[!grepl("\\Rik.1$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.2$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.3$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.4$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.5$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.6$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.7$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.8$", genecountsd14$Gene), ]
genecountsd14 <- genecountsd14[!grepl("\\.9$", genecountsd14$Gene), ]
# Remove rows where the gene name ends with ".1" or ends with "Rik.<digit>" (1-9)
genecountsd14 <- genecountsd14[!grepl("(\\.1$)|(Rik\\.[1-9]$)", genecountsd14$Gene), ]
rownames(genecountsd14) <- genecountsd14$Gene
# Remove duplicated gene rows (keeping the first occurrence)
genecountsd14 <- genecountsd14[!duplicated(genecountsd14$Gene), ]

# Now set the row names to the Gene column
rownames(genecountsd14) <- genecountsd14$Gene

# Optionally remove the Gene column since it's now in the row names
genecountsd14 <- genecountsd14[, -which(names(genecountsd14) == "Gene")]


# Set the gene names as row names
rownames(genecountsd14) <- genecountsd14$Gene



# Set gene names as row names and remove the gene column from the count matrix
rownames(genecountsd14) <- genecountsd14$Gene
countData <- genecountsd14

# Create sample metadata (colData)
# Here we assume sample names include "WT" or "CKO" to indicate their group.
sampleNames <- colnames(genecountsd14)
condition <- ifelse(grepl("^WT", sampleNames), "WT", "CKO")
colData <- data.frame(row.names = sampleNames,
                      condition = factor(condition))

# Create the DESeqDataSet object using the count matrix and metadata.
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# Optional: Pre-filter genes with very low counts to improve analysis speed.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run the DESeq2 pipeline to perform differential expression analysis.
dds <- DESeq(dds)

# Extract results comparing CKO vs WT (by default, the second level is compared to the first level)
res <- results(dds)

# Order the results by p-value (or adjusted p-value)
resOrdered <- res[order(res$padj), ]
summary(resOrdered)

# Optionally, write the results to a CSV file for further examination.
write.csv(as.data.frame(resOrdered), file = "DESeq2_resultsGSE264169.csv")

# Plotting an MA plot for visualization
plotMA(res, main="DESeq2 MA Plot", ylim=c(-2,2))
