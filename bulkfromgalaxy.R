### bulk rna seq analysis using tabular output from galaxy ###
### workflow https://usegalaxy.org/published/workflow?id=526b05cbffecbc85 ###
setwd()

# Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(edgeR)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)


# Read tabular files
control <- read.delim("1-control-no-dox_001.fastq.tabular", header=TRUE, sep="\t")
wnt1 <- read.delim("4-Wnt1-no-dox_001.fastq.tabular", header=TRUE, sep="\t")

# Inspect the data
head(control)
head(wnt1)


# Extract gene names from target_id
control$Gene_Name <- sapply(strsplit(control$target_id, "\\|"), function(x) x[6])

# Keep only relevant columns (Gene_Name and est_counts)
control <- control %>% select(Gene_Name, est_counts)

# Rename columns
colnames(control) <- c("Gene_ID", "Control")

# View cleaned data
head(control)

# Extract gene names from target_id
wnt1$Gene_Name <- sapply(strsplit(wnt1$target_id, "\\|"), function(x) x[6])

# Keep only relevant columns (Gene_Name and est_counts)
wnt1 <- wnt1 %>% select(Gene_Name, est_counts)

# Rename columns
colnames(wnt1) <- c("Gene_ID", "Wnt1")

# View cleaned data
head(wnt1)
write.csv(control, "cleancountscontrol.csv")
write.csv(wnt1, "cleancountswnt1.csv")

# Merge into a count matrix
counts <- merge(control, wnt1, by="Gene_ID")

# Group by Gene_ID and sum duplicate entries
counts <- counts %>%
  group_by(Gene_ID) %>%
  summarise_all(sum)


write.csv(counts, "counts.csv")
# Convert tibble to dataframe
counts <- as.data.frame(counts)

# Set Gene_ID as row names
rownames(counts) <- counts$Gene_ID  

# Remove the Gene_ID column
counts <- counts[, -1]  

# Convert all count values to numeric and round them
counts <- round(as.matrix(counts))

# Check that all values are integers
str(counts)# Convert all count values to numeric and round them
counts <- round(as.matrix(counts))

# Check that all values are integers
str(counts)



# Create metadata file
metadata <- data.frame(
  row.names = colnames(counts),
  condition = c("Control", "Wnt1")  # Define sample conditions
)

metadata$condition <- factor(metadata$condition, levels = c("Control", "Wnt1"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)

# Order results by significance
res <- res[order(res$padj), ]

# Save full results
write.csv(as.data.frame(res), file = "DEG_results.csv")

### As we dont have replicates for Wnt1, we will use edgeR for DEGs
library(edgeR)

# Define the group factor
group <- factor(c("Control", "Wnt1"))

# Create the DGEList object
dge <- DGEList(counts=counts, group=group)

# Manually set dispersion
dge$common.dispersion <- 0.1  # Assumed dispersion

# Create a model matrix
design <- model.matrix(~group)

# Fit the model using quasi-likelihood
fit <- glmFit(dge, design)

# Perform likelihood ratio test
lrt <- glmLRT(fit)

# Extract results
res <- topTags(lrt, n=Inf)
write.csv(res, file="edgeR_DEG_results.csv")

# View results
head(res)

# Rename the "X" column to "Gene_name"
colnames(res)[colnames(res) == "X"] <- "Gene_name"

# Save results to CSV
write.csv(res, file = "DEG_results.csv", row.names = FALSE)

# Load the CSV file
res_fixed <- read.csv("DEG_results.csv", header = TRUE)

# Rename the first column
colnames(res_fixed)[1] <- "Gene_name"

# Save the corrected file
write.csv(res_fixed, file = "DEG_results_fixed.csv", row.names = FALSE)




# Volcano Plot
ggplot(DEG_results_fixed, aes(x=logFC, y=-log10(PValue))) +
  geom_point(alpha=0.4) +
  geom_vline(xintercept = c(-1,1), col="red", linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), col="blue", linetype="dashed") +
  theme_minimal() + 
  ggtitle("Volcano Plot of DEGs")

#PLOT HORROROSO! Vamos de enhanced volcano!
# Create volcano plot
EnhancedVolcano(DEG_results_fixed,
                lab = DEG_results_fixed$Gene_name,        # Label gene names
                x = 'logFC',       # X-axis: Log2 Fold Change
                y = 'PValue',                 # Y-axis: Adjusted p-value
                xlab = bquote(~Log[2]~ "Fold Change"),  # X-axis label
                ylab = bquote(~-Log[10]~ "PValue"),  # Y-axis label
                title = "Volcano Plot of Differential Expression",
                pCutoff = 0.05,             # Adjusted p-value threshold
                FCcutoff = 1.0,             # Log2 Fold Change threshold
                pointSize = 2.0,            # Size of points
                labSize = 4.0,              # Size of labels
                colAlpha = 0.8,             # Transparency
                legendLabels = c("NS", "Log2 FC", "P-value", "Both"),
                legendPosition = "right",
                legendLabSize = 10,
                legendIconSize = 3.0,
                col = c("black", "blue", "red", "green")  # Color coding
)




