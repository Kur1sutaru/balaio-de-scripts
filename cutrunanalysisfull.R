#### new analysis cut and run
setwd("~/Downloads/cutrun")

# CUT&RUN / ChIP-seq Analysis Full Script

# ---------------------------
# Load Required Libraries
# ---------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # Change to appropriate TxDb if mouse
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(clusterProfiler)
library(readxl)
library(GenomicRanges)
library(ggplot2)
library(UpSetR)
library(rtracklayer)

# ---------------------------
# Preprocess BED File
# ---------------------------
bed_path <- "H3K4me3_ASL_vs_CTRL_gained_top100.bed"
peaks_df <- read.delim(bed_path, header = TRUE)

# Convert to GRanges
peaks <- GRanges(
  seqnames = peaks_df$chrom,
  ranges = IRanges(start = peaks_df$start, end = peaks_df$end),
  strand = "*"
)

# Fix seqlevels
seqlevelsStyle(peaks) <- "UCSC"
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# ---------------------------
# Annotate Peaks
# ---------------------------
peakAnno <- annotatePeak(
  peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

# Save annotation
annot_df <- as.data.frame(peakAnno)
write.csv(annot_df, "H3K4me3_ASL_vs_CTRL_gained_top100.csv", row.names = FALSE)

# ---------------------------
# Filter Peaks by DEG Gene Lists
# ---------------------------
# Load DEG Excel sheets
deg1 <- read_excel("UPDATED 20250427 gene list sig up ASL and sig down MMF.xlsx")
deg2 <- read_excel("UPDATED 20250427 gene list sig up ASL and down MMF.xlsx")
deg3 <- read_excel("UPDATED 20250427 gene list sig down ASL and sig up MMF.xlsx")
deg4 <- read_excel("UPDATED 20250427 gene list sig down ASL and up MMF.xlsx")

# Extract and label gene lists
deg1_genes <- toupper(as.character(deg1[[1]]))
deg2_genes <- toupper(as.character(deg2[[1]]))
deg3_genes <- toupper(as.character(deg3[[1]]))
deg4_genes <- toupper(as.character(deg4[[1]]))

# Create a named list of all gene sets
gene_lists <- list(
  "sig_up_ASL_and_sig_down_MMF" = deg1_genes,
  "sig_up_ASL_and_down_MMF" = deg2_genes,
  "sig_down_ASL_and_sig_up_MMF" = deg3_genes,
  "sig_down_ASL_and_up_MMF" = deg4_genes
)

# Filter and export peaks per list
for (name in names(gene_lists)) {
  glist <- gene_lists[[name]]
  filtered <- annot_df[annot_df$SYMBOL %in% glist, ]
  write.csv(filtered, paste0("filtered_", name, ".csv"), row.names = FALSE)
}

# ---------------------------
# Genomic Feature Distribution Plots
# ---------------------------
png("H3K4me3_ASL_vs_CTRL_gained_top100annotation_pie.png", width=1200, height=800, res=150)
plotAnnoPie(peakAnno)
dev.off()

png("H3K4me3_ASL_vs_CTRL_gained_top100annotation_barplot.png", width=1200, height=800, res=150)
plotAnnoBar(peakAnno)
dev.off()

png("H3K4me3_ASL_vs_CTRL_gained_top100tss_distribution.png", width=1200, height=800, res=150)
plotDistToTSS(peakAnno, title="Peak Distribution Around TSS")
dev.off()

# ---------------------------
# TSS Enrichment Heatmap
# ---------------------------
promoters <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peaks, windows = promoters)

png("H3K4me3_ASL_vs_CTRL_gained_top100tss_enrichment_heatmap.png", width=1200, height=800, res=150)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95)
dev.off()

# ---------------------------
# Functional Enrichment (GO & KEGG)
# ---------------------------
for (name in names(gene_lists)) {
  genes <- gene_lists[[name]]
  entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # GO
  go <- enrichGO(
    gene = entrez_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  write.csv(as.data.frame(go), paste0("GO_enrichment_", name, ".csv"), row.names = FALSE)
  
  # KEGG
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # Convert gene SYMBOLs to ENTREZ IDs
  entrez_ids <- bitr(
    filtered_peaks$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Run KEGG only if valid IDs are found
  if (nrow(entrez_ids) > 0) {
    kegg <- enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = 'hsa',  # Human
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
    write.csv(as.data.frame(kegg), "KEGG_enrichment_results.csv", row.names = FALSE)
    barplot(kegg, showCategory = 10, title = "KEGG Pathway Enrichment")
  } else {
    message("⚠️ No valid Entrez IDs found for KEGG enrichment.")
  }
  
}

# ---------------------------
# UpSet Plot for Multiple Comparisons
# ---------------------------
list_input <- fromList(gene_lists)
png("H3K4me3_ASL_vs_CTRL_gained_top100upset_plot.png", width=1600, height=1000, res=150)
upset(list_input, order.by="freq")
dev.off()

# ---------------------------
# Export to IGV-compatible BED
# ---------------------------
export(peaks, "H3K4me3_ASL_vs_CTRL_gained_top100_IGV.bed", format = "BED")

# ---------------------------
# Load required libraries
# ---------------------------
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  # use TxDb.Mmusculus.UCSC.mm10.knownGene for mouse
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggplot2)


# ---------------------------
# Ensure chromosome naming compatibility
# ---------------------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(peaks) <- seqlevelsStyle(txdb)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# ---------------------------
# Get promoter regions (±3kb from TSS)
# ---------------------------
promoters <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# ---------------------------
# Generate TagMatrix
# ---------------------------
tagMatrix <- getTagMatrix(peaks, windows = promoters)

# ---------------------------
# Plot 1: Average profile of signal around TSS
# ---------------------------
plotAvgProf(
  tagMatrix,
  xlim = c(-3000, 3000),
  conf = 0.95,
  resample = 500,
  xlab = "Distance to TSS (bp)",
  ylab = "Read Density"
)

# ---------------------------
# Plot 2: TagMatrix heatmap
# ---------------------------
tagHeatmap(tagMatrix)

# Just run the plotting functions directly without dev.off()
png("H3K4me3_ASL_vs_CTRL_gained_top100.bedavg_profile_plot.png", width = 800, height = 600)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95)
dev.off()

png("H3K4me3_ASL_vs_CTRL_gained_top100.bedtag_heatmap.png", width = 600, height = 800)
tagHeatmap(tagMatrix)
dev.off()


tagMatrix <- getTagMatrix(peaks, windows = promoters)

# Optional: check matrix content
summary(as.vector(tagMatrix))

# ✅ End of Analysis



