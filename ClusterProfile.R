#sextounossodecadasemana
# episodio de hoje: PathfindR

setwd("~/Documents/cristal")

install.packages("pathfindR")
pak::pkg_install("egeulgen/pathfindR")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("pathview", version = "3.8") 
BiocManager::install("AnnotationDbi", version = "3.8") 
BiocManager::install("org.Hs.eg.db", version = "3.8")
install.packages("pathfindR")
library(pak)
pak::pkg_install("pathfindR")
library(KEGGREST)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGREST")

#AGORA VAI
library(KEGGREST)
# Obtain list of M.musculus pathways
mmu_kegg_descriptions <- keggList("pathway", "mmu")

# Turn the identifiers of KEGGREST into KEGG-style pathway identifiers
kegg_ids <- sub("path:", "", names(mmu_kegg_descriptions))


# Obtain and parse genes per each pathway
mmu_kegg_genes <- sapply(kegg_ids, function(pwid){
  pw <- keggGet(pwid)
  pw <- pw[[1]]$GENE[c(FALSE, TRUE)] # get gene symbols, not descriptions
  pw <- sub(";.+", "", pw) # discard any remaining description
  pw <- pw[grep("^[A-Za-z0-9_-]+(\\@)?$", pw)] # remove mistaken lines that cannot be gene symbols
  pw <- unique(pw) # keep unique symbols
  pw
})

## Filter list to exclude terms with 0 genes (metabolic pathways)
mmu_kegg_genes <- mmu_kegg_genes[sapply(mmu_kegg_genes, length) != 0]

## Form the custom descriptions vector
names(mmu_kegg_descriptions) <- sub("path:", "", names(mmu_kegg_descriptions))
mmu_kegg_descriptions <- sub(" - Mus musculus \\(mouse\\)", "", mmu_kegg_descriptions)
mmu_kegg_descriptions <- mmu_kegg_descriptions[names(mmu_kegg_descriptions) %in% names(mmu_kegg_genes)]

## Save both as RDS files for later use
saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")

# So para recordar - The protein-protein interaction network (PIN)

###############################################################
# Setting the convert2alias arguments to FALSE                #
# (because the alias conversion only works on H.sapiens genes)#
###############################################################

library(pathfindR)

mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")
inputdf <- read.delim("~/Desktop/socorro/sos.txt",header=TRUE)

## Downloading the STRING PIN file to tempdir
url <- "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz"
path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
download.file(url, path2file)

## read STRING pin file
mmu_string_df <- read.table(path2file, header = TRUE)

## filter using combined_score cut-off value of 800
mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]

## fix ids
mmu_string_pin <- data.frame(Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
                             Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2))
head(mmu_string_pin, 2)
library(biomaRt)
mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_id"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(mmu_string_pin)),
                   mart = mmu_ensembl)
mmu_string_pin$Interactor_A <- converted$mgi_id[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_id[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
#mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
#mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(mmu_string_pin, 2)

# remove self interactions
self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) # this will return a matrix object

mmu_string_pin <- data.frame(A = mmu_string_pin[, 1],
                             pp = "pp",
                             B = mmu_string_pin[, 2])

path2SIF <- file.path("/Users/tiagofalcon/Desktop/socorro", "musculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
path2SIF <- normalizePath(path2SIF)



output_df <- run_pathfindR(input = inputdf,
                           gene_sets = "Custom",
                           custom_genes = mmu_kegg_genes,
                           custom_pathways = mmu_kegg_descriptions,
                           human_genes= FALSE, pin_name_path = path2SIF)


# We can also create a graphical summary of the top 20 enrichment results using enrichment_chart():
enrichment_chart(final_res[1:20, ])


#role novo 7x1
#Tentando utilizar o PIN do biogrid - https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.180/
# file   38.27 MBBIOGRID-ORGANISM-3.5.180.tab.zip

## read STRING pin file
mmu_biogrid_df <- read.delim("/Users/tiagofalcon/Desktop/cristal/biogrid_mmu_pin.txt")

## filter using combined_score cut-off value of 800
#mmu_biogrid_df <- mmu_biogrid_df[mmu_biogrid_df$combined_score >= 800, ]

## fix ids
mmu_biogrid_pin <- data.frame(Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
                              Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2))
head(mmu_string_pin, 2)
library(biomaRt)
mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_id"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(mmu_string_pin)),
                   mart = mmu_ensembl)
mmu_string_pin$Interactor_A <- converted$mgi_id[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_id[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
#mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
#mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(mmu_string_pin, 2)

# remove self interactions
self_intr_cond <- mmu_biogrid_pin$OFFICIAL_SYMBOL_A == mmu_biogrid_pin$OFFICIAL_SYMBOL_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_biogrid_pin <- unique(t(apply(mmu_biogrid_pin, 1, sort))) # this will return a matrix object

mmu_biogrid_pin <- data.frame(A = mmu_biogrid_pin[, 3],
                              pp = "pp",
                              B = mmu_biogrid_pin[, 4])

path2SIF <- file.path("/Documents/Luis/", "musculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
path2SIF <- normalizePath(path2SIF)



output_df <- run_pathfindR(input = inputdf,
                           gene_sets = "Custom",
                           custom_genes = mmu_kegg_genes,
                           custom_pathways = mmu_kegg_descriptions,
                           human_genes= FALSE, pin_name_path = path2SIF)


# We can also create a graphical summary of the top 20 enrichment results using enrichment_chart():
enrichment_chart(PathFindrantigoResult[1:20, ])

library(biomaRt)

mmu_ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
converted <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = Aldex2_111906_filtrado$row.names.Xald.,
                   mart = mmu_ensembl)
Aldex2_111906_filtrado$X <- converted$hgnc_symbol[match(Aldex2_111906_filtrado$row.names.Xald., converted$ensembl_gene_id)]
                                                        

stats <- read.csv("~/Desktop/cristal/16S_111906_filtrado_geneID.csv",stringsAsFactors = FALSE)
stats$GeneID <- as.character(entrez$NCBI.gene.ID[match(stats$X, entrez$Gene.stable.ID)])
stats<-stats%>% distinct(GeneID, .keep_all=TRUE)
stats<-na.omit(stats)

#Criar data frame com tabela transposta com o geneID / EntrezID e FDR 
inputCP<-t(data.frame("logFC"=stats$logFC,"GeneID" = stats$GeneID ))

#Vignett Like
geneList< -t(data.frame("logFC"=data$logFC,"GeneID" = data$GeneID ))
geneList <-(data[,2])
names(geneList) <-as.character(data[,1])
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)

#carregar as libraries tidyr dplyr
#WikiPathways enrichment 
wp2gene <- read.gmt("~/Desktop/cristal/wikipathways-20180810-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
data <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
pdf("troconojento.pdf")
barplot(data, showCategory = 15)
dev.off()
head(data)
ewp2 <- GSEA(geneList, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)

datazas <- data@result
is.data.frame(datazas)

#Enriquecimento ontológico com o cluster profile
CP_test<-enrichGO(inputCP, OrgDb = org.Hs.eg.db)
#Troca do nome de EntrezID p/ HGNC
CP_test_IDs<-setReadable(CP_test,OrgDb = org.Hs.eg.db)
#Salva tabela
write.csv(CP_test_IDs@result,"ClusterProfile_enrich_111906_16S.csv")
#Plot arvore de genes
barplot(CP_test)
typeof(CP_test)

#Dessa vez testando com o KEGG - enrichMKEGG, enrichmentKEGG
library(DOSE)
library(clusterProfiler)

data(geneList, package="DOSE")
head(geneList)
gene <- names(geneList)[abs(geneList) > 2]
head(gene)
de <- names(geneList) [1:100]
yy.df <- enrichKEGG(de, pvalueCutoff = 0.01)
head(yy.df)

satanas<- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.05,
                  verbose      = TRUE)
barplot(satanas@result, showCategory = 10)

satanas2<- setReadable(satanas, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
satanas2<-na.omit(satanas2)
satanas2<- setReadable(satanas, OrgDb = org.Hs.eg.db, keytype = "ENTREZID")
barplot(yy.df, showCategory= 10)
library(enrichplot)
