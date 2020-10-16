#sextounossodecadasemana
# episodio de hoje: PathfindR

setwd("~/Desktop/cristal")

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
library(pathfindR)
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
inputdf <- read.delim("~/Desktop/cristal/95224_16s.txt",header=TRUE)

mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
converted_input <- getBM(attributes = c("mgi_id", "mgi_symbol"),
                   filters = "mgi_id",
                   values = unique(unlist(inputdf$Data)),
                   mart = mmu_ensembl)
inputdf$Data <- converted_input$mgi_symbol[match(inputdf$Data, converted_input$mgi_id)]
write.table(inputdf,file="95224_16S_MGIsymbols.txt")


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
converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(mmu_string_pin)),
                   mart = mmu_ensembl)
mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
#mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
#mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(mmu_string_pin, 2)

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
inputdf <- na.omit(inputdf)


output_df <- run_pathfindR(input = inputdf,
                           gene_sets = "Custom",
                           custom_genes = mmu_kegg_genes,
                           custom_pathways = mmu_kegg_descriptions,
                           adj_method = "none",
                           score_quan_thr = 0.1,
                           p_val_threshold = 1,sig_gene_thr=0.5,
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
mmu_biogrid_pin <- data.frame(Interactor_A = sub("^10090\\.", "", mmu_biogrid_df$OFFICIAL_SYMBOL_A ),
                             Interactor_B = sub("^10090\\.", "", mmu_biogrid_df$OFFICIAL_SYMBOL_B))
head(mmu_biogrid_pin, 2)
library(biomaRt)
mart <- useDataset(dataset = "mmusculus_gene_ensembl", useMart("ensembl"))

mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
converted <- getBM(attributes = c("mgi_symbol", "mgi_id"),
                   filters = "mgi_symbol",
                   values = unique(unlist(mmu_biogrid_pin)),
                   mart = mmu_ensembl)
mmu_biogrid_pin$Interactor_A <- converted$mgi_id[match(mmu_biogrid_pin$Interactor_A, converted$mgi_symbol)]
mmu_biogrid_pin$Interactor_B <- converted$mgi_id[match(mmu_biogrid_pin$Interactor_B, converted$mgi_symbol)]
#mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
#mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(converted, 2)

# remove self interactions
self_intr_cond <- mmu_biogrid_pin$Interactor_A == mmu_biogrid_pin$Interactor_B
mmu_biogrid_pin <- mmu_biogrid_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_biogrid_pin <- unique(t(apply(mmu_biogrid_pin, 1, sort))) # this will return a matrix object

mmu_biogrid_pin <- data.frame(A = mmu_biogrid_pin[, 1],
                             pp = "pp",
                             B = mmu_biogrid_pin[, 2])

path2SIF <- file.path("/Users/tiagofalcon/Desktop/socorro", "musculusPIN.sif")
write.table(mmu_biogrid_pin,
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

# remove self interactions
self_intr_cond <- mmu_biogrid_pin$OFFICIAL_SYMBOL_A  == mmu_biogrid_pin$OFFICIAL_SYMBOL_B
mmu_biogrid_pin <- mmu_biogrid_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_biogrid_pin <- unique(t(apply(mmu_biogrid_pin, 1, sort))) # this will return a matrix object

mmu_biogrid_pin <- data.frame(A = mmu_biogrid_pin[, 3],
                             pp = "pp",
                             B = mmu_biogrid_pin[, 4])


output <- run_pathfindR(input = sos,
                                gene_sets = "Custom",
                                custom_genes = mmu_kegg_genes,
                                custom_descriptions = mmu_kegg_descriptions,
                                pin_name_path = path2SIF)



na.omit(zas)

pdf("zas_bp.pdf")
enrichment_chart(zas, top_terms)
dev.off()
is.character(zas$Pathway)

###########################
# Para seres humanos sapiens e pros non sapiens tbm
#######################################################
library(pathfindR)
Corrected_signal<- active_snw_search(Corrected_signal, pin_name_path = "KEGG")
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "hgnc_symbol"),
  values=Corrected_signal$GENE,
  mart=mart)

Corrected_signal$GENE<-genes$hgnc_symbol[match(Corrected_signal$GENE, genes$ensembl_gene_id)]
head(Corrected_signal)
Corrected_signal<- na.omit(Corrected_signal)
head(Corrected_signal)


write.table(Corrected_signal, file = "AA.csv", sep = "\t")

sos<- active_snw_search(sos, pin_name_path = "KEGG", search_method = "GR")
is.data.frame(sos)

RA_processed <- input_processing(input = sos, 
                                 p_val_threshold = 0.05, 
                                 pin_name_path  = "Biogrid",
                                 convert2alias = FALSE)

