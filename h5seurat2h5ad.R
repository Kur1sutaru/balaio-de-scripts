###convert h5seurat in h5ad

library(Seurat)
library(SeuratDisk)
library(hdf5r)

# Convert Assay5 to a standard Seurat Assay
Non_neuron_mouse_v3singlenuclei[["RNA"]] <- as(Non_neuron_mouse_v3singlenuclei[["RNA"]], "Assay")

# Ensure RNA assay is set as default
DefaultAssay(Non_neuron_mouse_v3singlenuclei) <- "RNA"

# Save the converted Seurat object
SaveH5Seurat(Non_neuron_mouse_v3singlenuclei, filename = "drg_non_neuron_sn.h5Seurat", overwrite = TRUE)

# Try converting to h5ad again
Convert("drg_non_neuron_sn.h5Seurat", dest = "h5ad", overwrite = TRUE)