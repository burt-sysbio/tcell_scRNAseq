require(Seurat)
require(here)
require(SeuratDisk)

path <- here("output", "cleaned_datasets")

#Convert(here("output", "cleaned_datasets", "anndata_filtered_ciucci_2019.h5ad"), dest = "h5seurat")

adata_ciucci <- LoadH5Seurat(here("output", "cleaned_datasets", "anndata_filtered_ciucci_2019.h5seurat"))
adata_ciucci <- NormalizeData(adata_ciucci)
adata_ciucci <- FindVariableFeatures(adata_ciucci, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(adata_ciucci), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(adata_ciucci)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

plot3 <-VlnPlot(adata_ciucci, features = top10, pt.size = 0)
plot3data <- plot3$data

topgenedata <- adata_ciucci.data[top10, 1:30]