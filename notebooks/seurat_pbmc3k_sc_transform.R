# load data sets from SeuratData to use with sc_transform

require(Seurat)
require(SeuratData)
require(dplyr)
require(here)

AvailableData()

data("pbmc3k")

# preprocesssing data using code from seurat vignette
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# subset based on mito content, gene counts and cell
pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# subset only T cells
tcells <- c("Memory CD4 T", "Naive CD4 T", "CD8 T")
pbmc3k_tcells <- subset(pbmc3k, subset = seurat_annotations %in% tcells)



# run SC transform and save output
study <- "pbmc3k_10x"
pbmc3k <- SCTransform(pbmc3k, vars.to.regress = "percent.mt", verbose = TRUE)
sct_pbmc3k <- pbmc3k[["SCT"]][[]]

savename <- paste0("variance_", study, "_seurat.csv")
write.csv(sct_pbmc3k, here("output", savename))

# run SC transform and save output
study <- "pbmc3k_10x_tcells"
pbmc3k_tcells <- SCTransform(pbmc3k_tcells, vars.to.regress = "percent.mt", verbose = TRUE)
sct_pbmc3k_cells <- pbmc3k_tcells[["SCT"]][[]]

savename <- paste0("variance_", study, "_seurat.csv")
write.csv(sct_pbmc3k_cells, here("output", savename))
