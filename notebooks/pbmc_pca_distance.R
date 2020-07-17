
require(Seurat)
require(SeuratData)
require(dplyr)
require(here)
require(pheatmap)

data("pbmc3k")


# preprocesssing data using code from seurat vignette
pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")

# subset based on mito content, gene counts and cell
pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# normalize
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc3k), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc3k)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


all.genes <- rownames(pbmc3k)
pbmc3k <- ScaleData(pbmc3k, features = all.genes)
pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))

# subset only T cells
tcells <- c("Memory CD4 T", "Naive CD4 T", "CD8 T")
pbmc3k <- subset(pbmc3k, subset = seurat_annotations %in% tcells)


pca_scores <- as.matrix(Embeddings(pbmc3k, reduction = "pca"))
cell_dist <- dist(pca_scores, diag = F, upper = F)
cell_dist <- as.matrix(cell_dist)

cells_trunc <- cell_dist[1:200, 1:200]
distplot <- pheatmap(cells_trunc, show_rownames = F, show_colnames = F, color = inferno(255),
                     filename = "cell_distances_pca.pdf")
