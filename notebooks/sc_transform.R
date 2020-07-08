require(Seurat)
require(SeuratDisk)
require(here)
require(ggplot2)

path <- here("output", "cleaned_datasets")

#Convert(here("output", "cleaned_datasets", "anndata_filtered_ciucci_2019.h5ad"), dest = "h5seurat")

adata_ciucci <- LoadH5Seurat(here("output", "cleaned_datasets", "anndata_filtered_ciucci_2019.h5seurat"))

adata_ciucci <- SCTransform(adata_ciucci, verbose = FALSE)

df_sct <- adata_ciucci[["SCT"]][[]]


p1 <- ggplot(df_sct, aes(sct.residual_mean))
p1 + geom_histogram(binwidth = 0.005) + xlim(-0.1,0.5)

p2 <- ggplot(df_sct, aes(sct.residual_variance))
p2 + geom_histogram(binwidth = 0.1) + xlim(0,7)

p3 <- ggplot(df_sct, aes(log10(sct.gmean), sct.residual_variance))

p3 + geom_point(alpha=0.3, shape=16)