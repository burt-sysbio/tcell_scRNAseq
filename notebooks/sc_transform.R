require(Seurat)
require(SeuratDisk)
require(here)

path <- here("output", "cleaned_datasets")

#Convert(here("output", "cleaned_datasets", "anndata_filtered_ciucci_2019.h5ad"), dest = "h5seurat")

study <- "xin_2018"
file <- paste0("anndata_filtered_", study, ".h5seurat")

adata <- LoadH5Seurat(here("output", "cleaned_datasets", file))
adata <- SCTransform(adata, verbose = FALSE)
df_sct <- adata[["SCT"]][[]]

savename <- paste0("variance_", study, "_seurat.csv")
write.csv(df_sct, here("output", savename))


#p1 <- ggplot(df_sct, aes(sct.residual_mean))
#p1 + geom_histogram(binwidth = 0.005) + xlim(-0.1,0.5)

#p2 <- ggplot(df_sct, aes(sct.residual_variance))
#p2 + geom_histogram(binwidth = 0.1) + xlim(0,7)

#p3 <- ggplot(df_sct, aes(log10(sct.gmean), sct.residual_variance))

#p3 + geom_point(alpha=0.3, shape=16)