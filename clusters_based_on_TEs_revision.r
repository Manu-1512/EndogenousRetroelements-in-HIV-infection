
dat <- read.delim("Merged_data_CD4_T_Elite_Nature_Cell.tsv", header=T, row.names=1, stringsAsFactors=F)

pbmc <- CreateSeuratObject(counts = dat, project = "pbmc3k", min.cells = 3, min.features = 200)
library(scater) # load the library
library(scran) #
pbmc.sce <- as.SingleCellExperiment(pbmc)
exprMat <- counts(pbmc.sce)
pbmc <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(exprMat),sparse = T), project = "Manu")
pbmc <- NormalizeData(pbmc, normalization.method = "RC", scale.factor = 1e6)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)

 my_ids <- c(rep(c("HC_U", "HC_S"), 6), rep(c("HC_S", "HC_U"), 10), rep("HC_S", 3), "HC_U",  rep("ART", 10), rep(c("EC_CM", "EC_EM", "EC_N", "EC_TM", "EC_Total"), 12), rep(c("EC_N", "EC_CM", "EC_TM", "EC_EM"), 12), 
  "EC_LN_N", "EC_LN_EM", "EC_N", "EC_CM", "EC_TM", "EC_EM", "EC_LN_N", "EC_LN_EM", "EC_N", "EC_CM", "EC_TM", "EC_EM", "EC_LN_N", "EC_LN_EM", "EC_N", "EC_CM", "EC_TM", "EC_EM", "EC_LN_N", "EC_LN_EM")
my_ids
pbmc@meta.data$celltypes <- my_ids
pbmc@meta.data
library(harmony)

study <- c(rep("NG", 36), rep("Nat", "70"), rep("Cell", 68))
pbmc@meta.data$study <- study
 library(harmony)
ifnb <- RunHarmony(pbmc, group.by.vars = "study")
 ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:15)
 ifnb <- FindClusters(object = ifnb, resolution = 0.5, graph.name="RNA_snn")
  
  ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:15)
ifnb@meta.data
 ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30)
 ifnb <- FindClusters(object = ifnb, resolution = 1, graph.name="RNA_snn")
  ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:15)

s_obj <- pbmc
pbmc <- ifnb
Idents(pbmc) <- pbmc@meta.data$celltypes


my_cols <-c("forestgreen", "red2", "black", "orange", "gold", "royalblue2","maroon", "midnightblue", "purple", "darkcyan")

png("Include_CD4_UMAP_EC_HC_ART.png",  width = 12.8, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(pbmc,reduction="umap", cols= my_cols, pt.size=3)
dev.off()


png("ViolinPlot_key_TEs_CD4_UMAP_EC_HC_ART.png",  width = 18.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
 VlnPlot(pbmc, c("L1PA2", "LTR7", "LTR12C", "THE1B"), ncol=2, col = my_cols,  raster=F, log = TRUE, add.noise=F, same.y.lims=FALSE, pt.size=2, alpha=0.3) 
dev.off()
pbmc <- s_obj
levels(pbmc)
pbmc@meta.data
my_ids
ja <- readRDS("~/Documents/Winter22/HIV/salmon/merged/Elite_Jiang_Boritz_merged_analysis.rds")

clusters <- ja@meta.data$seurat_clusters

pbmc@meta.data$seurat_clusters[37:174] <- clusters

new_id[1:36] <- c(rep(c("HC_U", "HC_S"), 6), rep(c("HC_S", "HC_U"), 10), rep("HC_S", 3), "HC_U")
new_id
new_id[37:174] = paste('Cluster', new_id[37:174], sep='_')
new_id
levels(pbmc)
pbmc@meta.data$seurat_clusters <- Idents(pbmc)
pbmc@meta.data$clusters <- new_id
Idents(pbmc) <- new_id
levels(pbmc)

my_levels <- c("HC_U","HC_S","Cluster_0","Cluster_1","Cluster_2","Cluster_3")

Idents(pbmc) <- pbmc@meta.data$clusters
levels(pbmc) <- my_levels

my_cols <- c("forestgreen", "purple", "red2", "gold", "royalblue2", "black")

png("ViolinPlot_key_TEs_CD4_UMAP_EC_HC_ART_V2.png",  width = 18.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
 VlnPlot(pbmc, c("L1HS", "HERVH", "MER41B", "LTR8"), ncol=2, col = my_cols,  raster=F, log = FALSE, add.noise=T, same.y.lims=FALSE, pt.size=2, alpha=0.3) 
dev.off()

ifnb <- FindVariableFeatures(object = ifnb, selection.method = "vst", nfeatures = 2000)
ifnb <- ScaleData(object = ifnb, features=all.genes)
ifnb <- RunPCA(object = ifnb)
ifnb <- RunHarmony(ifnb, group.by.vars = "study")
ifnb <- FindNeighbors(object = ifnb, dims = 1:20)
ifnb <- RunUMAP(object = ifnb, dims = 1:20)
png("Elite_first__clustering_HC_MVG_V2.png",  width = 12.8, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(ifnb,reduction="umap", cols= my_cols, pt.size=3)
dev.off()
savehistory("second_shot_analysis_final_revision.Rhistory")
