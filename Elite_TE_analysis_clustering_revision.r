Elite <- readRDS("~/Dropbox/Ongoing/elite/TE_loci_Seurat_object_Elite_Nat_LIR.rds")
colnames(Elite)
Elite_nat <- Elite[,1:60]
Elite_lir <- Elite[,61:128]
Elite_lir <- NormalizeData(Elite_lir, normalization.method = "CLR")
Elite_lir <- FindVariableFeatures(Elite_lir, selection.method = "vst", nfeatures = 2000)
all.genes <- VariableFeatures(object = Elite_lir)
Elite_lir <- ScaleData(Elite_lir, features = all.genes)
Elite_lir <- RunPCA(Elite_lir, npcs = 20, features = all.genes, approx=FALSE, verbose = FALSE)
Elite_nat <- NormalizeData(Elite_nat, normalization.method = "CLR")
Elite_nat <- FindVariableFeatures(Elite_nat, selection.method = "vst", nfeatures = 2000)
all.genes <- VariableFeatures(object = Elite_nat)
Elite_nat <- ScaleData(Elite_nat, features = all.genes)
Elite_nat <- RunPCA(Elite_nat, npcs = 20, features = all.genes, approx=FALSE, verbose = FALSE)
ifnb.list <- c(Elite_nat, Elite_lir)


######### Running through CCA integration
immune.combined <- IntegrateData(anchorset = immune.anchors, )
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, k.filter=70)
immune.combined <- IntegrateData(anchorset = immune.anchors)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures=200)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, k.filter=70)
immune.combined <- IntegrateData(anchorset = immune.anchors)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures=50)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, k.filter=70)
immune.combined <- IntegrateData(anchorset = immune.anchors)
ls()
dim(Elite)
te_obj = Elite


########### ######### Running through harmony integration


Elite <- readRDS("TE_loci_Seurat_object_Elite_Nat_LIR.rds")

Elite <- NormalizeData(Elite, normalization.method = "CLR")

Elite <- FindVariableFeatures(Elite, selection.method = "vst", nfeatures = 10000)

var.genes <- VariableFeatures(object = Elite)

Elite <- ScaleData(Elite, features = var.genes)

Elite <- RunPCA(Elite, features = VariableFeatures(object = Elite))

Elite@meta.data$class <- class

ifnb <- RunHarmony(Elite, group.by.vars = "class")
 
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:20)
 
ifnb <- FindClusters(object = ifnb, resolution = 1, graph.name="RNA_snn")
 
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:10)

ifnb <- FindClusters(object = ifnb, resolution = 1.2, graph.name="RNA_snn")


png("Proportion_studies_each_cluster_only_TEs_redone_class_August24.png",  width = 8.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
barplot(prop.table(table(Idents(ifnb), ifnb$orig.ident), margin = 2)*100, col=c("maroon","gold","forestgreen","slateblue"), las=2)
dev.off()


png("Elite_UMAP_reannotated_only_TEs_reclust_class_August24.png",  width = 10.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(ifnb, reduction = "umap", cols=c("maroon", "gold", "forestgreen", "slateblue"), pt.size=3)
 dev.off()
 


te_obj_clust4 <- ifnb



Elite <- readRDS("TE_loci_Seurat_object_Elite_Nat_LIR.rds")

Elite <- NormalizeData(Elite, normalization.method = "LogNormalize", scale.factor = 100000)

Elite <- FindVariableFeatures(Elite, selection.method = "vst", nfeatures = 10000)

var.genes <- VariableFeatures(object = Elite)

Elite <- ScaleData(Elite, features = var.genes)

Elite <- RunPCA(Elite, features = VariableFeatures(object = Elite))

Elite@meta.data$class <- class


 ifnb <- RunHarmony(Elite, group.by.vars = "class")
 
 ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:20)
 
 ifnb <- FindClusters(object = ifnb, resolution = 1, graph.name="RNA_snn")

png("Elite_UMAP_reannotated_only_TEs_reclust_class_August24_log.png",  width = 10.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(ifnb, reduction = "umap", cols=c("maroon", "gold", "forestgreen", "slateblue", "black"), pt.size=3)
 dev.off()


ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:10)


png("Proportion_studies_each_cluster_only_TEs_redone_class_August24_log_V2.png",  width = 8.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
barplot(prop.table(table(Idents(ifnb), ifnb$orig.ident), margin = 2)*100, col=c("maroon","darkorange2","forestgreen","slateblue", "black"), las=2)
dev.off()


 png("Elite_UMAP_reannotated_only_TEs_reclust_class_August24_V2.png",  width = 10.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(te_obj_clust4, reduction = "umap", cols=c("maroon", "darkorange2", "forestgreen", "slateblue", "black"), pt.size=3)
 dev.off()
 

 new <- c("cluster1","cluster1","cluster2","cluster3","cluster4")

 obj <- Rename_Clusters(seurat_object = ifnb, new_idents = new, meta_col_name = "New_Idents")

names(new) <- levels(ifnb)

obj <- RenameIdents(ifnb, new)


 png("Elite_UMAP_reannotated_only_TEs_reclust_class_August24_log_V2.png",  width = 10.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
 DimPlot(obj, reduction = "umap", cols=c("maroon", "orange", "forestgreen", "slateblue"), pt.size=3)
dev.off()

 mark_4 <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5) 
 mark_sub <- subset(mark_4, p_val_adj < 0.01)


mark_genes <- row.names(mark_sub)


png("TE_mark_EC_clusters_redone_CLR.png",  width = 34.7, height = 17.8, units = "cm", res = 600, pointsize = 12)
 DoHeatmap(subset(obj, downsample = 15), features = mark_genes,draw.lines=F,  size = 1) + scale_fill_gradientn(colours=c("steelblue", "lightblue", "white","red1", "red3"), values=rescale(c( -2, -0.5, 0, 1.5, 2)), guide="colorbar")
 dev.off()


mark <-  mark_4

mark$gene <-  gsub("-.*","", gsub("_DUP.*", "", (sapply(strsplit(mark$gene,"TE-"),"[[",2))))


cluster_table <- table(mark$cluster, mark$gene)

manu <- t(cluster_table)

manu <- manu[apply(manu, 1, function(row){any(row>5)}),]

sum <- apply(manu, 1, sum)

manu.pt <- manu*100/sum

manu.pt.candi <- manu.pt[grep("HERVL74|HERVK9|HERVH|HERVL|HERVE|HERVIP10F|HERVK22|HERVK|HERVK3|LTR8|ERV3|LTR12C|THE1B|L1HS|L1PA2|L1PA4|L1PA3|L1PA5|L1PA6|LTR5|SVA|LTR49|MER20|MER41|MLT1J", row.names(manu.pt)),]


mark_genes <- mark_sub %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


mark_genes <- mark_sub %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

 png("TE_mark_EC_clusters_redone_August_24.png",  width = 34.7, height = 17.8, units = "cm", res = 600, pointsize = 12)
 DoHeatmap(subset(obj, downsample = 15), features = mark_genes$gene, draw.lines=F,  size = 4) + scale_fill_gradientn(colours=c("steelblue", "lightblue", "white","red1", "red3"), values=rescale(c( -2, -0.5, 0, 1.5, 2)), guide="colorbar")
dev.off()

mark_4_genes <- mark_genes$gene


png("TE_mark_EC_clusters_redone_August_24.png",  width = 34.7, height = 17.8, units = "cm", res = 600, pointsize = 12)
 DoHeatmap(subset(obj, downsample = 15), features = mark_4_genes, draw.lines=F,  size = 4) + scale_fill_gradientn(colours=c("steelblue", "lightblue", "white","red1", "red3"), values=rescale(c( -2, -0.5, 0, 1.5, 2)), guide="colorbar")
dev.off()



png("TE_mark_EC_clusters_redone_August_24.png",  width = 34.7, height = 17.8, units = "cm", res = 600, pointsize = 12)
DoHeatmap(subset(obj, downsample = 15), features = mark_4_genes, draw.lines=F,  size = 4) + scale_fill_gradientn(colours=c("steelblue", "lightblue", "white","red1", "red3"), values=rescale(c( -2, -0.5, 0, 1.5, 2)), guide="colorbar")
dev.off()

write.table(manu.pt.candi, "Candidate_TE_families_proportion_in_each_cluster_4.tsv", row.names=T, col.names=NA, sep="\t", dec=".", quote=F)

savehistory("Final_TE_analysis_for_Elite_paper.Rhistory")

write.table(mark_sub, "TE_loci_at_family_annotations_marking_each_cluster_4.tsv", row.names=T, col.names=NA, sep="\t", dec=".", quote=F)

save.image("Elite_TE_analysis_clustering_August2024.RData")
savehistory("Elite_TE_analysis_clustering_August2024.Rhistory")
