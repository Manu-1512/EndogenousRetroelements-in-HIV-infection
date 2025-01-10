
pbmc_u <- subset(pbmc, idents = c("PBMC_N", "PBMC_EM", "LN_N", "LN_EM"))

pbmc_u <- FindVariableFeatures(pbmc_u, selection.method = "vst", nfeatures = 3000)
pbmc_u <- ScaleData(pbmc_u, verbose = FALSE, vars.to.regress = c("percent.mt"))
pbmc.markers <- FindAllMarkers(pbmc_u, only.pos = TRUE, min.pct = 0.45, logfc.threshold = 0.25)

ervs.mark <- pbmc.markers[grep("ERV|^LTR|^MER|^MLT|SVA|^L1|^THE1", pbmc.markers$gene),]

library(ggplot2)
library(scales)
DotPlot(pbmc_u, features=mark_fam, scale.min=50, scale.max=100) + scale_colour_gradient2(low = "white", mid = "lightblue", high = "red2") + RotatedAxis()
mark_fam
mark_fam <- mark_fam[c(2:20,1,21:40)]

png("DotPlot_PBMC_LIR_HIV_Infection.png",  width = 42.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
DotPlot(pbmc_u, features=rev(mark_fam), scale.min=50, scale.max=100) + scale_colour_gradient2(low = "slateblue1", mid = "slategray1", high = "red3") + RotatedAxis() 
dev.off()
dev.off()


dat <- read.delim("ASC_significant_matrix_TEs.tsv", row.names=1, stringsAsFactors=F)

dat <- dat[,apply(dat, 2, function(row){any(row>10)})]
dat <- log2(dat+1)

chr <- gsub(".*\\.", "", colnames(dat))
fa <- c("B_cells",rep("T_cells", 11), rep("NK",4), rep("T_cells", 3), "Myeloid_DC", rep("T_cells", 2), rep("TH_precursor", 4))


cat <- as.matrix(dat)
dim(cat)

png("Allele_specific_chromatin_TE_Immune_cells.png",  width = 22.7, height = 16.8, units = "cm", res = 300, pointsize = 12)
Heatmap(cat, name = "ATAC enrichment", col = colorRamp2(c(0,1, 2,3, 4), c("white","lightblue","purple", "violetred", "maroon")),rect_gp = gpar(col = "black", lwd = 1), column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 10),
top_annotation = HeatmapAnnotation(fa3 = fa, col = list(fa3 = c("B_cells" = "red2", "T_cells" = "blue", "NK" = "green", "Myeloid_DC" = "black", "TH_precursor" = "orange")),
                    annotation_legend_param = list(fa3 = list(at =  c("B_cells", "T_cells", "NK", "Myeloid_DC", "TH_precursor")))))
dev.off()
list.files()


de <- read.delim("DEGR_loci_Elite_vs_Control.tsv", row.names=1)
dim(de)
head(de)
tail(de)


loc_manu <- de[grep("-CHR", row.names(de)),]
dim(loc_manu)
gen_manu <- de[-grep("-CHR", row.names(de)),]
dim(gen_manu)


loc_mans <- subset(loc_manu, pvalue < 0.05)

loc_mans$fam <- gsub("-.*", "", row.names(loc_mans))


manu <- de

loc_mans <- subset(loc_manu, pvalue < 0.01)


loc_mans$fam <- gsub("-.*", "", row.names(loc_mans))
loc_df <- as.data.frame(table(loc_mans$fam))
loc_df <- loc_df[order(loc_df$Freq),]

loc_mans <- subset(loc_manu, pvalue < 0.02)
loc_mans$fam <- gsub("-.*", "", row.names(loc_mans))
loc_df <- as.data.frame(table(loc_mans$fam))
loc_df <- loc_df[order(loc_df$Freq),]


png("DEGR_locus_Elite_Control_CD4_T_cells_V2.png",  width = 14.8, height = 16.8, units = "cm", res = 600, pointsize = 12)
with(de, plot(log2FoldChange, sqrt(-log10(padj)), pch=21, main="", col="grey70", xlim=c(-10,10)))
with(subset(de, log2FoldChange > 1 & padj < 0.05), points(log2FoldChange, sqrt(-log10(pvalue)), pch=21, col="pink"))
with(subset(de, log2FoldChange < -1 & padj < 0.05), points(log2FoldChange, sqrt(-log10(pvalue)), pch=21, col="lightblue"))
with(subset(loc_mans, log2FoldChange > 1 & pvalue < 0.05), points(log2FoldChange, sqrt(-log10(pvalue)), pch=18, cex=1, col="maroon"))
with(subset(loc_mans, log2FoldChange < 1 & pvalue < 0.05), points(log2FoldChange, sqrt(-log10(pvalue)), pch=18, cex=1, col="slateblue"))
abline(h = 1, col = "black", lty = 2, lwd = 1)
dev.off()


loc_df <- subset(loc_df, Freq > 5)

loc_df <- as.data.frame(table(loc_mans$fam))

loc_df <- loc_df[order(loc_df$Freq),]

 png("Barplot_DER_loci_Elite_vs_Control_V3.png",  width = 28.7, height = 14.8, units = "cm", res = 300, pointsize = 12)
par(mar=c(6,4,4,4))
my_bar <- barplot(as.matrix(loc_df$Freq), border=T , beside=T, names.arg=loc_df$Var1, cex.names=0.2, ylim=c(0,45), las=2 , space = 0.4,
                  col=rev(topo.colors(230)))
dev.off()


can <- read.delim("Candi_TEs_with_ATAC_TSS_genes.bed", header=F)

mans_can <- mans[which(row.names(mans) %in% can$V5),]
