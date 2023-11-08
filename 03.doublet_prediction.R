#Doublet prediction
library(scDblFinder)
library(scater)
dblOutdir <- 'dblFinder/'
#First merge un-stable clusters C12 and C13 to C0
Ident_final = 'Cluster_final'
SeuratObject@meta.data[, Ident_final] = as.character(SeuratObject@meta.data[, Ident1])
SeuratObject@meta.data[SeuratObject@meta.data[, Ident_final] %in% c('12', '13'), Ident_final] = '0'
SeuratObject@meta.data[, Ident_final] <- SeuratClusterLableResort(SeuratObject@meta.data[, Ident_final])
SeuratObject@meta.data[, 'seurat_clusters'] = Idents(SeuratObject)
Idents(SeuratObject) <- Ident_final
Ident1 <- Ident_final
#
set.seed(10010101)
dbl_classifi <- NULL
for(s1 in unique(SeuratObject@meta.data$Sample)){
	cat('Processing', s1, '\n')
	sub1 <- SeuratObject[, SeuratObject@meta.data$Sample == s1]
	dbl_classifi0 <- scDblFinder(GetAssayData(sub1, slot = 'counts'), clusters = as.character(sub1@meta.data[, Ident1]))
	dbl_classifi0 <- colData(dbl_classifi0)
	dbl_classifi <- rbind(dbl_classifi, df0(Sample = s1, Barcode = rownames(dbl_classifi0), dbl_classifi0))
}
write.table(dbl_classifi, file = paste0(dblOutdir, 'dbl.classification.', Ident1, '.tsv'), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
#
SeuratObject@meta.data$scDblFinder.score <- dbl_classifi[match(SeuratObject@meta.data$Barcode, dbl_classifi$Barcode), 'scDblFinder.score']
SeuratObject@meta.data$scDblFinder.class <- dbl_classifi[match(SeuratObject@meta.data$Barcode, dbl_classifi$Barcode), 'scDblFinder.class']
p1 <- FeaturePlot(SeuratObject, reduction = "umap", features = 'scDblFinder.score', pt.size = 0.2, label = TRUE)
p2 <- DimPlot(SeuratObject, reduction = "umap", group.by = "scDblFinder.class")
print(p1 | p2)
SeuratObject@meta.data$scDblFinder.class <- factor(SeuratObject@meta.data$scDblFinder.class, levels = c('singlet', 'doublet'))
table(SeuratObject@meta.data$scDblFinder.class)
#Distribution of doublet scores for each cluster
p3 <- ggplot(SeuratObject@meta.data, aes(x = seurat_clusters, y = scDblFinder.score)) + 
		geom_jitter(aes(color = scDblFinder.class), size = 0.05) + 
		geom_violin(trim=FALSE, scale = "width", position="dodge") + theme_bw2()
print(p3)
#
p4 <- DimPlot(SeuratObject[, SeuratObject@meta.data$scDblFinder.class == 'singlet'], reduction = "umap", split.by = "Dx", label = TRUE, label.size = 3) + guides(color = guide_legend(override.aes = list(size = 2), ncol = 2) )
print(p4)
#
#
#**********************************************
# Create a clean version
SeuratObject <- SeuratObject[, SeuratObject@meta.data[, 'scDblFinder.class'] == 'singlet']
save.rds(SeuratObject, paste0(outdir, 'dbl.classification.SeuratObject.clean.', Ident1, '.RDS'))
