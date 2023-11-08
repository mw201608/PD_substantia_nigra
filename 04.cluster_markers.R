#******************************************************************
#Plot and compute markers
#
genes1 <- c('SLC17A6', 'GAD1', "SLC6A3", "TH", "SLC18A2", 
			'DCC', 'GALNTL6', 'RIT2', 'RBFOX3', "AQP4", "C3", "MOG", 'VCAN', 'FLT1')
genes1 <- genes1[genes1 %in% genes$Symbol]
features1 <- genes$Geneid[match(genes1, genes$Symbol)]
#
p1 <- VlnPlot2(VlnPlot(SeuratObject[, SeuratObject@meta.data$Dx == 'Control'], features = features1, pt.size = 0, ncol = 1, same.y.lims = TRUE, combine = FALSE), genes, genes1)
p2 <- VlnPlot2(VlnPlot(SeuratObject[, SeuratObject@meta.data$Dx == 'PD'], features = features1, pt.size = 0, ncol = 1, same.y.lims = TRUE, combine = FALSE), genes, genes1)
print(p1 | p2)
#
logfc.threshold = log(1.2)
markers <- computeSeuratClusterMarkers(x = SeuratObject, logfc.threshold = logfc.threshold, method = 'FindMarkers', test.use = 'wilcox', grouping.var = NULL)
#
markers <- df0(genes[markers$Geneid, c('Geneid', 'Symbol')], markers[, -1])
write.table(markers, file = paste0(outdir, 'cluster.markers.', Ident1, '.tsv'), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
topn = markers %>% group_by(factor(Cluster, levels = levels(Idents(SeuratObject)))) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice(1:5) %>% ungroup()
d1 = SeuratObject[rownames(SeuratObject) %in% topn$Geneid, SubsetByGroup(Idents(SeuratObject), n = 500, seed = 12345)]
p1 = DoHeatmap3(d1, feature_df = as.data.frame(topn), 
	group.col = 'Cluster', id.col = 'Geneid', symbol.col = 'Symbol',
	assay = 'RNA', slot = 'data', size = 1, row_text_size = 4) + NoLegend()
print(p1)
#

# Control only
#
logfc.threshold = log(1.2)
markers <- computeSeuratClusterMarkers(x = SeuratObject[, SeuratObject@meta.data$Dx == 'Control'], logfc.threshold = logfc.threshold, method = 'FindMarkers', test.use = 'wilcox', grouping.var = NULL)
#
markers <- df0(genes[markers$Geneid, c('Geneid', 'Symbol')], markers[, -1])
write.table(markers, file = paste0(outdir, 'control.cluster.markers.', Ident1, '.tsv'), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
topn = markers %>% group_by(factor(Cluster, levels = levels(Idents(SeuratObject)))) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice(1:5) %>% ungroup()
d1 = SeuratObject[rownames(SeuratObject) %in% topn$Geneid, SeuratObject@meta.data$Dx == 'Control']
d1 = d1[, SubsetByGroup(Idents(d1), n = 500, seed = 12345)]
p1 = DoHeatmap3(d1, feature_df = as.data.frame(topn), 
	group.col = 'Cluster', id.col = 'Geneid', symbol.col = 'Symbol', 
	assay = 'RNA', slot = 'data', size = 1, row_text_size = 4) + NoLegend()
print(p1)
