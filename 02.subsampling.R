# Cluster stability and robustness analysis by sub-sampling
subClusterDir = 'scclust/'
source('subsampling.R')
vars.to.regress = c('Sex', 'Age', 'PMI')
so1 <- CreateSeuratObject(counts = GetAssayData(SeuratObject, slot = 'counts'), min.cells = 0, min.features = 0, names.field = 2, names.delim = '.', meta.data = SeuratObject@meta.data)
so1$RNA_snn_res.0.2 <- so1$seurat_clusters <- NULL
#may take some time to run the sub-sampling analysis
subClusts <- runSubSampleClustering(so1, split.by = 'Sample', rate = 0.8, nSamples = 100, k.param = 20, resolution = 0.2, dims = 1:20, vars.to.regress = vars.to.regress)
saveRDS(subClusts, file = paste0(subClusterDir, 'subClusts.RDS'))
#
subOrgins <- lapply(subClusts, function(x, meta, Ident1) setNames(meta[names(x), Ident1], names(x)), meta = SeuratObject@meta.data, Ident1 = Ident1)
jrcp1 <- JaccardRainCloudPlot(subOrgins, subClusts)
print(jrcp1)
#
stableClusters <- AssignStableCluster(subOrgins, subClusts, jaccard_cutoff = 0.6, method = "jaccard_percent", percent_cutoff = 0.6)
save.rds(stableClusters, file = paste0(subClusterDir, 'stableClusters.jaccard_0.6.', Ident1, '.RDS'))
#
cocp <- lapply(levels(Idents(SeuratObject)), function(x, clusters, subSampleClusters){
	cat('Compute cluster', x, '...\n')
	CoClusterProbability(x = names(clusters)[clusters == x], clusters = clusters, subSampleClusters = subSampleClusters)
}, clusters = Idents(SeuratObject), subSampleClusters = subClusts)

save.rds(cocp, file = paste0(subClusterDir, 'stableClusters.cocp.', Ident1, '.RDS'))
#
cocp_mean = lapply(cocp, function(cocp1) lapply(cocp1, function(x) rowMeans(x, na.rm = TRUE)))
cocp_mean = lapply(cocp_mean, function(cocp1) do.call(rbind, lapply(names(cocp1), function(x) df0(Cluster = x, Cell = names(cocp1[[x]]), Proportion = cocp1[[x]]))))
names(cocp_mean) = levels(Idents(SeuratObject))
cocp_mean = do.call(rbind, lapply(names(cocp_mean), function(x) df0(Target = x, cocp_mean[[x]])))
write.tsv(cocp_mean, file = paste0(subClusterDir, "CoClusterDist.", Ident1, ".tsv"))
cocp_mean = read.tsv(file = paste0(subClusterDir, "CoClusterDist.", Ident1, ".tsv"))
#
cocp_mean$Target = factor(cocp_mean$Target, levels = sort(as.integer(unique(as.character(cocp_mean$Target)))))
cocp_mean$Cluster = factor(cocp_mean$Cluster, levels = sort(as.integer(unique(as.character(cocp_mean$Cluster)))))
p1 <- ggplot(cocp_mean, aes(x = Cluster, y = Proportion)) + geom_jitter(aes(color = Cluster), size = 0.05) + 
		facet_grid(Target ~ .) + guides(color = 'none') + theme_bw2() + theme(strip.background = element_blank())
print(p1)
#
p1 <- ggplot(cocp_mean[cocp_mean$Target %in% c(12, 13), ], aes(x = Cluster, y = Proportion)) + geom_jitter(aes(color = Cluster), size = 0.05) + 
	facet_wrap(~Target, ncol = 1) + guides(color = 'none') + theme_bw2() + theme(strip.background = element_blank())
print(p1)
