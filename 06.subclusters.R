#06. subclustering
for(Clust1 in c('6', '7')){
    so1 = SeuratObject[, SeuratObject@meta.data[, Ident1] == Clust1]
    so1 <- RunPCA(so1, verbose = FALSE)
    sds = Stdev(object = so1, reduction = "pca")
    print(cumsum(sds^2) / sum(sds^2))
    dims = 1:which(cumsum(sds^2) / sum(sds^2) >= 0.9)[1]
    dims = 1:10
    so1 <- RunUMAP(so1, dims = dims, verbose = FALSE)
    so1 <- FindNeighbors(so1, dims = dims, verbose = FALSE)
    so1 <- FindClusters(so1, resolution=0.1, verbose = FALSE)
    saveRDS(so1, file = paste0(outdir, "subclusters.c", Clust1, ".RDS"))
    p1 <- DimPlot(so1, reduction = "umap", label = TRUE) + ggtitle(Ident1)
    png(paste0(outdir, "subclusters.umap.c", Clust1, ".png"), width = 3000, height = 1000, res = 300)
    print(p1)
    dev.off()
    #
    markers = computeSeuratClusterMarkers(so1, logfc.threshold = log(1.2), method = 'FindMarkers')
    write.table(markers, file = paste0(outdir, "subclusters.markers.c", Clust1, ".tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}
