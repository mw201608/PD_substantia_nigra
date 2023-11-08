#05. compute DEGs
library(MAST)
# MAST test
SeuratObject@meta.data$Sex2 = ifelse(SeuratObject@meta.data$Sex == 'F', 0 , 1)
latent.vars = c('Age', 'PMI', 'Sex2')
logfc.threshold <- log2(1.2)
result <- list()
for(c1 in levels(Idents(SeuratObject))){
	r1 <- try(FindMarkers(object = SeuratObject, slot = 'data', assay = 'RNA', ident.1 = 'PD', ident.2 = 'Control', group.by = 'Dx', subset.ident = c1, logfc.threshold = logfc.threshold, test.use = 'MAST', latent.vars = latent.vars), silent = TRUE)
	if (inherits(r1, "try-error")){
	    print(r1)
	    next
	}
	if (nrow(r1) == 0) next 
	result[[c1]] <- data.frame(Geneid = rownames(r1), Symbol = NA, Cluster = c1, Contrast = 'PD_vs_Control', r1, stringsAsFactors = FALSE)
}
result <- do.call(rbind, result)
rownames(result) <- NULL
result$Symbol <- genes[result$Geneid, "Symbol"]
write.table(result, file = paste0(outdir, "cluster_specific.DEGs.MAST.", Ident1, ".tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
