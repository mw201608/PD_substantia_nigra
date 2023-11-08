library(Seurat)
library(scclusteval)
library(harmony)
#https://github.com/crazyhottommy/scclusteval
#https://crazyhottommy.github.io/EvaluateSingleCellClustering/5k_pbmc.html
RandomSubsetData <- function(object, rate, random.subset.seed = NULL, ...){
        ncells <- nrow(object@meta.data)
        ncells.subsample <- round(ncells * rate)
        set.seed(random.subset.seed)
        selected.cells <- sample(colnames(object), ncells.subsample)
        object <- subset(object, cells =  selected.cells, ...)
        return(object)
}
runSubSampleClustering <- function(object, normalization.method = "LogNormalize", integration.method = 'harmony', rate = 0.8, nSamples = 100, k.param = 20, resolution = 0.8, dims = 1:20, split.by = NULL, vars.to.regress = NULL, seed = 1234){
	#object is a Seurat object
    #return a list of data.frame, each contains the cell id and cluster assignment.
    set.seed(seed)
	results <- lapply(1:nSamples, function(i, dims = NULL, split.by, normalization.method, integration.method, vars.to.regress){
        cat('Sub sample', i, '...\n')
		s1 = RandomSubsetData(object, rate = rate)
        s1 <- NormalizeData(s1, normalization.method = normalization.method, scale.factor = 10000)
        s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 2000)
        s1 <- ScaleData(s1, features = NULL, vars.to.regress = vars.to.regress)
		s1 <- RunPCA(s1, verbose = FALSE)
		if(is.null(dims)){
    		sds = Stdev(object = s1, reduction = "pca")
            dims = 1:which(cumsum(sds^2) / sum(sds^2) >= 0.9)[1]
		}
        if(integration.method == 'harmony'){
    		s1 = RunHarmony(s1, group.by.vars = split.by, reduction = "pca", dims.use = dims, plot_convergence = FALSE)
			s1 <- RunUMAP(s1, dims = dims, reduction = "harmony", verbose = FALSE)
			s1 <- RunTSNE(s1, dims = dims, reduction = "harmony", verbose = FALSE)
			s1 <- FindNeighbors(s1, dims = dims, reduction = "harmony", k.param = k.param, verbose = FALSE)
		}else{
			s1 <- RunUMAP(s1, dims = dims, reduction = "pca", verbose = FALSE)
			s1 <- RunTSNE(s1, dims = dims, reduction = "pca", verbose = FALSE)
			s1 <- FindNeighbors(s1, dims = dims, reduction = "pca", k.param = k.param, verbose = FALSE)
		}
		s1 <- FindClusters(s1, resolution = resolution, verbose = FALSE)
        Idents(s1)
	}, dims = dims, split.by = split.by, normalization.method = normalization.method, integration.method = integration.method, vars.to.regress = vars.to.regress)
}
CoClusterProbability <- function(x, clusters, subSampleClusters) {
    #x, a vector of cell ids in a cluster of interest
    #y, a named vector of clustering from the full data
    #subSampleClusters is a list of clusterings from repeated subsamples; each entry is a named vector (factor) of clustering result
    result = lapply(levels(clusters), function(X, x, nc){
        temp = matrix(0, nrow = length(x), ncol = nc)
        rownames(temp) = x
        temp[] = NA
        temp
    }, x = x, nc = length(subSampleClusters))
    names(result) = levels(clusters)
    i = 0
    for(s in subSampleClusters){
        i = i + 1
        cat('Process subsample', i, '...\n')
        x_s = x[x %in% names(s)]
        y_s = clusters[names(s)]
        if(length(x_s) > 0) for(lvls in levels(s)){
            cells = names(s)[s == lvls]
            cells_x = intersect(cells, x_s)
            if(length(cells_x) == 0) next
            cells_y = intersect(cells, names(y_s))
            cells_y = droplevels(y_s[cells_y])
            for(lvly in levels(cells_y)){
                z = cells_y[cells_y == lvly]
                result[[lvly]][cells_x, i] = length(z)/sum(y_s == lvly)
            }
        }
    }
    result
}
