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
VlnPlot2 <- function(pp, genes, varList, axis.text.x.size=7, axis.text.y.size=7, axis.title.x.size = 7, axis.title.y.size = 7, ncol = 1, n.breaks.y = NULL, feature.label.position = 'left', strip.text.x.angle = NULL, strip.position = 'top'){
        stopifnot(is.list(pp))
        if(is(pp, 'ggplot')){ #for stacked VlnPlot
                pp$data[, 'feature'] = genes[as.character(pp$data[, 'feature']), 'Symbol']
                pp$data[, 'feature'] = factor(pp$data[, 'feature'], levels = varList)
                if(!is.null(strip.text.x.angle)) pp$theme$strip.text.x$angle <- strip.text.x.angle
                pp$theme$axis.text.x$size <- axis.text.x.size
                pp$theme$axis.text.y$size <- axis.text.y.size
                pp$theme$axis.title.x$size <- axis.title.x.size
                pp$theme$axis.title.y$size <- axis.title.y.size
                if(strip.position == 'bottom') pp$facet$params$switch <- 'x'
                return(pp)
        }
        for(i in 1:length(pp)){
                ppn=genes[pp[[i]]$labels$title, 'Symbol']
                if(is.na(ppn)){
                        cat('Try to rename', pp[[i]]$labels$title, '\n')
                        if(grepl('^rna_', pp[[i]]$labels$title, ignore.case = TRUE)) pp[[i]]$labels$title = sub('^rna_', '', pp[[i]]$labels$title, ignore.case = TRUE)
                        if(grepl('^sct_', pp[[i]]$labels$title, ignore.case = TRUE)) pp[[i]]$labels$title = sub('^sct_', '', pp[[i]]$labels$title, ignore.case = TRUE)
                        ppn = genes[pp[[i]]$labels$title, 'Symbol']
                }
                pp[[i]] <- pp[[i]] + theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(size=axis.text.x.size), 
                        axis.text.y=element_text(size=axis.text.y.size), axis.title.y=element_text(size=axis.title.y.size)) + 
                        guides(color='none', fill='none')
                if(!is.null(n.breaks.y)) pp[[i]] <- pp[[i]] + scale_y_continuous(n.breaks = n.breaks.y)
                if(feature.label.position == 'left'){
                        pp[[i]] <- pp[[i]] + labs(title=NULL, y=ppn)
                }else{
                        annotations <- data.frame(xpos = -Inf, ypos = Inf, annotateText = paste0('  ', ppn), hjustvar = 0, vjustvar = 1)
                        pp[[i]] <- pp[[i]] + labs(y = NULL, title = NULL) + geom_text(data = annotations, aes(x=xpos, y=ypos, hjust=hjustvar,
                                vjust=vjustvar, label=annotateText), inherit.aes = FALSE)
                }
                if(is.na(ppn)){
                        cat('Can not map symbol for', pp[[i]]$labels$title, '\n')
                        next
                }
                names(pp)[i]=ppn
                if(ppn!=varList[length(varList)]) pp[[i]] <- pp[[i]] + theme(axis.text.x = element_blank())
        }
        pp = pp[varList [varList %in% names(pp)] ]
        patchwork:::wrap_plots(pp, ncol = ncol)
}