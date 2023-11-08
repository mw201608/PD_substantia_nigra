library(ggplot2)
library(Seurat)
library(harmony)
library(dplyr)
library(Matrix)
library(patchwork)
library(future)
options(future.globals.maxSize= 1024^100)
outdir <- 'output/'
#Read in pre-merged count matrix data
SeuratObject <- readRDS("PD.paper.seuratObject.RDS") #data available at https://www.synapse.org/#!Synapse:syn52911948
SeuratObject@assays$RNA@data <- SeuratObject@assays$RNA@counts #to minimize the file size, data slot is empty in the saved file so it needs to recovered.
genes <- SeuratObject@assays$RNA@meta.features
# run integration/normalization
vars.to.regress = c('Sex', 'Age', 'PMI')
SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)
SeuratObject <- ScaleData(SeuratObject, features = NULL, vars.to.regress = vars.to.regress)
#
SeuratObject <- RunPCA(SeuratObject, verbose = FALSE)
ElbowPlot(SeuratObject, ndims = 50)
sds = Stdev(object = SeuratObject, reduction = "pca")
print(cumsum(sds^2) / sum(sds^2))
dims = 1:which(cumsum(sds^2) / sum(sds^2) >= 0.9)[1]
SeuratObject = RunHarmony(SeuratObject, group.by.vars = "Sample", reduction = "pca", dims.use = dims, plot_convergence = TRUE)
p1 <- DimPlot(object = SeuratObject, reduction = "harmony", pt.size = .1, group.by = "Dx")
print(p1)
SeuratObject <- RunUMAP(SeuratObject, dims = 1:20, reduction = "harmony", verbose = FALSE)
SeuratObject <- RunTSNE(SeuratObject, dims = 1:20, reduction = "harmony", verbose = FALSE)
SeuratObject <- FindNeighbors(SeuratObject, dims = 1:20, reduction = "harmony", verbose = FALSE)
SeuratObject <- FindClusters(SeuratObject, resolution = 0.2, verbose = FALSE)
