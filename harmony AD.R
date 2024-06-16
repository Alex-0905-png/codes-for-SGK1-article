library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)


dir = c('0AD1/', "0AD2/", "0AD3/", "0AD4/", "0AD5/", "0AD6/", "0AD7/")
names(dir) = c('0AD1', "0AD2", "0AD3", "0AD4", "0AD5", "0AD6", "0AD7")     

counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)


data1<-subset(scRNA,subset = nFeature_RNA > 500)


data1 <- NormalizeData(data1,normalization.method="LogNormalize", scale.factor=10000)
data1<-FindVariableFeatures(data1,selection.method = "vst",nfeatures = 2000)

dir = c("ad1/", "ad2/")
names(dir) = c("ad1", "ad2")     

counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

data2<-subset(scRNA,subset = nFeature_RNA > 500)

data2 <- NormalizeData(data2,normalization.method="LogNormalize", scale.factor=10000)
data2<-FindVariableFeatures(data2,selection.method = "vst",nfeatures = 2000)

seurat1<-data1
seurat2<-data2
seurat1$dataset <- 'dataset1'
seurat2$dataset <- 'dataset2'

seurat_combined <- merge(seurat1, y = list(seurat2), add.cell.ids = c("set1", "set2"), project = "combined")

seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined)
seurat_combined <- ScaleData(seurat_combined)

seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))

seurat_combined <- RunHarmony(seurat_combined, group.by.vars = "dataset", project.dim = FALSE)


scRNA1<-seurat_combined


ElbowPlot(scRNA1,ndims=20,reduction = "pca")
pc.num=1:20
scRNA1<-FindNeighbors(scRNA1,dims = pc.num)
scRNA1<-FindClusters(scRNA1,resolution = 1.0)
scRNA1<-BuildClusterTree(scRNA1)
PlotClusterTree(scRNA1)
scRNA1=RunUMAP(scRNA1,dims = pc.num)
embed_umap<-Embeddings(scRNA1,"umap")
DimPlot(scRNA1, reduction = "umap",label = T)