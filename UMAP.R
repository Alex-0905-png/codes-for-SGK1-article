Sys.setenv(LANGUAGE="en")
options(stringsAsFactors = FALSE)
rm(list = ls())
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)

###pso
dir = c('PSO1/', "PSO2/", "PSO3/", "PSO4/", "PSO5/",
        'PSO6/', "PSO7/", "PSO8/", "PSO9/", "PSO10/",
        'PSO11/', "PSO12/", "PSO13/")
names(dir) = c('PSO1', "PSO2", "PSO3", "PSO4", "PSO5",
               'PSO6', "PSO7", "PSO8", "PSO9", "PSO10",
               'PSO11', "PSO12", "PSO13") 

###con
dir = c('CO1/', "CO2/", "CO3/", "CO4/", "CO5/")
names(dir) = c('CO1', "CO2", "CO3", "CO4", "CO5") 


counts <- Read10X(data.dir =dir)
scRNA=CreateSeuratObject(counts,min.cells = 3,project="os",min.features = 300)

seurat_combined <- NormalizeData(scRNA)
seurat_combined <- FindVariableFeatures(seurat_combined)
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))

scRNA1<-seurat_combined


ElbowPlot(scRNA1,ndims=20,reduction = "pca")
pc.num=1:20
scRNA1<-FindNeighbors(scRNA1,dims = pc.num)
scRNA1<-FindClusters(scRNA1,resolution = 1.0)
scRNA1<-BuildClusterTree(scRNA1)
PlotClusterTree(scRNA1)
scRNA1=RunUMAP(scRNA1,dims = pc.num)

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)

DimPlot(scRNA1,reduction = "umap",label = TRUE)+theme_dr(xlength = 0.22, ylength = 0.22,
                                                         arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=5,height=4)
DimPlot(scRNA1, reduction = "umap",
        pt.size = 0.5,
        label = T,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

FeaturePlot(scRNA1,features = "CD3D")+theme_dr(xlength = 0.22, ylength = 0.22,
                                               arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(scRNA1,features = "TIGIT")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

FeaturePlot(scRNA1,features = "GZMB")+theme_dr(xlength = 0.22, ylength = 0.22,
                                               arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(scRNA1,features = "GZMB")+theme_dr(xlength = 0.22, ylength = 0.22,
                                               arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

FeaturePlot(scRNA1,features = "CD68")+theme_dr(xlength = 0.22, ylength = 0.22,
                                               arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(scRNA1,features = "CD68")+theme_dr(xlength = 0.22, ylength = 0.22,
                                               arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

FeaturePlot(scRNA1,features = "KRT6A")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(scRNA1,features = "KRT6A")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

FeaturePlot(scRNA1,features = "HLA-DQB1")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                   arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(scRNA1,features = "HLA-DQB1")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                   arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27),
  celltype = c("T_cells", "Keratinocytes", "Keratinocytes","T_cells","Keratinocytes","T_cells",
               "Keratinocytes","DC","Macrophages","DC","T_cells",
               "Keratinocytes","T_cells","T_cells","DC","Keratinocytes",
               "Keratinocytes","T_cells","Keratinocytes","T_cells","Keratinocytes",
               "Keratinocytes","T_cells","Keratinocytes","Keratinocytes","Keratinocytes",
               "T_cells","Keratinocytes")
)

scRNA1@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  scRNA1$celltype[scRNA1$seurat_clusters == cluster_id] <- cell_type
}

outFile="umap.pdf"  
pdf(file=outFile,width=5,height=4)
DimPlot(scRNA1, reduction = "umap",group.by = "celltype",
        pt.size = 0.5,
        label = T,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

markers <- FindAllMarkers(object = scRNA1, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 
write.csv(markers,file = "pso_marker.csv")