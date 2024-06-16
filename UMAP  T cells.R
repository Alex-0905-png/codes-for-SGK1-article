pbmc1 <- scRNA1
sub_pbmc02_matrix<- pbmc1@assays$RNA@counts[,which(pbmc1@meta.data$seurat_clusters%in%c(0,5,3,13,19,12,17,10,26))]
sub_pbmc02 <- CreateSeuratObject(counts = sub_pbmc02_matrix,project = 'sub_02')
sub_pbmc02 <- NormalizeData(sub_pbmc02,normalization.method = "LogNormalize",scale.factor = 10000)
sub_pbmc02 <- FindVariableFeatures(sub_pbmc02,selection.method = 'vst',nfeatures = 2000)
#
sub_pbmc02 <- ScaleData(sub_pbmc02)
sub_pbmc02<-RunPCA(sub_pbmc02,features = VariableFeatures(scRNA1))
DimPlot(sub_pbmc02,reduction = "pca",group.by = "orig.ident")
ElbowPlot(sub_pbmc02,ndims=20,reduction = "pca")
sub_pbmc02 <- FindNeighbors(sub_pbmc02, dims = 1:5) 
sub_pbmc02 <- FindClusters(sub_pbmc02, resolution = 1) 
head(Idents(sub_pbmc02),8)
sub_pbmc02<- RunUMAP(sub_pbmc02, dims = 1:10) 

DimPlot(sub_pbmc02,reduction = "umap",label = TRUE)+theme_dr(xlength = 0.22, ylength = 0.22,
                                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
DimPlot(sub_pbmc02, reduction = "umap",
        pt.size = 0.5,
        label = T,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()


FeaturePlot(sub_pbmc02,features = "TIGIT")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())

outFile="umap.pdf"  
pdf(file=outFile,width=4,height=4)
FeaturePlot(sub_pbmc02,features = "TIGIT")+theme_dr(xlength = 0.22, ylength = 0.22,
                                                    arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()

cluster_types <- data.frame(
  cluster = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14),
  celltype = c("CD4", "CD4_Treg", "CD4","CD4","CD4_Treg","CD4",
               "CD8","CD4","CD8","CD4_Treg","NKT",
               "CD8","NKT","NKT","CD8")
)

sub_pbmc02@meta.data$celltype ="NA"
for(i in 1:nrow(cluster_types)) {
  cluster_id <- cluster_types$cluster[i]
  cell_type <- cluster_types$celltype[i]
  sub_pbmc02$celltype[sub_pbmc02$seurat_clusters == cluster_id] <- cell_type
}

outFile="Rplots.pdf"      
pdf(file=outFile,width=4,height=4)
DimPlot(sub_pbmc02, reduction = "umap",group.by = "celltype",
        pt.size = 0.5,
        label = F,label.box = F
)+theme_dr(xlength = 0.22, ylength = 0.22,
           arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())
dev.off()


###########
cell_counts <- sub_pbmc02@meta.data %>%
  group_by(celltype) %>%
  summarise(count = n())

print(cell_counts)



