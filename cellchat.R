library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(cowplot)
library(tidyverse)
library(CellChat)

af<-scRNA1
af$type <-Idents(af)
identity <- subset(af@meta.data, select = "celltype")
cellchat = createCellChat(object = af,meta =identity, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB = CellChatDB.human

showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use = CellChatDB
cellchat@netP$pathways
cellchat@DB = CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(object=cellchat,raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

mycolor <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE")

interaction_count <- cellchat@net$count
interaction_weight <- cellchat@net$weight

p1<-chordDiagram(
  x = interaction_count,
  grid.col = mycolor,                    
  directional = 1,                             
  direction.type = c("arrows", "diffHeight"),   
  diffHeight = -0.01,                            
  annotationTrack = c("name", "grid", "axis"),  
  annotationTrackHeight = c(0.05, 0.08),     
  link.arr.type = "big.arrow",                 
  link.sort = TRUE,                              
  link.largest.ontop = TRUE,                  
  transparency = 0.25                
)

p2<-chordDiagram(
  x = interaction_weight,
  grid.col = mycolor,                    
  directional = 1,                             
  direction.type = c("arrows", "diffHeight"),   
  diffHeight = -0.01,                            
  annotationTrack = c("name", "grid", "axis"),  
  annotationTrackHeight = c(0.05, 0.08),     
  link.arr.type = "big.arrow",                 
  link.sort = TRUE,                              
  link.largest.ontop = TRUE,                  
  transparency = 0.25                
)

cellchat@netP$pathways 
pathways.show <- c("TIGIT") 

table(cellchat@idents) 
vertex.receiver = seq(1,2) 
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.size = 5,layout = "hierarchy",vertex.receiver = vertex.receiver,pt.title=20,vertex.label.cex = 0.4)
