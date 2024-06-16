library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

dge.cluster <- FindMarkers(sub_pbmc02,ident.1 = c(0:7,11:13,15,14,16),ident.2 = c(8,9,10,17))

deg<-dge.cluster

colnames(deg)
logFC=0
P.Value = 0.05
type1 = (deg$p_val< P.Value)&(deg$avg_log2FC< -logFC)
type2 = (deg$p_val< P.Value)&(deg$avg_log2FC> -logFC)
deg$Group2 = ifelse(type1,"Down",ifelse(type2,"Up","Not-Sig"))
table(deg$Group2)


diff.genes<-rownames(subset(deg,Group2!='Not-Sig'))


diff.df <- bitr(diff.genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
?bitr
go.diff2 <- enrichGO(gene = diff.df$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     pAdjustMethod = 'BH',
                     pvalueCutoff =0.01,
                     qvalueCutoff = 0.05,
                     ont="all",
                     readable =T)

deg$SYMBOL<-diff.genes
Enrichment<-deg%>%
  mutate(SYMBOL=rownames(.))%>%
  inner_join(diff.df, by = "SYMBOL")%>%
  dplyr::select(SYMBOL,ENTREZID,avg_log2FC)

go <- data.frame(Category = go.diff2$ONTOLOGY,
                 ID = go.diff2$ID,
                 Term = go.diff2$Description,
                 Genes = gsub("/", ", ", go.diff2$geneID),
                 adj_pval = go.diff2$p.adjust)
write.csv(go,'pso_GO pathway.csv')
write.csv(go,'ad_GO pathway.csv')


genelist <- data.frame(ID = Enrichment$SYMBOL,
                       logFC = Enrichment$avg_log2FC)

circ <- circle_dat(go, genelist)
write.csv(circ,file = "pso_GO.csv")
write.csv(circ,file = "ad_GO.csv")


head(circ)

reduced_circ <- reduce_overlap(circ, overlap = 0.60)
GOBubble(reduced_circ, labels = 0.5)

termNum = 10                           
geneNum = nrow(genelist)                      

chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[c(19,34,196)])#选择感兴趣的gene和Term

#AD
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[c(27,11,89)])#选择感兴趣的gene和Term


GOChord(chord)

outFile="Rplot.pdf"      
pdf(file=outFile,width=6,height=7)
GOChord(chord)
dev.off()

res<-deg
library(org.Hs.eg.db)
entrezid_all = mapIds(x = org.Hs.eg.db,  
                      keys = rownames(res), 
                      keytype = "SYMBOL", 
                      column = "ENTREZID") 
entrezid_all = na.omit(entrezid_all)  
entrezid_all = data.frame(entrezid_all) 
head(entrezid_all)
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], 
                         keyType = "kegg",
                         organism= "human",  
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = 1,  
                         qvalueCutoff =1)  
KEGG_enrich  = data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'pso_KEGG_enrich.csv')
write.csv(KEGG_enrich,'ad_KEGG_enrich.csv')


display_number = 30置
KEGG_enrich = as.data.frame(KEGG_enrich)[1:display_number[1], ]
KEGG_enrich = as.data.frame(KEGG_enrich)[c(3,19,22,24,44,46,55,56), ]
KEGG_enrich = as.data.frame(KEGG_enrich)[c(13,14,30,34,36,38,39,53), ]

kk = as.data.frame(KEGG_enrich)
rownames(kk) = 1:nrow(kk)
kk$order=factor(rev(as.integer(rownames(kk))),labels = rev(kk$Description))
#柱状图#
p1<-ggplot(kk,aes(y=order,x=Count,fill=pvalue))+
  geom_bar(stat = "identity",width=0.8)+ 
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "KEGG Pathways Enrichment",  
       x = "Gene number",
       y = "Pathway")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()

outFile="Rplot.pdf"      
pdf(file=outFile,width=6,height=3)
p1
dev.off()

library(enrichplot)

df<-dge.cluster
df$SYMBOL<-rownames(df)
df_id<-bitr(rownames(df), 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = "org.Hs.eg.db")
df_all<-merge(df,df_id,by="SYMBOL",all=F)
head(df_all) 
dim(df_all) 

df_all_sort <- df_all[order(df_all$avg_log2FC, decreasing = T),]
gene_fc = df_all_sort$avg_log2FC 
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID 
head(gene_fc)
KEGG <- gseKEGG(gene_fc, organism = "human") 
write.csv(KEGG,'pso_KEGG_GSEA_ALL.csv')
write.csv(KEGG,'ad_KEGG_GSEA_ALL.csv')


p2<-gseaplot2(KEGG, "hsa04657", color = "firebrick", rel_heights=c(1, .2, .6))
sortKEGG<-KEGG[order(KEGG$enrichmentScore, decreasing = T),]

outFile="Rplot.pdf"      
pdf(file=outFile,width=6,height=4)
p2
dev.off()