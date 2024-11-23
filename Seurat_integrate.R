library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

ctr=Read10X("data/ctrl/filtered_feature_bc_matrix")
Ctrl = CreateSeuratObject(counts=ctr,min.cells = 3,min.features = 200,project="ctrl")
Ctrl$dataset <- "Ctrl"
Ctrl[["percent.mt"]] <- PercentageFeatureSet(Ctrl, pattern = "^mt-")
Ctrl <- subset(Ctrl, subset = ((percent.mt<10))&(nCount_RNA>1000)&(nCount_RNA<15000))

trt = Read10X("data/treat/filtered_feature_bc_matrix")
treat = CreateSeuratObject(counts=trt,min.cells = 3,min.features = 200,project="treat")
treat$dataset <- "treat"
treat[["percent.mt"]] <- PercentageFeatureSet(treat, pattern = "^mt-")
treat <- subset(treat, subset = ((percent.mt<10))&(nCount_RNA>1000)&(nCount_RNA<15000))

Obj <- merge(Ctrl,treat)
Obj <- NormalizeData(Obj) |> 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

obj <- ScaleData(Obj)
obj <- RunPCA(obj)

obj <- IntegrateLayers(object = obj, 
                       method = CCAIntegration,
                       orig.reduction = "pca",
                       new.reduction = "integrated.cca",
                       verbose = F)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- split(obj,f = obj$dataset)

integrate <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:28)
integrate <- FindClusters(integrate, resolution = 0.2)

integrate <- RunTSNE(integrate, dims = 1:28, reduction = "integrated.cca")

VlnPlot(integrate,features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = 0)

ggsave("integrate/vlnplot.pdf",width = 12,height=12)

DimPlot(integrate, reduction = "tsne", group.by = "seurat_clusters",pt.size = 2,label=T,label.size=10)
ggsave("integrate/tsne.pdf",width = 12,height = 12)

markers = c("Epcam","Krt19","Myh11","Des","Rgs5","Pdgfrb","Pecam1","Lyve1","S100b","Gfap","Cd52","Ptprc","Acta2","Dpt","Col1a2","Col6a2")

DotPlot(integrate,features = markers,dot.scale =8)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1,vjust=1,angle = 45))+
  labs(x=NULL,y=NULL)

ggsave("integrate/markers_integrate.pdf",width = 8,height = 6)

current.cluster.ids <- 0:13
new.cluster.type <- c("Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Lymphatic","Fibroblasts","Endothelial","Fibroblasts","Immune","Pericytes","Glial","Fibroblasts")
integrate@active.ident <- plyr::mapvalues(integrate@active.ident, from = current.cluster.ids, to = new.cluster.type)
integrate@meta.data$clusterType <- integrate@active.ident

DimPlot(integrate, reduction = "tsne", group.by = "clusterType",pt.size = 3)
ggsave("integrate/tsne_cluster_type.pdf",width = 16,height = 16)

markers = c("Epcam","Rgs5","Pecam1","Lyve1","S100b","Gfap","Cd52","Ptprc","Col1a2","Pdpn")

for(i in markers){
FeaturePlot(integrate, features = i, reduction = "tsne", pt.size=3)
ggsave(paste0("integrate/featureplot_",i,".pdf"),width = 16,height = 16)
}

markers = c("Epcam","Krt19","Myh11","Des","Rgs5","Pdgfrb","Pecam1","Lyve1","S100b","Gfap","Cd52","Ptprc","Acta2","Dpt","Col1a2","Col6a2")
DotPlot(integrate,features = markers,dot.scale =8)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1,vjust=1,angle = 45))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+ 
  scale_color_gradientn(values = seq(0,1,0.2),
  colours = c("#03071e","#370617","#6a040f","#9d0208","#d00000","#dc2f02","#e85d04","#f48c06","#faa307","#ffba08"))
ggsave("integrate/markers_type2.pdf",width = 8,height = 3)

fib = subset(integrate,clusterType=="Fibroblasts")

fib <- FindNeighbors(fib, reduction = "integrated.cca", dims = 1:30)
fib <- FindClusters(fib, resolution = 0.1)

DimPlot(fib, reduction = "tsne", group.by = "seurat_clusters",pt.size = 2)

ggsave("integrate/tsne_fib.pdf",width = 10,height = 10)

current.cluster.ids <- 0:6

new.cluster.ids <- c("Fibroblasts 1","Fibroblasts 2","Fibroblasts 3","Fibroblasts 4","Fibroblasts 5","Myofibroblasts","Fibroblasts 6")

fib@active.ident <- plyr::mapvalues(fib@active.ident, from = current.cluster.ids, to = new.cluster.ids)

fib@meta.data$clusterType <- fib@active.ident

DimPlot(fib, reduction = "tsne", group.by = "clusterType",pt.size = 2)

ggsave("integrate/tsne_fib_type.pdf",width = 10,height = 10)

fib[["RNA"]] <- JoinLayers(fib[["RNA"]])
markers <- FindAllMarkers(fib,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

for(i in top10$gene){
FeaturePlot(integrate, features = i, reduction = "tsne", pt.size=3)
ggsave(paste0("integrate/featureplot_",i,".pdf"),width = 16,height = 16)
}

DimPlot(fib, reduction="tsne",label=F,pt.size = 1,split.by = "dataset")
ggsave("integrate/dataset_split_fib.pdf",width = 14,height = 7)

DoHeatmap(fib,features=top10$gene)
ggsave("integrate/heatmap_fib.pdf",width=10,height=8)

meta = fib@meta.data
meta$clusters <- as.numeric(meta$seurat_clusters)
meta_ctr = meta[meta$dataset=="Ctrl",]
meta_trt = meta[meta$dataset=="treat",]

num_ctr <- meta_ctr %>% group_by(clusters) %>% tally()
colnames(num_ctr) <- c("cluster","ctr")

num_trt <- meta_trt %>% group_by(clusters) %>% tally()
colnames(num_trt) <- c("cluster","trt")

num <- merge(num_ctr,num_trt, by="cluster")
num$ctr_percentage = num_ctr$ctr*100/colSums(num[,2:3])[[1]]
num$trt_percentage = num_trt$trt*100/colSums(num[,2:3])[[2]]

num$cluster <- factor(num$cluster,level=x_name)

ggplot(num, aes(cluster, ctr_percentage, fill = cluster)) +
  geom_bar(stat = "identity")
ggsave("integrate/cell_number_percentage_ctr.pdf",width=12,height=8)

ggplot(num, aes(cluster, trt_percentage, fill = cluster)) +
  geom_bar(stat = "identity")
ggsave("integrate/cell_number_percentage_trt.pdf",width=12,height=8)



