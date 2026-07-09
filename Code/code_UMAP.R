###scRNA-seq
.libPaths(c("~/SeuratV4",.libPaths()))
library(Seurat)
packageVersion("Seurat")

setwd("~/Single_cell_APEX1/scRNA-seqGSE207493")

library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(dplyr)

K81.data <- Read10X(data.dir = "RCC81/")
kid81 <- CreateSeuratObject(counts = K81.data, project = "mRCC81", min.cells = 8, min.features = 500)
kid81[["percent.mt"]] <- PercentageFeatureSet(kid81, pattern = "^MT-")
VlnPlot(kid81, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid81 <- subset(kid81, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

K84.data <- Read10X(data.dir = "RCC84/")
kid84 <- CreateSeuratObject(counts = K84.data, project = "mRCC84", min.cells = 7, min.features = 500)
kid84[["percent.mt"]] <- PercentageFeatureSet(kid84, pattern = "^MT-")
VlnPlot(kid84, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid84 <- subset(kid84, subset = nFeature_RNA > 500 & nFeature_RNA < 4700 & percent.mt < 10)

K86.data <- Read10X(data.dir = "RCC86/")
kid86 <- CreateSeuratObject(counts = K86.data, project = "mRCC86", min.cells = 9, min.features = 500)
kid86[["percent.mt"]] <- PercentageFeatureSet(kid86, pattern = "^MT-")
VlnPlot(kid86, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid86 <- subset(kid86, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K87.data <- Read10X(data.dir = "RCC87/")
kid87 <- CreateSeuratObject(counts = K87.data, project = "mRCC87", min.cells = 9, min.features = 500)
kid87[["percent.mt"]] <- PercentageFeatureSet(kid87, pattern = "^MT-")
VlnPlot(kid87, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid87 <- subset(kid87, subset = nFeature_RNA > 500 & nFeature_RNA < 2900 & percent.mt < 10)

K94.data <- Read10X(data.dir = "RCC94/")
kid94 <- CreateSeuratObject(counts = K94.data, project = "mRCC94", min.cells = 8, min.features = 500)
kid94[["percent.mt"]] <- PercentageFeatureSet(kid94, pattern = "^MT-")
VlnPlot(kid94, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid94 <- subset(kid94, subset = nFeature_RNA > 500 & nFeature_RNA < 3700 & percent.mt < 10)

K96.data <- Read10X(data.dir = "RCC96/")
kid96 <- CreateSeuratObject(counts = K96.data, project = "mRCC96", min.cells = 8, min.features = 500)
kid96[["percent.mt"]] <- PercentageFeatureSet(kid96, pattern = "^MT-")
VlnPlot(kid96, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid96 <- subset(kid96, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K99.data <- Read10X(data.dir = "RCC99/")
kid99 <- CreateSeuratObject(counts = K99.data, project = "mRCC99", min.cells = 9, min.features = 500)
kid99[["percent.mt"]] <- PercentageFeatureSet(kid99, pattern = "^MT-")
VlnPlot(kid99, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid99 <- subset(kid99, subset = nFeature_RNA > 500 & nFeature_RNA < 2800 & percent.mt < 10)

K100.data <- Read10X(data.dir = "RCC100/")
kid100 <- CreateSeuratObject(counts = K100.data, project = "mRCC100", min.cells = 10, min.features = 500)
kid100[["percent.mt"]] <- PercentageFeatureSet(kid100, pattern = "^MT-")
VlnPlot(kid100, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid100 <- subset(kid100, subset = nFeature_RNA > 500 & nFeature_RNA < 3800 & percent.mt < 10)

K101.data <- Read10X(data.dir = "RCC101/")
kid101 <- CreateSeuratObject(counts = K101.data, project = "mRCC101", min.cells = 10, min.features = 500)
kid101[["percent.mt"]] <- PercentageFeatureSet(kid101, pattern = "^MT-")
VlnPlot(kid101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid101 <- subset(kid101, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K103.data <- Read10X(data.dir = "RCC103/")
kid103 <- CreateSeuratObject(counts = K103.data, project = "mRCC103", min.cells = 10, min.features = 500)
kid103[["percent.mt"]] <- PercentageFeatureSet(kid103, pattern = "^MT-")
VlnPlot(kid103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid103 <- subset(kid103, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K104.data <- Read10X(data.dir = "RCC104/")
kid104 <- CreateSeuratObject(counts = K104.data, project = "mRCC104", min.cells = 10, min.features = 500)
kid104[["percent.mt"]] <- PercentageFeatureSet(kid104, pattern = "^MT-")
VlnPlot(kid104, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid104 <- subset(kid104, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K106.data <- Read10X(data.dir = "RCC106/")
kid106 <- CreateSeuratObject(counts = K106.data, project = "mRCC106", min.cells = 10, min.features = 500)
kid106[["percent.mt"]] <- PercentageFeatureSet(kid106, pattern = "^MT-")
VlnPlot(kid106, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid106 <- subset(kid106, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K112.data <- Read10X(data.dir = "RCC112/")
kid112 <- CreateSeuratObject(counts = K112.data, project = "mRCC112", min.cells = 10, min.features = 500)
kid112[["percent.mt"]] <- PercentageFeatureSet(kid112, pattern = "^MT-")
VlnPlot(kid112, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid112 <- subset(kid112, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K113.data <- Read10X(data.dir = "RCC113/")
kid113 <- CreateSeuratObject(counts = K113.data, project = "mRCC113", min.cells = 9, min.features = 500)
kid113[["percent.mt"]] <- PercentageFeatureSet(kid113, pattern = "^MT-")
VlnPlot(kid113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid113 <- subset(kid113, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K114.data <- Read10X(data.dir = "RCC114/")
kid114 <- CreateSeuratObject(counts = K114.data, project = "mRCC114", min.cells = 7, min.features = 500)
kid114[["percent.mt"]] <- PercentageFeatureSet(kid114, pattern = "^MT-")
VlnPlot(kid114, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid114 <- subset(kid114, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 10)

K115.data <- Read10X(data.dir = "RCC115/")
kid115 <- CreateSeuratObject(counts = K115.data, project = "mRCC115", min.cells = 9, min.features = 500)
kid115[["percent.mt"]] <- PercentageFeatureSet(kid115, pattern = "^MT-")
VlnPlot(kid115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid115 <- subset(kid115, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K116.data <- Read10X(data.dir = "RCC116/")
kid116 <- CreateSeuratObject(counts = K116.data, project = "mRCC116", min.cells = 6, min.features = 500)
kid116[["percent.mt"]] <- PercentageFeatureSet(kid116, pattern = "^MT-")
VlnPlot(kid116, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid116 <- subset(kid116, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K119.data <- Read10X(data.dir = "RCC119/")
kid119 <- CreateSeuratObject(counts = K119.data, project = "mRCC119", min.cells = 10, min.features = 500)
kid119[["percent.mt"]] <- PercentageFeatureSet(kid119, pattern = "^MT-")
VlnPlot(kid119, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid119 <- subset(kid119, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K120.data <- Read10X(data.dir = "RCC120/")
kid120 <- CreateSeuratObject(counts = K120.data, project = "mRCC120", min.cells = 10, min.features = 500)
kid120[["percent.mt"]] <- PercentageFeatureSet(kid120, pattern = "^MT-")
VlnPlot(kid120, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid120 <- subset(kid120, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

mRCC <- merge(x = kid81, y = list(kid84, kid86, kid87, kid94, kid96, kid99, kid100, kid101, kid103, kid104, kid106, kid112, kid113, kid114, kid115, kid116, kid119, kid120))
mRCC[["percent.mt"]] <- PercentageFeatureSet(mRCC, pattern = "^MT-")
VlnPlot(mRCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size= 0)
plot1 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
mRCC <- NormalizeData(mRCC, normalization.method = "LogNormalize", scale.factor = 10000)
mRCC <- NormalizeData(mRCC)
mRCC <- FindVariableFeatures(mRCC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mRCC), 10)
plot1 <- VariableFeaturePlot(mRCC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
mRCC <- CellCycleScoring(mRCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mRCC <- ScaleData(mRCC, vars.to.regress = c("S.Score", "G2M.Score"))
mRCC <- RunPCA(mRCC, pc.genes = mRCC@var.genes, npcs = 30, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
mRCC <- mRCC %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(mRCC, 'harmony')
harmony_embeddings[1:5, 1:5]
mRCC <- mRCC %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity()
DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
mRCC.markers <- FindAllMarkers(mRCC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(mRCC.markers,sep="\t",file="./scRNA_marker.csv")
new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
#names(new.cluster.ids) <- levels(mRCC)
#mRCC <- RenameIdents(mRCC, new.cluster.ids)
DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
DimPlot(mRCC, reduction = "umap", split.by = "orig.ident",label = TRUE, pt.size = .5, ncol = 5)
##save object
saveRDS(mRCC, file= "./scRNA.RDS")

library(qs)
# 保存单细胞数据为 qs 格式
qsave(mRCC, file = "./scRNA.qs")
# 读取 qs 格式的单细胞数据
scRNA <- qread("./scRNA.qs")
identical(scRNA,mRCC)


markers = c(  'CLDN4','CLDN7','CD24','KRT7','KRT18','KRT19',#Tumor cell
              'EPCAM','GPX3',   #肾小管上皮细胞 #上皮细胞
              'CD3D', 'CD3E', 'IL7R','NKG7', 'KLRD1','CCL5',#NK cell #T cell
              'APOE' , 'C1QA', 'C1QB','CD14',#Macrophage
              'RGS5', 'PLN','HIGD1B',       #Mesangial cell
              'COL1A1', 'PDGFRB', 'ACTA2',#Fibroblast
              'VWF', 'CD34','PECAM1','PLVAP',#Endothelial
              'APOBEC3A',  'FCN1', 'S100A9', 'S100A8', #Monocyte
              'CLEC10A','CD1C','CD1E','FCER1A', #Dendric
              'MS4A1','CD79A', 'JCHAIN','IGHG3', #B cell
              'CPA3','TPSAB1','TPSB2','MS4A2'#Mast cell
)

top_5=markers
#免疫细胞marker表达
library(ggplot2)
scRNA=mRCC
DotPlot(scRNA, features = top_5,cols = c('#ECBC43','red')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(scRNA, reduction = "umap", label = TRUE, pt.size = .5)

table(scRNA@meta.data$seurat_clusters,scRNA@meta.data$celltype)
Idents(scRNA)=scRNA@meta.data$celltype

library(scRNAtoolVis)
scRNAtoolVis::averageHeatmap(object = scRNA,markerGene = unique(top_5))



#细胞比例
library(scRNAtoolVis)
cellRatioPlot(object = scRNA,
              sample.name = "orig.ident",
              celltype.name = "celltype")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
#celltype

cellRatioPlot(object = scRNA,
              sample.name = "tissue_type",
              celltype.name = "celltype")

table(scRNA$seurat_clusters)
### 注释
#NK_T_cells Endothelial_cells Fibroblasts Epithelial_cells Macrophages Monocytes
#B_cells Dendritic_cells Mast_cells Unknown Malignant_cells
celltype <- c('Epithelial_cells' ,         #0 
              'NK_T_cells' ,        #1
              'NK_T_cells' ,         #2
              'NK_T_cells' ,              #3
              'Macrophage' ,              #4
              'Fibroblasts' ,             #5
              'Endothelial_cells' ,             #6
              'Monocytes' ,                #7  
              'NK_T_cells' ,          #8
              'Epithelial_cells' ,             #9 
              'Macrophage' ,          #10
              'Dendritic_cells' ,              #11
              'Epithelial_cells' ,         #12
              'NK_T_cells' ,                 #13
              'Endothelial_cells' ,            #14
              'Epithelial_cells' ,          #15 
              'B_cells' ,           #16 
              'B_cells' ,              #17
              'Mast_cells' ,        #18
              'B_cells' ,                #19 
              'NK_T_cells' )           #20


Idents(scRNA) <- scRNA@meta.data$seurat_clusters
names(celltype) <- levels(scRNA)
scRNA<- RenameIdents(scRNA, celltype)
scRNA@meta.data$celltype <- Idents(scRNA)
#!!!!!!
Idents(scRNA)=scRNA@meta.data$celltype
table(scRNA@meta.data$celltype)


colors=c('#ebc03e','#D0552B','#C0F4CA','#377b4c','#88C7ED',
         '#7bc7cd','#5d84a4','#313c63','#6DA3E5')

p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=5,
             #  split.by = 'tissue_type',
             cols = colors, reduction='umap',pt.size = 1)
p2

saveRDS(scRNA,file ='scRNA_anno.RDS')