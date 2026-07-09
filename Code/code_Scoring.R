library(Seurat)
library(SeuratData)
library(UCell)
library(irGSEA)
library(AUCell) ## 版本需要大于1.14

#BiocManager::install('AUCell')

## 可再次读取文件
#scRNA=readRDS('./scRNA_freash_anno_1v1.RDS')

set.seed(123)
a=sample(1:ncol(scRNA),25000)
scRNA=scRNA[,a]
table(scRNA$gene_group)


## 1. AUCell 
cells_rankings <- AUCell_buildRankings(scRNA@assays$RNA@data,  nCores=1, plotStats=TRUE) 

gc()


cells_AUC <- AUCell_calcAUC(geneset, cells_rankings,nCores =1, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- as.numeric(getAUC(cells_AUC)['Ferroptosis', ])


## 2. Ucells and singscore
scRNA <- irGSEA.score(object = scRNA, assay = "RNA",
                      slot = "data", seeds = 123, ncores = 1,msigdb = F,
                      custom = T, geneset = geneset,
                      method = c('UCell','singscore'),
                      kcdf = 'Gaussian')

uc=as.data.frame(scRNA@assays$UCell@counts)
uc=as.data.frame(t(uc))

singscore=as.data.frame(scRNA@assays$singscore@counts)
singscore=as.data.frame(t(singscore))


## 3. GSVA
exp<- as.matrix(scRNA@assays$RNA@data)


library(GSVA)

# 第一步：创建 ssGSEA 方法的参数对象
# 注意：新版 ssgseaParam 不需要 kcdf 参数，且默认处理方式可能已包含你之前的需求
gsva_par <- ssgseaParam(exprData = exp, 
                        geneSets = geneset)

# 第二步：运行 GSVA 分析
matrix <- gsva(gsva_par)

ssgsea=as.data.frame(t(matrix))


# 5. addmodulescore自带函数
scRNA=AddModuleScore(scRNA,features = geneset,name = 'Add')



# 组合
score=data.frame(AUCell=aucs,UCell=uc$Ferroptosis,singscore=singscore$Ferroptosis,ssgsea=ssgsea$Ferroptosis,Add=scRNA$Add1)

# scale 标准化
score<- scale(score)

# 0-1标准化

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

score=apply(score, 2, normalize)
score=as.data.frame(score)
score$Scoring=rowSums(score)
colnames(scRNA@meta.data)

scRNA$Add1=NULL

# 将score添加入metadata，需要充分理解 
scRNA@meta.data=cbind(scRNA@meta.data,score)

library(RColorBrewer) 
library(viridis)
library(wesanderson)

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 12, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:12]


dev.off()
library(ggplot2)
DotPlot(scRNA,features = colnames(score)) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45,  hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")



VlnPlot(scRNA,features = colnames(score))

VlnPlot(scRNA,features = 'SLC7A11')

library(viridis)
FeaturePlot(scRNA,features = 'Scoring',
            cols = magma(10),label = T,reduction = 'umap')
Idents(scRNA)=scRNA$gene_group

## 作显著性分析
rt=scRNA@meta.data
rt=rt[,1:21]
rt=rt[rt$gene_group=='APEX1+Epi'|rt$gene_group=='APEX1-Epi',]

library(ggsci)
library(ggplot2)
library(ggpubr)
ggviolin(
  rt,
  x = "gene_group",
  y = "Scoring",
  color = "black",
  fill = col_vector[2],
  xlab = "",
  ylab = "Scoring"
) +
  stat_compare_means(
    aes(group = gene_group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = F,
    size = 4.5
  ) +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))+scale_fill_jco()

saveRDS(scRNA,'scRNA_APEX1_score.RDS')
