

library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(progeny)
library(dorothea)
library(viper)
library(decoupleR)

#######一、转录因子活性#########
#BiocManager::install("OmnipathR")
#####单细胞转录因子活性评分、通路活性评分
scRNA <- readRDS("~/Single_cell_APEX1/scRNA-seqGSE207493/scRNA_APEX_divi.RDS")

#1 加载数据集合----## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
##如果是小鼠，就用
##dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))

#1 过滤数据集的置信度----## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%  
  dplyr::filter(confidence %in% c("A","B","C"))
#类似结果，作者现在更推荐collectri的结果！
net <- decoupleR::get_collectri(split_complexes = T,organism='human' )
head(net)
regulon

set.seed(123)
scRNA=subset(scRNA,celltype %in% 'Epithelial_cells')
#a=sample(1:ncol(scRNA),20000)
#scRNA=scRNA[,a]
Idents(scRNA)
scRNA$gene_group=factor(scRNA$gene_group,levels = c('APEX1-Epi','APEX1+Epi'))
table(scRNA$gene_group)

#2计算细胞的TF活性----## We compute Viper Scores 
scRNA <- dorothea::run_viper(scRNA, regulon,                  
                             options = list(method = "scale", minsize = 4,                                  
                                            eset.filter = FALSE, cores = 10,                                  
                                            verbose = FALSE))

AA=scRNA@assays$dorothea@data


#3可视化方案1
#3 对dorothea矩阵进行降维聚类分群-----## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object =  scRNA) <- "dorothea"
scRNA <- ScaleData( scRNA)
scRNA <- RunPCA( scRNA, features = rownames( scRNA), verbose = FALSE)
scRNA <- FindNeighbors( scRNA, dims = 1:10, verbose = FALSE)
scRNA <- FindClusters( scRNA, resolution = 0.5, verbose = FALSE)
scRNA <- RunUMAP( scRNA, dims = 1:10, umap.method = "uwot", metric = "cosine")
scRNA.markers <- FindAllMarkers( scRNA, only.pos = TRUE, min.pct = 0.25,
                                 logfc.threshold = 0.25, verbose = FALSE)
top10 <-  scRNA.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC);top10
scRNA@assays$dorothea@data[1:4,1:4]
top10 =top10[top10$avg_log2FC > 0.3,] 
DoHeatmap( scRNA,top10$gene,size=3,slot='scale.data')


####可视化方案三
#4获取Viper得分矩阵，评价细胞群的TF活性----
## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
DefaultAssay(object =  scRNA) <- "dorothea"
scRNA <- ScaleData(scRNA)
viper_scores_df <- GetAssayData( scRNA, slot = "scale.data",
                                 assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t();viper_scores_df[1:4,1:4]

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents( scRNA)),
                            cell_type = as.character(Idents( scRNA)),
                            check.names = F) #也可以使用其他的分类信息
CellsClusters[1:4, ]
## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>%
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters);viper_scores_clusters[1:4,]

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity));summarized_viper_scores

table(summarized_viper_scores$cell_type)

#5细胞群间变化最大的20个TFs进行可视化----
## We select the 20 most variable TFs. (20*9 populations = 180)  9个细胞亚群
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%  mutate(var = var(avg))  %>% # 计算变异度 calculates the variance of the "avg" column. So, for each transcription factor, it computes the variance of its average score.  
  ungroup() %>%
  top_n(90, var) %>%
  distinct(tf);highly_variable_tfs


## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) ;summarized_viper_scores_df



palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                   length.out=ceiling(palette_length/2) + 1),  
               seq(max(summarized_viper_scores_df)/palette_length,  
                   max(summarized_viper_scores_df),            
                   length.out=floor(palette_length/2)))


#summarized_viper_scores_df=summarized_viper_scores_df[c(5,6,15,16,9,4,3,10,7,14,8,1,2,11,12,13),]
#rownames(summarized_viper_scores_df)=c('C1','C2','C3','C4', 'C5','C6','C7','C8', 'C9','C10','C11','C12','C13','C14','C15','C16')

viper_hmap <- pheatmap::pheatmap(t(summarized_viper_scores_df),fontsize=14,    
                                 fontsize_row = 10,               
                                 color=my_color, breaks = my_breaks,   
                                 cluster_cols = F,cluster_rows = T,
                                 main = "DoRothEA (ABC)", angle_col = 45,   
                                 #  treeheight_col = 0,  
                                 border_color = NA)

pdf_file_path <- './Figure2/TF_pathway/Endo_TF.pdf' # Using a cleaner file path
# Create the directory if it doesn't exist
dir.create(dirname(pdf_file_path), showWarnings = FALSE) 
cairo_pdf(pdf_file_path, width = 6, height = 9, family = "sans")
viper_hmap
dev.off()


#可视化方案四
#遗憾的是在这种可视化方案中，效果似乎不太好。
#这也可以理解，毕竟这是转录因子活性分析，看的是转录因子下游基因的表达情况，而不是转录因子本身。所以对转录因子本身进行可视化似乎不是一个好的选择
Idents( scRNA)='RNA'
scRNA@meta.data$GATA1=  scRNA@assays$dorothea@scale.data['GATA1',] 
FeaturePlot( scRNA,features = "GATA1",reduction = "umap" )| DimPlot( scRNA,label=T,group.by = 'cell.type')
FeaturePlot( scRNA,features = "PBX2",reduction = "umap" )| DimPlot( scRNA,label=T,group.by = 'cell.type')

#######二、通路活性#########

# We create a data frame with the specification of the cells that belong to ## each cluster to match with the Progeny scores. 
CellsClusters <- data.frame(Cell = names(Idents( scRNA)), 
                            CellType = as.character(Idents( scRNA)),  
                            stringsAsFactors = FALSE) 
head(CellsClusters)
DimPlot( scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


#2-------## We compute the Progeny activity scores and add them to our Seurat object 
## as a new assay called Progeny. 
scRNA <- progeny( scRNA, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)  #"Human" Mouse
scRNA@assays$progeny 
scRNA@assays$progeny %>%dim()
scRNA@assays$progeny@data[,1:19]
# Assay data with 14 features for 2638 cells # First 10 features: 
# Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
scRNA <- Seurat::ScaleData(scRNA, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(scRNA, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

#三种可视化
library(viridis)
progeny_hmap = pheatmap::pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                                  fontsize_row = 10, 
                                  color=myColor, breaks = progenyBreaks, 
                                  main = "PROGENy (500)", angle_col = 45,
                                  treeheight_col = 0,  border_color = NA)

summarized_progeny_scores_df=summarized_progeny_scores_df[c(5,6,15,16,9,4,3,10,7,14,8,1,2,11,12,13),]
rownames(summarized_progeny_scores_df)=c('C1','C2','C3','C4',
                                         'C5','C6','C7','C8',
                                         'C9','C10','C11','C12',
                                         'C13','C14','C15','C16')
progeny_hmap = pheatmap::pheatmap(t(summarized_progeny_scores_df),
                                  fontsize=12,
                                  fontsize_row = 10, color=myColor,  
                                  cluster_cols = F,cluster_rows = T,
                                  breaks = progenyBreaks, main = "PROGENy",       
                                  angle_col = 45, treeheight_col = 0, border_color = NA)

pdf_file_path <- './Figure2/TF_pathway/Endo_pathway.pdf' # Using a cleaner file path
# Create the directory if it doesn't exist
dir.create(dirname(pdf_file_path), showWarnings = FALSE) 
cairo_pdf(pdf_file_path, width = 6, height = 5, family = "sans")
progeny_hmap
dev.off()

progeny_hmap = pheatmap(t(summarized_progeny_scores_df), 
                        fontsize=12,           
                        fontsize_row = 10, color=turbo(90), #"inferno" "magma" "cividis" "viridis"         
                        #  breaks = progenyBreaks,                  
                        main = "PROGENy",        
                        angle_col =90, treeheight_col = 0, border_color = NA)
#7*4.7


