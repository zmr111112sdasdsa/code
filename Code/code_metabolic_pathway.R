#首先下载KEGG分级注释
#构建KEGG代谢通路文件 https://www.kegg.jp/kegg-bin/get_htext?br08901.keg
setwd("~/Single_cell_APEX1")
scRNA=readRDS('./scRNA-seqGSE207493/scRNA_APEX_divi.RDS')
set.seed(123)
a=sample(1:ncol(scRNA),20000)
scRNA=scRNA[,a]
Idents(scRNA)
table(scRNA$gene_group)

setwd("./Code/代谢活性分析/")
# library(rjson)
# KEGG_pathway <- rjson::fromJSON(file = "./br08901.json",simplify=F)
# KEGG_pathway <- KEGG_pathway[[2]]#结果在2中
# KEGG_pathway <- KEGG_pathway[[1]]#只取代谢
# KEGG_pathway <- KEGG_pathway[[2]]
# 
# 
# cate1 <- c()
# for (i in 1:length(KEGG_pathway)) {
#   length_names <- length(KEGG_pathway[[i]][[2]])
#   cate1 <- append(cate1, rep(KEGG_pathway[[i]][[1]], length_names))
#   
# }
# 
# cate2 <- c()
# for (i in 1:length(KEGG_pathway)) {
#   
#   cate2 <- append(cate2, KEGG_pathway[[i]][[2]])
#   
# }
# 
# cate3 <- c()
# for (i in 1:length(cate2)) {
#   
#   cate3 <- append(cate3, as.character(cate2[[i]]))
#   
# }
# 
# #去除数字
# library(stringr)
# vec = cate3
# vec <- str_replace(vec, "[0-9]+", "")
# vec <- str_replace(vec, "  ", "")
# 
# KEGG_pathways_Metabolism_anno <- data.frame(cat1  = cate1,
#                                        cat2  = cate3,
#                                        vec  = vec)
# 
# rownames(KEGG_pathways_Metabolism_anno) <- KEGG_pathways_Metabolism_anno$vec

#================================================================================
#================================================================================
#================================================================================
library(stringr)
library(reshape2)
library(scales)
library(scater)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)

#1. Loading the data
library(fgsea)#读入KEGG gmt
pathway_file <- "./KEGG注释文件及通路gmt文件/KEGG_metabolism_nc.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)

#data1: your celltype or sample group
all_cell_types <- as.vector(scRNA$gene_group)
cell_types <- unique(all_cell_types)
#data2 metabolism pathway gene set
pathway_file <- "./KEGG注释文件及通路gmt文件/KEGG_metabolism_nc.gmt"
pathways <- gmtPathways(pathway_file)


Pathway_act_Score <- function(data,
                              pathways,
                              assay,
                              filterGene=F,
                              Mean_cut=NULL,
                              percent_cut=NULL,
                              all_cell_types,
                              cell_types
){
  
  
  DefaultAssay(data) <- assay
  norm_tpm <- GetAssayData(data, layer = "data")
  norm_tpm <- as.matrix(norm_tpm)
  norm_tpm <- as.data.frame(norm_tpm)
  
  ###filter genes
  if(filterGene==F){
    
    norm_tpm = norm_tpm
    
  }else{
    
    pctExpr <- function(mat){
      nexpr <- apply(mat, 1, function(x) sum(x > 0))
      nexpr / ncol(mat)
    }
    
    norm_tpm$percent_exp = pctExpr(norm_tpm)
    norm_tpm$mean_exp = rowMeans(norm_tpm[,-ncol(norm_tpm)])
    
    norm_tpm <- norm_tpm[which(norm_tpm$mean_exp >Mean_cut & norm_tpm$percent_exp >percent_cut),]
    norm_tpm <- norm_tpm[,-ncol(norm_tpm)]
    norm_tpm <- norm_tpm[,-ncol(norm_tpm)]
  }
  
  
  ##calculate how many pathways of one gene involved.
  num_of_pathways <- function (pathway_file,overlapgenes){
    pathway_names <- names(pathway_file)
    filter_pathways <- list()
    for (p in pathway_names){
      genes <- pathway_file[[p]]
      common_genes <- intersect(genes,overlapgenes)
      if(length(common_genes>=5)){
        filter_pathways[[p]] <- common_genes
      }
    }
    
    all_genes <- unique(as.vector(unlist(filter_pathways)))
    gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
    for(p in pathway_names){
      for(g in filter_pathways[[p]]){
        gene_times[g,"num"] = gene_times[g,"num"]+1
      }
    }
    gene_times
  } 
  
  
  #some genes occur in multiple pathways.
  gene_pathway_number <- num_of_pathways(pathways,rownames(norm_tpm))
  
  ##Calculate the pathway activities
  #mean ratio of genes in each pathway for each cell type
  pathway_names <- names(pathways)
  mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
  mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
  ###calculate the pvalues using shuffle method
  pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))
  
  
  for(p in pathway_names){
    genes <- pathways[[p]]
    genes_comm <- intersect(genes, rownames(norm_tpm))
    if(length(genes_comm) < 5) next
    
    pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
    pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
    
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    
    #remove genes which are zeros in any celltype to avoid extreme ratio value
    keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]
    
    if(length(keep)<3) next
    
    #using the loweset value to replace zeros for avoiding extreme ratio value
    pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
    pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))
    
    
    pathway_number_weight = 1 / gene_pathway_number[keep,]
    #
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    #exclude the extreme ratios
    col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
    col_q1 <- col_quantile["25%",]
    col_q3 <- col_quantile["75%",]
    col_upper <- col_q3 * 3
    col_lower <- col_q1 / 3
    outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
    
    if(sum(!outliers) < 3) next
    
    keep <- names(outliers)[!outliers]
    pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
    pathway_number_weight = 1 / gene_pathway_number[keep,]
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
    mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
    mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
    
    ##shuffle 5000 times:  
    ##define the functions 
    group_mean <- function(x){
      sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
    }
    column_weigth_mean <- function(x){
      apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
    }
    #####  
    times <- 1:5000
    weight_values <- pathway_number_weight/sum(pathway_number_weight)
    shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
    names(shuffle_cell_types_list) <- times
    mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
    ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
    mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
    
    shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
    rownames(shuffle_results) <- times
    colnames(shuffle_results) <- cell_types
    for(c in cell_types){
      if(is.na(mean_expression_shuffle[p,c])) next
      if(mean_expression_shuffle[p,c]>1){
        pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
      }else if(mean_expression_shuffle[p,c]<1){
        pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
      }
      if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
      pvalues_mat[p,c] <- pval
    }
  }
  
  
  analysis_res <- list(mean_expression_shuffle,
                       mean_expression_noshuffle,
                       pvalues_mat)
  
  names(analysis_res) <- c("mean_expression_shuffle",
                           "mean_expression_noshuffle",
                           "pvalues_mat")
  return(analysis_res)
  
}


#================================================================================
#analysis------------------------------------------------------------------------
scRNA$newcelltype <- paste0(scRNA$orig.ident,"_", scRNA$gene_group)
#data1: your celltype or sample group
all_cell_types <- as.vector(scRNA$newcelltype)
cell_types <- unique(all_cell_types)
#data2 metabolism pathway gene set
pathway_file <- "./KEGG注释文件及通路gmt文件/KEGG_metabolism_nc.gmt"
pathways <- gmtPathways(pathway_file)

metabolism_activaty <- Pathway_act_Score(scRNA,
                                         pathways=pathways,
                                         assay = "RNA",
                                         filterGene=T,
                                         Mean_cut = 0.001,
                                         percent_cut =0.1,
                                         all_cell_types = all_cell_types,
                                         cell_types = cell_types)

saveRDS(metabolism_activaty,'./metabolism_activaty.RDS')
#plot1----------------------------------------------------------------------------
#通路活性数据就是mean_expression_shuffle
pt <- metabolism_activaty[[1]]
#去除所有列都是NA的数据，为什么会有NA，原因是我们前面筛选表达矩阵了，那么
#有些通路的基因在我们的数据中不存在，所以评分就会是NA呀！我们去除这样的数据
all_NA <- rowAlls(is.na(pt))
pt <- pt[!all_NA,]


#然后就可以进行可视化了，首先使用ggplot进行作图
dat <- pt
#以下内容是原作者提供的数据排序，因为我们要添加注释，所以排序这里对我们没用
#不过感兴趣的还是可以学习以下这个代码，应用在其他的内容中
# sort_row <- c()
# sort_column <- c()
# 
# for(i in colnames(dat)){
#   select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
#   tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
#   sort_row <- c(sort_row,tmp)
# }
# sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
# sort_column <- names(sort_column)
# dat[is.na(dat)] <- 1
# dat <- as.data.frame(dat)

#读入KEGG annotation文件，挑选代谢通路的注释
KEGG_anno <- read.csv("KEGG注释文件及通路gmt文件/kegg_annotation.csv", header = T)
rownames(KEGG_anno) <- KEGG_anno$directory3
Anno_pathway <- KEGG_anno[rownames(dat),]
identical(rownames(dat),rownames(Anno_pathway))
dat=as.data.frame(dat)
dat$pathway_anno <- Anno_pathway$directory2


#宽数据转化为长数据
dat$pathway <- rownames(dat)
melt_dat = melt(dat,
                id.vars = c("pathway","pathway_anno"),
                measure.vars = 1:(length(dat)-2),
                variable.name = c('celltype'),
                value.name = 'score')
#按照分组排个序
melt_dat <- melt_dat %>%arrange(pathway_anno, desc(score))

melt_dat$celltype <- factor(melt_dat$celltype, levels = c("APEX1+Epi","APEX1-Epi","NK_T_cells",
                                                          "Fibroblasts","Dendritic_cells",
                                                          "Macrophage","Endothelial_cells",
                                                          "B_cells","Monocytes","Mast_cells"))
melt_dat$pathway <- factor(melt_dat$pathway, levels = unique(melt_dat$pathway))

p = ggplot(melt_dat,aes(x=celltype,y=pathway, fill=score))+
  geom_tile(color='black')+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,color = "black", size=8),
        axis.text.y = element_text(size=8,vjust=1,color = "black"),
        legend.margin = margin(-0.2,-0.2,0,0,'cm'))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))+
  scale_y_discrete(position = 'right',expand = c(0,0)) + 
  scale_x_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=colorRampPalette(c("#1AA3FF","white","#FF6B67"))(100),
                       guide = guide_colorbar(ticks.colour = "black",
                                              frame.colour = "black"),
                       name = "Pathway activity score",
                       values = c(0,0.15,0.3,0.5,1))

p
