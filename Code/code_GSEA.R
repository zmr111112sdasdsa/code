## gsea富集-------------------------
gc()
setwd('e:/work/IMcompre2/')
load('scRNA_T_receptor.Rdata')
library(Seurat)

## findmarker寻找差异基因
df=FindMarkers(scRNA_T,ident.1 = 'Receptor_high',ident.2 = 'Receptor_low',logfc.threshold = 0)


gc()
df$mRNAs=rownames(df)

## 保存
save(df,file ='deg.Rdata')

#4.制作genelist
gene <- df$mRNAs
## 转换
library(clusterProfiler)
gene = bitr(gene, fromType="SYMBOL", 
            toType="ENTREZID",
            OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,
                        SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(logFC=df$avg_log2FC,
                      SYMBOL = df$mRNAs)
gene_df <- merge(gene_df,
                 gene,by="SYMBOL")

## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, 
                decreasing = TRUE)

#5.运行GSEA分析
library(clusterProfiler)
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("c5.go.bp.v2023.1.Hs.entrez.gmt")
# 需要网络

y <- GSEA(geneList,
          TERM2GENE =hallmarks,pvalueCutoff = 0.25,seed = 123)


y=setReadable(y,OrgDb="org.Hs.eg.db",keyType = 'ENTREZID')


rt=y@result

## 显著的通路
rt=rt[rt$p.adjust<0.25,]

## NES>0代表正相关
rt1=rt[rt$NES>0,]
rt1=rt1[order(rt1$NES,decreasing = T),]
rt1=rt1[1:10,]
rt2=rt[rt$NES<0,]
rt2=rt2[order(rt2$NES),]
rt2=rt2[1:10,]

# 重新组合用来画图
rt=rbind(rt1,rt2)
library(ggplot2)
library(tidyverse)
rt$group <- ''
rt$group[which(rt$NES>0)]='up'
rt$group[which(rt$NES <0)]='down'
colnames(rt)
ggplot(rt,aes(reorder(Description,NES),NES,fill=group))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'none')+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=18.5))+
  geom_text(data = rt[which(rt$NES>0),],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=4)+
  geom_text(data = rt[which(rt$NES<0),],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=4)+
  geom_text(data = rt[which(rt$NES>0),],aes(label=p.adjust),
            hjust=-0.1, size=4, color='red')+
  geom_text(data = rt[which(rt$NES<0),],aes(label=p.adjust),
            hjust=1.1, size=4, color="red")+
  scale_fill_manual(values = c("#1084A4",
                               "#8D4873"))+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  labs(x='', y='NES')


