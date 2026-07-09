
gc()
scRNA=readRDS('./scRNA.RDS')
scRNA=subset(scRNA,tissue_type =='Tumor')  ## 一般要运行这行，此处细胞少不运行
scRNA_epi=subset(scRNA,celltype =='Epithelial_cells')
scRNA_epi$APEX1_group=ifelse(as.numeric(scRNA_epi@assays$RNA@counts['APEX1',])>0,'APEX1+Epi','ALG-Epi')

Idents(scRNA_epi)=scRNA_epi$APEX1_group
df=FindAllMarkers(scRNA_epi,only.pos = T,logfc.threshold =0.5 )

genes <- !grepl(pattern = "^RP[L|S]|MT", x = df$gene)
genes

df=df[genes,]

df_pos=df[df$cluster=='APEX1+Epi',]
df_neg=df[df$cluster=='APEX1-Epi',]

write.csv(df_pos,file ='APEX1+markers.csv',quote = F)
write.csv(df_neg,file ='APEX1-markers.csv',quote = F)


## bulk评估浸润---------------------------------------
load('./bulk/KIRC_tpm.Rdata')


## 肿瘤和正常
group=sapply(strsplit(colnames(exprSet_tcga_mRNA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))

exprSet=exprSet_tcga_mRNA[,group_list=='tumor']

max(exprSet)
## 不需Log
#exprSet=log2(exprSet+1)

gene_set=read.csv('./APEX1+markers.csv')
#gene_set=read.csv('./APEX1-markers.csv')

gene_set=gene_set$gene

gene_set=list(gene_set)

names(gene_set)='APEX1+Epi'
#names(gene_set)='APEX1-Epi'

library(genefilter)
library(GSVA)
library(Biobase)

gsva_matrix<- gsva(as.matrix(exprSet), gene_set,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

gsva_matrix=as.data.frame(gsva_matrix)


data=as.data.frame(t(gsva_matrix))
data$ID=rownames(data)


### APEX1阳性的上皮与生存预后 -------------------------------------------------------
library(data.table)
# 读取DFS数据
suv=fread('./bulk/KM_Plot _Disease_Free__Survival_(months).txt',data.table = F)
colnames(suv)=c('project','ID','fustat','futime')
rownames(suv)=suv$ID
cli=dplyr::select(suv,'futime','fustat')
# 去掉尾巴
cli=tidyr::separate(cli,col = 'fustat',sep = ':',into = 'fustat')
cli$fustat=as.integer(cli$fustat)

colnames(cli)=c("futime", "fustat")
cli=na.omit(cli)

# 去掉01a尾巴
data$ID=stringr::str_sub(data$ID,1,12)
# 可能有重复，需去重
data=data[!duplicated(data$ID),]
# 重新设置行名
rownames(data)=NULL
rownames(data)=data$ID


# 取交集
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,]

out=cbind(cli,data)
out=cbind(id=row.names(out),out)

# 输出
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)
#write.table(out,file="expTime2.txt",sep="\t",row.names=F,quote=F)


### 生存分析-------

library(survminer)
library(survival)
library(dplyr)
library(caret)


rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     


library(survival)
library(survminer)


i='APEX1+Epi'
#i='APEX1-Epi'


## 选择最优cutoff
res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(rt[,i]<= cutoff, "Low", "High")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
## 构建生存函数
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
set.seed(1234)
bioCol=col_vector[sample(1:length(col_vector),2)]
p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title=i,
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.88),
             font.legend=12,
             xlab="Time(Months)",
             palette = bioCol,
             #surv.median.line = "hv",
             risk.table=F,
             cumevents=F,
             risk.table.height=.25)
p
