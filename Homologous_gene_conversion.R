#物种间同源转换脚本，只需要替换物种名即可
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})
#读入RDS文件
cow<-readRDS("")

DefaultAssay(cow)<-'RNA'
dim(cow)
#提取seurat的count矩阵
count <- cow@assays$RNA@counts
dim(count)
#转换为human的同源基因
common_one2one<-read.csv("/home/11-placenta_RNA/start_csv/ten_common_genes.csv")
head(common_one2one)
one2one<-common_one2one[, c("cow_gene", "human_gene")]
#取表达矩阵与one2one物种相同的基因
data2 <- count[which(row.names(count)%in% common_one2one$cow_gene),]
data2[1:4,1:4]
dim(data2)
#使用 one2one 的 "cow_gene" 列中的基因来匹配 data2 数据框的行名
datam2 <- one2one[one2one[,"cow_gene"] %in% rownames(data2),] 
dim(datam2)
####去重
datam <- subset(datam2, !duplicated(datam2$cow_gene))
datam <- subset(datam2, !duplicated(datam2$human_gene))
dim(datam)
####
data2 <- as.data.frame(data2)
dim(data2)
#添加一列
data2[,ncol(data2)+1] =rownames(data2)
###
colnames(data2)[ncol(data2)] ="cow_gene"
data1=merge(datam, data2, by='cow_gene', all=F)
data1[1:4,1:4]
rownames(data1) = data1$human_gene
data1[1:4,1:4]
#删除前两列
data1=data1[,3:ncol(data1)]
data1[1:4,1:4]
dim(data1)
saveRDS(data1,"cow_mat.RDS")
############利用新的表达矩阵创建seurat#############
cellmetadata <-cow@meta.data
saveRDS(cellmetadata,"cow_meta.RDS")
asthma=CreateSeuratObject(counts = data1,project = "cow")
asthma=AddMetaData(asthma,metadata = cellmetadata)
head(asthma)
asthma = asthma %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20)
  
DimPlot(asthma,group.by="subclass2.1",label=T)
asthma$orig.ident<-NULL
asthma$nCount_RNA<-NULL
asthma$nFeature_RNA<-NULL
dim(asthma)
saveRDS(asthma,"cow_cross_TenSpecies.RDS")
marmoset %>%
                                as_tibble() %>%
                                filter(gene %in% orthologous_genes$marmoset_gene) %>%
                                mutate(species = "marmoset") %>%
                                left_join(
                                    orthologous_genes %>% as_tibble() %>%
                                        select(human_gene, marmoset_gene) %>%
                                        set_names("new_gene", "gene"),
                                    by = "gene"
                                ) %>%
                                select(-gene) %>%
                                rename("gene" = "new_gene")
