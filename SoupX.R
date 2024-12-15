###############################################################
###################SoupX去除环境中的RNA污染####################
###############################################################
#install.packages('devtools')
#devtools::install_github('JiekaiLab/dior')
#install.packages('SoupX')
rm(list=ls())
#https://github.com/constantAmateur/SoupX
library(SoupX)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})
library(dior)


setwd("H:\\Matrix")
#devtools::instsc_github("constantAmateur/SoupX",ref='devel')
library(SoupX)
library(Seurat)
library(DropletUtils)
#全矩阵，未经过QC的数据
cow_row <-readRDS("./C4-4-TP/Cow_placenta.rds")
#经过QC过的数据
cow_filter<-readRDS("./C4-4-TP/Cow_placenta_raw.RDS")
tod<-cow_row@assays$RNA@counts  #全矩阵
toc<-cow_filter@assays$RNA@counts #经过QC过的矩阵
#保证基因名一致
tod <- tod[rownames(toc),]
##SoupX帮助文档建议提供分析矩阵的聚类亚群分组，因此这里利用分析矩阵做一个简单聚类
all <- toc
all <- CreateSeuratObject(all)
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor =1e6)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)
all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:20)
#提取聚类后的meta.data信息
matx <- all@meta.data
#第二种方法计算
sc = SoupChannel(tod, toc, calcSoupProfile = FALSE)

library(Matrix)
toc = sc$toc
scNoDrops = SoupChannel(tod, toc, calcSoupProfile = FALSE)
# Calculate soup profile
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops, soupProf)

sc<-scNoDrops
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = setContaminationFraction(sc, 0.2)
sc = autoEstCont(sc)
#校正矩阵
out = adjustCounts(sc)
#保存校正后的矩阵，输出为10X格式
DropletUtils:::write10xCounts("./",out,version="3")
#重新读取数据，再降维分群看看数据有什么区别~
srat = CreateSeuratObject(out)
table(srat$orig.ident)
library(harmony)

#####
cow <- NormalizeData(object = srat, scale.factor = 1e6)
cow <- FindVariableFeatures(object = cow, nfeatures = 2000, verbose = FALSE)
cow <- ScaleData(object = cow)
dim(GetAssayData(cow, assay = "RNA", slot = "scale.data"))
cow <- RunPCA(object = cow, assay = "RNA")
pdf(file="cow_FeaturePlot.pdf",width=8,height=7)
FeaturePlot(cow, reduction = "pca", c("nCount_RNA", "nFeature_RNA"))
dev.off()
ElbowPlot(cow, reduction = "pca")
n.pcs = 20
# CLustering and UMAP dimesntion reduction 
cow <- FindNeighbors(object = cow, assay = "RNA", reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
cow <- FindClusters(object = cow, resolution =0.5)
table(Idents(cow))
# To visualize
cow <- RunUMAP(object = cow, assay = "RNA", reduction = "pca", dims = 1:n.pcs)
DimPlot(cow)
