# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat)
packageVersion("Seurat");
library(Matrix); library(stringr); library(readr)
library(here); library(fitdistrplus)
library(dplyr); library(monocle); library(reticulate)
library(ggplot2)
setwd("/work/home/luo_funong/Tan_Guang-hui/single_cell_RNA")

load(file = "cc.genes.rda") 
mito_genes = c('ND1','ND2',
               'COX1','COX2','ATP8','ATP6',
               'COX3','ND3','ND4L','ND4','ND5','ND6','CYTB')
samples = c("GSM6241484_19_5_3-PD", "GSM6241484_19_5_2-PD", "GSM6241484_19_5_1-PD",
            "GSM6241481_15_5_3-PD", "GSM6241481_15_5_2-PD", "GSM6241481_15_5_1-PD")

cc.genes <- NULL

prepare_datasets <- function(sample_name, cc.genes, mito_genes){
  # Read data
  data <- Read10X(data.dir = paste("Rattus_norvegicus/", sample_name, "/filtered_feature_bc_matrix/", sep = ""))
  seurat.object <- CreateSeuratObject(counts = data, min.cells = 1, min.features =200, project = sample_name)
  # Mitrocondria
  kp=mito_genes %in% rownames(seurat.object)
  table(kp)
  mito_genes=mito_genes[kp]
  C<-GetAssayData(object = seurat.object, slot = "counts")
  percent_mito <- Matrix::colSums(C[mito_genes,])/Matrix::colSums(C)*100
  seurat.object <- AddMetaData(seurat.object, percent_mito, col.name = "percent_mito")
  return(seurat.object)
}

day19_3 = prepare_datasets(samples[1], cc.genes, mito_genes)
day19_2 = prepare_datasets(samples[2], cc.genes, mito_genes)
day19_1 = prepare_datasets(samples[3], cc.genes, mito_genes)
day15_3 = prepare_datasets(samples[4], cc.genes, mito_genes)
day15_2 = prepare_datasets(samples[5], cc.genes, mito_genes)
day15_1 = prepare_datasets(samples[6], cc.genes, mito_genes)
dim(day19_3)
dim(day19_2)
dim(day19_1)
dim(day15_3)
dim(day15_2)
dim(day15_1)
#保存原始RDS
# List of objects to save
objects_list <- list(day19_3, day19_2, day19_1,day15_3,day15_2,day15_1)
# Names for the RDS files
file_names_1 <- c("Rattus_norvegicus_E19_3_row.rds", "Rattus_norvegicus_E19_2_row.rds", 
                  "Rattus_norvegicus_E19_1_row.rds", "Rattus_norvegicus_E15_3_row.rds", "Rattus_norvegicus_E15_2_row.rds",
                  "Rattus_norvegicus_E15_1_row.rds")
# Loop through the objects and save them as RDS files
for (i in 1:length(objects_list)) {
  saveRDS(objects_list[[i]], file = file_names_1[i])
}


# save.image("robjs/all.objs.RData")
# load("robjs/all.objs.RData")

Rattus_norvegicus_placenta = merge(day19_3, y = c(day19_2, day19_1, day15_3, day15_2, day15_1), add.cell.ids = samples, project = "Rattus_norvegicus_placentaEmbryo")
dim(Rattus_norvegicus_placenta)
table(Rattus_norvegicus_placenta$orig.ident)

# Analysis starts here
vln_plot <-VlnPlot(object = Rattus_norvegicus_placenta, features = c("nCount_RNA", "nFeature_RNA", "percent_mito"), pt.size = 0.01, group.by = "orig.ident")
# 保存为PDF
ggsave("Rattus_norvegicus_vln_plot.pdf", plot = vln_plot, device = "pdf")
FeatureScatter(object = Rattus_norvegicus_placenta, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", use.raw=T)
FeatureScatter(object = Rattus_norvegicus_placenta, feature1 = "nFeature_RNA", feature2 = "percent_mito", do.return=TRUE) + geom_hline(yintercept = 20) + geom_vline(xintercept = 200)

# Filter out cells with few reads and few genes.
Rattus_norvegicus_placenta <- subset(Rattus_norvegicus_placenta, subset = nFeature_RNA >= 500 & percent_mito <= 20)
dim(Rattus_norvegicus_placenta)

# Assign cell cycle score
Rattus_norvegicus_placenta <- CellCycleScoring(object = Rattus_norvegicus_placenta, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Normalise, scale, and run PCA
Rattus_norvegicus_placenta <- NormalizeData(object = Rattus_norvegicus_placenta, scale.factor = 1e6)

Rattus_norvegicus_placenta <- FindVariableFeatures(object = Rattus_norvegicus_placenta, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Rattus_norvegicus_placenta <- ScaleData(object = Rattus_norvegicus_placenta)
dim(GetAssayData(Rattus_norvegicus_placenta, assay = "RNA", slot = "scale.data"))
Rattus_norvegicus_placenta <- RunPCA(object = Rattus_norvegicus_placenta, assay = "RNA")
FeaturePlot(Rattus_norvegicus_placenta, reduction = "pca", c("nCount_RNA", "nFeature_RNA"))
ElbowPlot<-ElbowPlot(Rattus_norvegicus_placenta, reduction = "pca")
ggsave("Rattus_norvegicus_ElbowPlot.pdf", plot = ElbowPlot, device = "pdf")
#n.pcs = 20

# CLustering and UMAP dimesntion reduction 
#Rattus_norvegicus_placenta <- FindNeighbors(object = Rattus_norvegicus_placenta, assay = "RNA", reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
#Rattus_norvegicus_placenta <- FindClusters(object = Rattus_norvegicus_placenta, resolution = 0.5)
#table(Idents(Rattus_norvegicus_placenta))

# To visualize
#Rattus_norvegicus_placenta <- RunUMAP(object = Rattus_norvegicus_placenta, assay = "RNA", reduction = "pca", dims = 1:n.pcs)
#DimPlot(Rattus_norvegicus_placenta, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.5")
#DimPlot(Rattus_norvegicus_placenta, reduction = "umap", group.by = "orig.ident")

# save(Rattus_norvegicus_placenta, file="robjs/Rattus_norvegicus_placenta_raw.Robj")
# load(here("robjs", "Rattus_norvegicus_placenta_raw.Robj"))
