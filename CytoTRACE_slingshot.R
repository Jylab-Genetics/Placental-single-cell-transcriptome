
setwd("......") #设置路径

library(tradeSeq)
library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(gtools)
library(furrr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(monocle) # version 2.20
library(slingshot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(psych)
library(stats)
#library(cartography)
library(colorRamps)
library(scales)
library(colorspace)
#library(rcartocolor)
args <- commandArgs(trailingOnly=TRUE)
set.seed(42)
## Read RDS
## ---------------------------------- #
R <- readRDS("/home/pig/PIG_TR.rds")

#================================================================================
#1、CytoTRACE分析推断细胞分化起点
#================================================================================
#CytoTRACE分析,输入用单细胞count矩阵，
exp1 <- as.matrix(R @assays$RNA@counts)
exp1 <- exp1[apply(exp1 > 0,1,sum) >= 5,]#筛选一下，表达量太低的也没什么意义，也是缩小数据
#分析很简单，CytoTRACE函数一句话,CytoTRACE函数可多线程
#如果你的数据比较大，那么分析很慢，建议开多线程
results <- CytoTRACE(exp1,ncores = 1)

#meta数据,用于可视化
phenot <- R $cell_type
phenot <- as.character(phenot)
names(phenot) <- rownames(R@meta.data)
emb <- R@reductions[["umap"]]@cell.embeddings
#plotCytoTRACE函数可视化结果，outputDir设置为当然文件夹，可视化结果就会保存在相应文件夹
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = './')
plotCytoGenes(results, numOfGenes = 30, outputDir = './')#可视化基因

#它也会返回一个每个细胞评分的data，我们可以将这个数据读入作图
cyto <- read.table("CytoTRACE_plot_table.txt");head(cyto)
cutoff <- quantile(cyto$CytoTRACE, 0.75)
cyto$diff <- "low"
cyto[cyto$CytoTRACE > cutoff, ]$diff <- "high"
ggplot(cyto, aes(Component1, Component2, color = diff)) +
  geom_point(size = 1.5, alpha = 1.0)




R$cell_type <- Idents(R) %>% as.character()
R_ALL <- R

## 设置Colors，有多少cluster就设置多少颜色
## ---------------------------------- #
real_colors <- c("UNC1" = "#33A02C",
                 "UNC2" = "#B2DF8A",
                 "UNC3" = "#55A1B1",
                 "UNC4" = "#8DD3C7",
                 "UNC5" = "#A6761D",
                 "UNC6" = "#E6AB02") 
umap_colors <- c("UNC1" = "#DC050C",
                 "UNC2" = "#FB8072",
                 "UNC3" = "#1965B0",
                 "UNC4" = "#7BAFDE",
                 "UNC5" = "#882E72",
                 "UNC6" = "#B17BA6")


## Subset 细胞亚群
## ---------------------------------- #
#R <- subset(R, cell_type %in% include_celltype )
root_cell <- "UNC3"
#Idents(R) <- R$cell_type %>% as.character()


## colors
## ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = 8, name = "Dark2")

C <- palette_celltype[as.factor(R$cell_type)]
names(C) <- as.factor(R$cell_type)

color <- unique(C)
names(color) <- unique(names(C))

## run Slingshots 
## ---------------------------------- #
start.clus <- root_cell
reduction = 'umap' #这里需要根据自己的数据更改
sds     = slingshot(Embeddings(R, reduction), clusterLabels = Idents(R), start.clus = start.clus )
R@tools[['slingshot']] = SlingshotDataSet(sds)
pseudotime = slingPseudotime(sds)

sds_all = slingshot(Embeddings(R_ALL, reduction), clusterLabels = Idents(R_ALL), start.clus = start.clus )

## Plot slingshot curves#
curves = colnames(pseudotime)

## Plot slingshot curves
## ---------------------------------- #
sds_all$reducedDim %>% rownames() -> cell_id
umap_colors[R_ALL$cell_type] -> cell_color
names(cell_color) <- colnames(R_ALL)
cell_color[cell_id] -> U_C

## Plot slingshot curve II by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds

R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
R_META <- R_ALL[["cell_type"]] %>% dplyr::mutate("cell_id" = rownames(.))
R_UMAP <- Embeddings(object = R_ALL, reduction = "umap") %>% data.frame() %>% 
  dplyr::mutate("cell_id" = rownames(.)) %>% 
  left_join(., R_META, by="cell_id") %>% 
  left_join(., R_TIME, by="cell_id") 

#画一个UMAP
ggplot(R_UMAP, aes(x=UMAP_1, y=UMAP_2, group = eval(parse(text="cell_type")), colour = eval(parse(text="cell_type")))) +
  geom_point(size=0.8, alpha=0.7) + 
  theme_light() + theme_classic() +
  scale_colour_manual(values=umap_colors) + 
  theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none')

#这一步很重要，画一个UMAP,把需要的cluster及谱系标注出来，
#可以修改Lineage这个参数
library(tidyr)
R_UMAP_long <- R_UMAP %>%
  gather(key = "Lineage", value = "Value", Lineage1:Lineage2)
# 绘制图形
ggplot(R_UMAP_long, aes(x = UMAP_1, y = UMAP_2, colour = Value)) +
  geom_point(size = 0.8, alpha = 0.7) +
  theme_light() + theme_classic() +
  scale_colour_viridis_c(na.value = "#D3D3D3", option = "C") +
  theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_alpha(guide = 'none') +
  facet_wrap(~ Lineage)
ggsave(filename = "R_UMAP.svg", width= 5, height = 5)
#make_umap_c(R_UMAP, "Lineage1")
#make_umap_c(R_UMAP, "Lineage2")

#给谱系标记上曲线，这一步的结果很重要
## Plot slingshot curve by cell type
## ---------------------------------- #
pseudotime_orig <- pseudotime
sds_orig <- sds
pdf(file = 'slingshot_curves.separate.colored.by.time.all.umap.NEW.pdf', width = 7.67, height=4.23)
par(mfrow = c(1, 2))
for ( c_num in seq(1, length(curves))) {
  
  sds <- sds[,curves[c_num]]
  
  sds$reducedDim %>% rownames() -> cell_id
  sds_all$reducedDim %>% rownames() %>% as.data.frame() -> all_cell_id
  colnames(all_cell_id) <- "cell_id"
  
  pseudotime = slingPseudotime(sds)
  colors = palette[cut(pseudotime[,1], breaks = 100)]
  names(colors) <- cell_id
  colors %>% as.data.frame() %>% rownames_to_column() -> colors_df  
  colnames(colors_df) <- c("cell_id", "color")
  
  all_cell_id %>% 
    left_join(., colors_df, by="cell_id") %>% 
    mutate(color = ifelse(is.na(.[["color"]])==TRUE, "#D3D3D3", color)) -> new_colors 
  
  print(plot(sds_all$reducedDim, col = new_colors$color, pch = 16, cex = 0.5, main = curves[c_num] ) +
          lines(SlingshotDataSet(sds), linInd = c_num, lwd = 2, col = 'black'))
  
  sds <- sds_orig
  pseudotime <- pseudotime_orig
  
}
dev.off()
## Add pseudotimes to meta data of R object 
## ---------------------------------- #
## > R[[c("Lineage1","Lineage2")]] %>% as.data.frame() -> b
## > pseudotime %>% as.data.frame() -> a
## > all.equal(a,b)
## [1] TRUE
## ---------------------------------- #
for ( curve in curves ) {
  pseudotime_sub <- pseudotime[colnames(R),curve]
  R <- AddMetaData(object = R,
                   metadata = pseudotime_sub,
                   col.name = curve
  )
}


## Condition density along pseudotime
## ---------------------------------- #
df <- data.frame(R[["cell_type"]], R[[curve]]) 
colnames(df) <- c("cell_type", "Lineage1")
na.omit(df) -> df
ggplot(df, aes(x=Lineage1, fill=cell_type)) +
  geom_density(alpha=0.4) + theme_classic()+
  scale_fill_manual(values=C)
ggsave("cell_density.pdf")
######################################################
######################################################
#以下代码绘制拟时间的热图
L <- R[["Lineage1"]] %>% deframe()
names(L) <- rownames(R[["Lineage1"]])
L[!is.na(L)] %>% names() -> L2_cell
#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
R[, L2_cell] -> new_R 
Idents(new_R) <- new_R$cell_type %>% as.character()

# Transfer into cds
# ---------------------------------- #
cds <- as.CellDataSet(new_R)
# Estimate size factor
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# select superset of feature genes as genes expressed in at least 5% of all the cells.
# ---------------------------------- #
cds <- detectGenes(cds, min_expr = 0.1)
fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
cds_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))

clustering_DEG_genes <- differentialGeneTest(cds[cds_genes,],
                                             fullModelFormulaStr = '~sample',
                                             cores = 20)

cds_ordering_df <- clustering_DEG_genes %>% filter(qval < 0.05 & use_for_ordering == TRUE) %>% arrange(qval)
cds_ordering_df[1:200, ] %>% select(gene_short_name) %>% deframe() %>% as.character() -> cds_ordering_genes

print(head(pData(cds)))

pData(cds)[["Lineage1"]] -> pData(cds)$Pseudotime
diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=1)", 
                                      cores = 10)

diff_test_res %>% filter(qval <0.01) %>% arrange(qval)

sig_gene_names <- row.names(diff_test_res)
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 2,
                        cores = 10,
                        show_rownames = T)

#########################
#########################
#下面需要先加载get_pseudotime_matrix自定义函数，这个函数在节末提供
# make plots
# ---------------------------------- #
DE_N_Set <- c(50) # c(20, 30, 50)

for ( DE_N in DE_N_Set ) { 
  
  diff_test_res %>% filter(qval <0.01) %>% arrange(qval) %>% top_n(., DE_N, wt=-qval) %>% rownames() -> sig_gene_names
  
  
  hm <- get_pseudotime_matrix(cds[sig_gene_names,],  
                              cluster_rows = TRUE,
                              hclust_method = "ward.D",
                              num_clusters = 6,
                              hmcols = NULL,
                              add_annotation_row = NULL,
                              add_annotation_col = NULL,
                              show_rownames = FALSE,
                              use_gene_short_name = TRUE,
                              norm_method = "log",
                              scale_max=3,
                              scale_min=-3,
                              trend_formula = "~sm.ns(Pseudotime, df=1)", 
                              return_heatmap=TRUE,
                              cores=1)
  
  
  bks = c(seq(min(hm), 0, length.out=ceiling(200/2) + 1),
          seq(max(hm)/200, max(hm),length.out=floor(200/2)))
  
  my_color4 = plasma(length(bks))
  my_color5 = colorRampPalette(rev(rcartocolor::carto_pal(7, "Sunset")))(length(bks))
  my_color6 = colorRampPalette(rcartocolor::carto_pal(7, "ag_Sunset"))(length(bks))
  my_color7 = colorRampPalette(rev(rcartocolor::carto_pal(7, "SunsetDark")))(length(bks))
  
  my_color_set <- list(my_color4, my_color5, my_color6, my_color7)
  my_color_name <- c("plasma", "Sunset", "ag_Sunset", "SunsetDark")
  
  
  
  # cluster and re-order rows
  # ---------------------------------- #
  # ALL_HCS <- c( "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  # ---------------------------------- #
  
  ALL_HCS <- c("ward.D")
  
  for ( sub_color in seq(1,length(my_color_set)))  {
    for ( ALL_HC in c(ALL_HCS) ) {  
      print(monocle::plot_pseudotime_heatmap(cds[sig_gene_names,],
                                             ##add_annotation_col = your_column,
                                             cluster_rows = TRUE,
                                             trend_formula = "~sm.ns(Pseudotime, df=1)",
                                             hclust_method = ALL_HC, 
                                             num_clusters = 1,
                                             hmcols = my_color_set[sub_color][[1]],
                                             scale_max = 3, 
                                             scale_min = -3,
                                             cores = 1,
                                             show_rownames = T,
                                             return_heatmap = FALSE))
      dev.off()
    }
    
    
  }
  
}

##以下代码画图
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        ##add_annotation_col = your_column,
                        cluster_rows = TRUE,
                        trend_formula = "~sm.ns(Pseudotime, df=1)",
                        hclust_method = ALL_HC, 
                        num_clusters = 1,
                        hmcols = my_color_set[sub_color][[1]],
                        scale_max = 3, 
                        scale_min = -3,
                        cores = 20,
                        show_rownames = T,
                        return_heatmap = FALSE)


######################
######################
######################
#画单个基因随着伪时间的变化图
lung_genes <- row.names(subset(fData(cds),gene_short_name %in% c("Pdgfra","Itga9","Wnt5a","Ryk")))
lung_genes_subset <- cds[lung_genes,]
plot_genes_in_pseudotime(lung_genes_subset,color_by = "cell_type",ncol = 2)


######################
######################
######################
###get_pseudotime_matrix自定义函数

setwd("./")

library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(gtools)
library(furrr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(monocle) # version 2.20
library(slingshot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(stats)
args <- commandArgs(trailingOnly=TRUE)
set.seed(42)



## this script is from 
## https://github.com/cole-trapnell-lab/monocle-release/blob/master/R/plotting.R
## ------------------------------- ## 



get_pseudotime_matrix <- function(cds_subset, 
                                  
                                  cluster_rows = TRUE,
                                  hclust_method = "ward.D2", 
                                  num_clusters = 6,
                                  
                                  hmcols = NULL, 
                                  
                                  add_annotation_row = NULL,
                                  add_annotation_col = NULL,
                                  show_rownames = FALSE, 
                                  use_gene_short_name = TRUE,
                                  
                                  norm_method = c("log", "vstExprs"), 
                                  scale_max=3, 
                                  scale_min=-3, 
                                  
                                  trend_formula = '~sm.ns(Pseudotime, df=3)',
                                  
                                  return_heatmap=FALSE,
                                  cores=1){
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), max(pData(cds_subset)$Pseudotime),length.out = 100)) 
  
  m <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,  
                       relative_expr = T, new_data = newdata)
  
  
  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]
  
  norm_method <- match.arg(norm_method)
  
  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds_subset, expr_matrix=m)
  }     
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }
  
  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min
  
  heatmap_matrix <- m
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  } 
  
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=cluster_rows, 
                 show_rownames=F, 
                 show_colnames=F, 
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 border_color = NA,
                 color=hmcols)
  
  if(cluster_rows) {
    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  } else {
    annotation_row <- NULL
  }
  
  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])  
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }
  
  if(!is.null(add_annotation_col)) {
    if(nrow(add_annotation_col) != 100) {
      stop('add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!')
    }
    annotation_col <- add_annotation_col
  } else {
    annotation_col <- NA
  }
  
  if (use_gene_short_name == TRUE) {
    if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      
      row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    if(!is.null(annotation_row))
      row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  if(!is.null(annotation_row))
    row.names(annotation_row) <- row_ann_labels
  
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  
  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols = FALSE, 
                     cluster_rows = cluster_rows, 
                     show_rownames=show_rownames, 
                     show_colnames=F, 
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     treeheight_row = 20, 
                     breaks=bks,
                     fontsize = 6,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE,
                     filename=NA
  )
  
  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(heatmap_matrix)
  }
}
